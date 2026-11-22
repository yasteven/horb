#!/usr/bin/env python3
"""Resumable SPARC signed-density fits; invoked by ./run.sh --release do sparc-fit."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import threading
import time

import numpy as np
import scipy
from horb_sparc import (ROOT, DATA, G, STATES, digest, save_json, load_data,
                        select_galaxies, mass_to_light, prepare_galaxy, inventory)
from horb_sparc.gallery import publish, overview


def parser():
    p = argparse.ArgumentParser(description=__doc__, epilog="All source data are bundled. Same command resumes verified completed galaxies. Change --name when changing scientific/plot settings.")
    p.add_argument("--name", default="default", help="run directory name under reports/sparc and plots/sparc")
    p.add_argument("--galaxies", default="all", help="all or comma-separated exact SPARC names")
    p.add_argument("--jobs", type=int, default=1, help="parallel galaxies; 1 is conservative on Jetson")
    p.add_argument("--force", action="store_true", help="rerun selected galaxies even if their results verify")
    p.add_argument("--search", choices=["beam", "anneal", "exhaustive"], default="beam")
    p.add_argument("--budget", type=int, default=2000)
    p.add_argument("--beam-width", type=int, default=6)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--max-components", type=int, default=4)
    p.add_argument("--max-abs-weight", type=float, default=20)
    p.add_argument("--criterion", choices=["chi2", "bic"], default="chi2")
    p.add_argument("--ridge", type=float, default=0)
    p.add_argument("--weights", choices=["signed", "nonnegative"], default="signed")
    p.add_argument("--ml-disk", type=float, default=0.5)
    p.add_argument("--ml-bulge", type=float, default=0.7)
    p.add_argument("--states", default=",".join(STATES))
    scale = p.add_mutually_exclusive_group()
    scale.add_argument("--scale-factors", help="a0/Rdisk values; default: 0.125,0.25,0.5,1,2")
    scale.add_argument("--scales", help="absolute a0 values in kpc, instead of scaling by each disk")
    p.add_argument("--ref-mass", type=float, help="override per-galaxy Rmax*max(Vobs)^2/G normalization in Msun")
    p.add_argument("--n-radial", type=int, default=128)
    p.add_argument("--grid-points", type=int, default=129)
    p.add_argument("--plots", choices=["both", "rotation", "none"], default="both")
    p.add_argument("--density-bins", type=int, default=201)
    p.add_argument("--curve-bin", type=Path, help=argparse.SUPPRESS)
    p.add_argument("--validate-only", action="store_true", help=argparse.SUPPRESS)
    return p


def validate(args, catalog, curves):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}", args.name):
        raise ValueError("--name must be 1-80 letters, digits, underscores or hyphens, starting with a letter/digit")
    names = select_galaxies(args.galaxies, catalog)
    mass_to_light(args.ml_disk, args.ml_bulge)
    states = args.states.split(",")
    if len(states) != len(set(states)) or not set(states) <= set(STATES) or not states:
        raise ValueError(f"--states requires distinct names from {','.join(STATES)}")
    values = [float(x) for x in (args.scales or args.scale_factors or "0.125,0.25,0.5,1,2").split(",")]
    if not values or not np.isfinite(values).all() or min(values) <= 0 or len(values) != len(set(values)):
        raise ValueError("scales/factors must be distinct, positive and finite")
    if (args.jobs < 1 or args.budget < 1 or args.beam_width < 1 or args.seed < 0
            or args.max_components < 1 or not 8 <= args.n_radial <= 512
            or not 2 <= args.grid_points <= 10000 or not 3 <= args.density_bins <= 1001):
        raise ValueError("invalid jobs/search/grid settings; see --help (radial order 8..512, grid 2..10000, density bins 3..1001)")
    if (not np.isfinite([args.max_abs_weight, args.ridge]).all() or args.max_abs_weight <= 0
            or args.ridge < 0 or args.ref_mass is not None and (not np.isfinite(args.ref_mass) or args.ref_mass <= 0)):
        raise ValueError("invalid weight bound, ridge or reference mass")
    count = len(states)*len(values)
    if count > 300:
        raise ValueError("at most 300 basis columns are supported per batch")
    if args.budget < count+1:
        raise ValueError(f"--budget must be >= {count+1} to evaluate the null model and every single orbital")
    if args.search == "exhaustive":
        k = min(args.max_components, count, max(len(curves[n])-1 for n in names))
        needed = sum(math.comb(count, i) for i in range(k+1))
        if args.budget < needed:
            raise ValueError(f"exhaustive search needs --budget >= {needed} for the selected galaxies")
    config = {k: v for k, v in vars(args).items()
              if k not in {"name", "galaxies", "jobs", "force", "curve_bin", "validate_only", "scales", "scale_factors"}}
    config.update(scale_mode="absolute_kpc" if args.scales else "disk_scale_multiple", scale_values=values)
    return names, config


def galaxy_parameters(config, row, curve):
    scales = [x*(row["disk_scale_kpc"] if config["scale_mode"] == "disk_scale_multiple" else 1)
              for x in config["scale_values"]]
    mass = config["ref_mass"] or float(curve[-1, 0]*max(curve[:, 1])**2/G)
    if mass <= 0 or not np.isfinite([mass, *scales]).all():
        raise ValueError("invalid derived basis normalization")
    return {"scales_kpc": scales, "ref_mass_msun": mass,
            "max_components_effective": min(config["max_components"], len(curve)-1, len(scales)*len(config["states"].split(","))),
            "r_min_kpc": float(curve[0, 0]), "r_max_kpc": float(curve[-1, 0])}


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def software_fingerprint(binary):
    paths = sorted({*ROOT.glob("scripts/horb_sparc/*.py"), *ROOT.glob("scripts/horb_signed/*.py"),
                    *ROOT.glob("scripts/horb_common/*.py"), *ROOT.glob("scripts/plot_signed*.py"),
                    ROOT/"scripts/fit_signed_density.py", ROOT/"scripts/fit_sparc.py"})
    plotting = {}
    for name in ("matplotlib", "Pillow"):
        try:
            plotting[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            plotting[name] = None
    return {"files": {str(p.relative_to(ROOT)): digest(p) for p in paths},
            "curve_binary_sha256": digest(binary), "numpy": np.__version__, "scipy": scipy.__version__,
            "python": sys.version, "plotting": plotting}


def completed(status_path, key):
    try:
        status = json.loads(status_path.read_text())
        if status.get("status") != "complete" or status.get("cache_key") != key or not status.get("outputs"):
            return None
        for relative, expected in status["outputs"].items():
            path = (ROOT/relative).resolve()
            if not path.is_relative_to(ROOT) or digest(path) != expected:
                return None
        return status
    except (OSError, ValueError, KeyError, TypeError):
        return None


@contextmanager
def run_lock(directory):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory/".lock").open("a") as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as e:
            raise ValueError("this SPARC run is already active; use another --name or wait for it") from e
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


class Children:
    def __init__(self):
        self.cancelled = threading.Event()
        self.lock = threading.Lock()
        self.active = set()

    def run(self, argv, log):
        env = {**os.environ, "MPLBACKEND": "Agg", "OPENBLAS_NUM_THREADS": "1",
               "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"}
        argv = [str(x) for x in argv]
        log.write("$ "+shlex.join(argv)+"\n"); log.flush()
        with self.lock:
            if self.cancelled.is_set():
                raise RuntimeError("batch interrupted")
            child = subprocess.Popen(argv, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
            self.active.add(child)
        try:
            code = child.wait()
        finally:
            with self.lock:
                self.active.discard(child)
        if code:
            raise RuntimeError(f"{Path(argv[0]).name} exited with status {code}; see run.log")

    def cancel(self):
        self.cancelled.set()
        with self.lock:
            for child in self.active:
                try:
                    child.terminate()
                except ProcessLookupError:
                    pass  # A child may finish between registration and cancellation.


def fit_one(name, args, config, catalog, curves, source, key, reports, plots, children):
    directory = reports/name
    directory.mkdir(parents=True, exist_ok=True)
    status_path = directory/"status.json"
    parameters = galaxy_parameters(config, catalog[name], curves[name])
    state = {"galaxy": name, "cache_key": key, "status": "running", "parameters": parameters}
    save_json(status_path, state)
    started = time.monotonic()
    try:
        with (directory/"run.log").open("w") as log:
            prepare_galaxy(directory, catalog[name], curves[name], source, args.ml_disk, args.ml_bulge)
            lo, hi = parameters["r_min_kpc"], parameters["r_max_kpc"]
            children.run([args.curve_bin, "signed_basis", "--out", directory/"basis.csv",
                          "--states", args.states, "--scales", ",".join(format(v, ".17g") for v in parameters["scales_kpc"]),
                          "--ref-mass", parameters["ref_mass_msun"], "--r-min", lo, "--r-max", hi,
                          "--r-step", (hi-lo)/(args.grid_points-1), "--rc-csv", directory/"rotation_curve.csv",
                          "--n-radial", args.n_radial], log)
            command = [sys.executable, ROOT/"scripts/fit_signed_density.py", "--basis-csv", directory/"basis.csv",
                       "--rc-csv", directory/"rotation_curve.csv", "--baryons-csv", directory/"baryons.csv",
                       "--dataset-meta", directory/"dataset.json", "--out-prefix", directory/"fit",
                       "--r-min", lo, "--r-max", hi, "--max-components", parameters["max_components_effective"]]
            for key_name in ("search", "budget", "beam_width", "seed", "max_abs_weight", "criterion", "ridge", "weights"):
                command.extend(["--"+key_name.replace("_", "-"), config[key_name]])
            children.run(command, log)
            outputs = [directory/f for f in ("rotation_curve.csv", "baryons.csv", "dataset.json", "basis.csv", "fit.json",
                       "fit_weights.csv", "fit_curve.csv", "fit_observations.csv", "fit_candidates.csv")]
            if args.plots != "none":
                image = plots/f"{name}_rotation.png"
                children.run([sys.executable, ROOT/"scripts/plot_signed_fit.py", "--fit-prefix", directory/"fit", "-o", image], log)
                outputs.append(image)
            if args.plots == "both":
                density = directory/"density_xz.csv"
                image = plots/f"{name}_density.png"
                children.run([args.curve_bin, "signed_density", "--weights", directory/"fit_weights.csv",
                              "--out", density, "--extent", hi, "--bins", args.density_bins], log)
                children.run([sys.executable, ROOT/"scripts/plot_signed_density.py", density, "-o", image,
                              "--title", f"SPARC {name}: signed orbital density"], log)
                outputs.extend([density, image])
            fit = json.loads((directory/"fit.json").read_text())
            fields = ("k", "chi2", "chi2_per_point", "rms_kms", "selection_score", "evaluated_subsets", "failed_subsets",
                      "budget_reached", "matrix_rank", "condition_number", "coefficient_bound_reached",
                      "net_mass_msun", "positive_mass_msun", "negative_mass_msun", "components")
            state.update(status="complete", metrics={k: fit[k] for k in fields},
                         outputs={str(p.relative_to(ROOT)): digest(p) for p in outputs})
    except Exception as e:
        state.update(status="interrupted" if children.cancelled.is_set() else "failed", error=str(e))
    state["seconds"] = time.monotonic()-started
    save_json(status_path, state)
    return state


def check_run_config(path, config):
    if path.exists():
        old = json.loads(path.read_text())
        if old.get("schema_version") != 1 or old.get("config") != config:
            raise ValueError("run settings differ; choose a new --name (even --force does not mix settings)")


def main():
    args = parser().parse_args()
    catalog, curves, source = load_data()
    names, config = validate(args, catalog, curves)
    reports, plots = ROOT/"reports/sparc"/args.name, ROOT/"plots/sparc"/args.name
    check_run_config(reports/"run.json", config)
    for name in names:
        galaxy_parameters(config, catalog[name], curves[name])
    if args.validate_only:
        print(f"SPARC validated: {len(names)} galaxies, {sum(len(curves[n]) for n in names)} points; no quality cuts")
        return 0
    if args.curve_bin is None or not args.curve_bin.is_file() or not os.access(args.curve_bin, os.X_OK):
        raise ValueError("run through ./run.sh --release do sparc-fit (or supply an executable --curve-bin)")
    args.curve_bin = args.curve_bin.resolve()
    software = software_fingerprint(args.curve_bin)
    key = signature({"config": config, "source": source, "software": software})
    children = Children()
    with run_lock(reports):
        check_run_config(reports/"run.json", config)
        plots.mkdir(parents=True, exist_ok=True)
        old = json.loads((reports/"run.json").read_text()) if (reports/"run.json").exists() else {}
        selected = sorted(set(old.get("selected_galaxies", [])) | set(names))
        run = {"schema_version": 1, "name": args.name, "config": config, "source": source, "software": software,
               "cache_key": key, "selected_galaxies": selected,
               "notes": ["CPU azimuthal-mean orbital gravity; signed density weights.",
                         "Best candidate found by the configured search, not a claim of a global optimum.",
                         "chi2_per_point is in-sample chi2/N, not reduced chi2 or evidence for a physical mechanism."]}
        save_json(reports/"run.json", run)
        inventory(reports, catalog, curves, sorted(catalog))
        statuses = {}
        for name in selected:
            status = completed(reports/name/"status.json", key)
            if status is not None and not (args.force and name in names):
                statuses[name] = status
            else:
                statuses[name] = {"galaxy": name, "status": "pending"}
        todo = [name for name in names if statuses[name]["status"] != "complete"]
        publish(reports, plots, run, catalog, curves, statuses)
        print(f"SPARC: {len(names)-len(todo)} verified cached; {len(todo)} to fit; jobs={args.jobs}", flush=True)
        pool = ThreadPoolExecutor(max_workers=args.jobs)
        futures = {pool.submit(fit_one, name, args, config, catalog, curves, source, key, reports, plots, children): name for name in todo}
        interrupted = False
        try:
            for i, future in enumerate(as_completed(futures), 1):
                name = futures[future]
                try:
                    status = future.result()
                except Exception as e:
                    status = {"galaxy": name, "status": "failed", "error": str(e)}
                statuses[name] = status
                detail = f"RMS={status['metrics']['rms_kms']:.3f} km/s" if status["status"] == "complete" else status.get("error", "")
                print(f"[{i}/{len(todo)}] {name}: {status['status']} {detail}", flush=True)
                publish(reports, plots, run, catalog, curves, statuses)
        except KeyboardInterrupt:
            interrupted = True
            children.cancel()
            for future in futures:
                future.cancel()
            print("Interrupted; stopping children. Run the same command to resume.", flush=True)
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
            for future, name in futures.items():
                if future.done() and not future.cancelled():
                    try:
                        statuses[name] = future.result()
                    except Exception as e:
                        statuses[name] = {"galaxy": name, "status": "failed", "error": str(e)}
            publish(reports, plots, run, catalog, curves, statuses)
        if args.plots != "none":
            overview(plots/"overview.png", args.name, statuses)
        good = sum(statuses[name]["status"] == "complete" for name in names)
        print(f"Completed {good}/{len(names)} requested galaxies. Summary: {reports/'summary.csv'}")
        print(f"Gallery: {plots/'index.html'}")
        return 130 if interrupted else (0 if good == len(names) else 1)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        print(f"SPARC: {error}", file=sys.stderr)
        sys.exit(1)

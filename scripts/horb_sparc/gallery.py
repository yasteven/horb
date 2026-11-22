"""Human-readable batch inventory; keep failed/unselected galaxies visible."""
from collections import Counter
import html
import os
from pathlib import Path
from fit_signed_density import atomic_write, write_csv
from horb_sparc import save_json


def publish(reports, plots, run, catalog, curves, statuses):
    rows = []
    for name in sorted(catalog):
        status = statuses.get(name, {"status": "not_selected"})
        c, curve = catalog[name], curves[name]
        row = {"galaxy": name, "status": status["status"], "quality": c["quality"],
               "inclination_deg": c["inclination_deg"], "distance_mpc": c["distance_mpc"],
               "points": len(curve), "r_min_kpc": float(curve[0, 0]), "r_max_kpc": float(curve[-1, 0]),
               "max_components_effective": status.get("parameters", {}).get("max_components_effective", "")}
        metrics = status.get("metrics", {})
        for k in ("k", "rms_kms", "chi2", "chi2_per_point", "evaluated_subsets", "failed_subsets",
                  "budget_reached", "matrix_rank", "condition_number", "coefficient_bound_reached",
                  "positive_mass_msun", "negative_mass_msun", "net_mass_msun"):
            row[k] = metrics.get(k, "")
        row["components"] = "; ".join(f"{t['state']}@{t['a0_star_kpc']:.6g}:{t['weight']:+.8g}"
                                      for t in metrics.get("components", []))
        row["error"] = status.get("error", "")
        rows.append(row)
    write_csv(reports/"summary.csv", list(rows[0]), [list(row.values()) for row in rows])
    counts = dict(Counter(row["status"] for row in rows))
    save_json(reports/"summary.json", {"schema_version": 1, "run": run["name"], "counts": counts, "galaxies": rows})
    esc = lambda x: html.escape(str(x), quote=True)
    links = lambda path: esc(os.path.relpath(path, plots))
    parts = [f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HORB · SPARC · {esc(run['name'])}</title><style>
body{{font:16px system-ui;max-width:1500px;margin:32px auto;padding:0 20px;color:#e5edf5;background:#101923}}
a{{color:#73d7eb}}input{{font:inherit;padding:10px;width:320px;max-width:90%}}table{{border-collapse:collapse;width:100%;font-size:14px}}
th,td{{text-align:left;padding:9px;border-bottom:1px solid #344350}}th{{position:sticky;top:0;background:#182733}}
.muted{{color:#b3c2d0}}.error{{color:#ffa598}}img{{max-width:100%}}.scroll{{overflow:auto}}</style>
<h1>SPARC rotation-curve fits</h1><p>Run <b>{esc(run['name'])}</b> · {esc(', '.join(f'{n} {k}' for k,n in sorted(counts.items())))}</p>
<p>All 175 galaxies remain in this inventory. Published gas, disk and bulge; fixed M/L = {run['config']['ml_disk']:g}, {run['config']['ml_bulge']:g}.
{esc(run['config']['search'])} search, up to {run['config']['max_components']} signed terms; small datasets use at most N−1.</p>
<p class="muted">Scores use the fitted points (χ²/N, not reduced χ²). Basis choice, baryons, quality and inclination affect interpretation.
Q=1 high, Q=2 medium, Q=3 low. No quality or inclination cuts. Negative orbital weights are allowed; no jet data are used.</p>
<p><a href="{links(reports/'summary.csv')}">Summary CSV</a> · <a href="{links(reports/'run.json')}">Run settings and provenance</a>''']
    if run["config"]["plots"] != "none":
        parts.append(' · <a href="overview.png">Overview plot (written at batch end)</a>')
    parts.append('</p><label>Filter by galaxy or status <input id="filter" placeholder="e.g. NGC2403"></label><div class="scroll"><table><thead><tr><th>Galaxy</th><th>Status</th><th>Q</th><th>i (°)</th><th>N</th><th>Terms</th><th>RMS (km/s)</th><th>χ²/N</th><th>Results</th></tr></thead><tbody>')
    for row in rows:
        name, status = row["galaxy"], row["status"]
        result = ""
        if status == "complete":
            result = f'<a href="{links(reports/name/"fit.json")}">Report</a> · <a href="{links(reports/name/"fit_weights.csv")}">Weights</a>'
            if run["config"]["plots"] != "none":
                result += f' · <a href="{esc(name)}_rotation.png">Rotation</a>'
            if run["config"]["plots"] == "both":
                result += f' · <a href="{esc(name)}_density.png">Density</a>'
            if row["coefficient_bound_reached"]:
                result += ' <span class="error">weight bound reached</span>'
        elif row["error"]:
            result = f'<span class="error">{esc(row["error"])}</span> · <a href="{links(reports/name/"run.log")}">Log</a>'
        fmt = lambda v: f"{v:.3f}" if isinstance(v, (float, int)) else esc(v)
        parts.append(f'<tr><td>{esc(name)}</td><td>{esc(status)}</td><td>{row["quality"]}</td><td>{row["inclination_deg"]:g}</td><td>{row["points"]}</td><td>{row["k"]}</td><td>{fmt(row["rms_kms"])}</td><td>{fmt(row["chi2_per_point"])}</td><td>{result}</td></tr>')
    parts.append('''</tbody></table></div><script>
document.querySelector('#filter').addEventListener('input',function(){const q=this.value.toLowerCase();document.querySelectorAll('tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});
</script></html>''')
    atomic_write(plots/"index.html", lambda f: f.write("\n".join(parts)))


def overview(path, name, statuses):
    import numpy as np
    import matplotlib.pyplot as plt
    good = [s for s in statuses.values() if s["status"] == "complete"]
    if not good:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    rms = [s["metrics"]["rms_kms"] for s in good]
    axes[0].hist(rms, bins=min(25, max(5, len(good)//5)), color="#257b8b", edgecolor="white")
    axes[0].set(xlabel="In-sample RMS [km/s]", ylabel="Galaxies", title=f"Median RMS: {np.median(rms):.2f} km/s")
    terms = Counter(s["metrics"]["k"] for s in good)
    axes[1].bar(sorted(terms), [terms[k] for k in sorted(terms)], color="#bb6941")
    axes[1].set(xlabel="Selected orbital terms", ylabel="Galaxies", xticks=sorted(terms), title="Configured sparse search")
    fig.suptitle(f"HORB / SPARC · {name} · {len(good)} completed fits")
    fig.text(.5, .015, "Fixed published baryonic components; all quality flags included. Scores describe the fitted data.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.95)); fig.savefig(path, dpi=180); plt.close(fig)

#!/usr/bin/env python3
"""Prepare all 175 SPARC curves and their published gas/disk/bulge force tables offline."""
import argparse
from pathlib import Path
from horb_sparc import (ROOT, load_data, select_galaxies, mass_to_light,
                        prepare_galaxy, inventory)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, default=ROOT/"reports/sparc_inputs")
    p.add_argument("--galaxies", default="all")
    p.add_argument("--ml-disk", type=float, default=0.5)
    p.add_argument("--ml-bulge", type=float, default=0.7)
    args = p.parse_args()
    mass_to_light(args.ml_disk, args.ml_bulge)
    catalog, curves, source = load_data()
    names = select_galaxies(args.galaxies, catalog)
    for name in names:
        prepare_galaxy(args.out_dir/name, catalog[name], curves[name], source, args.ml_disk, args.ml_bulge)
    inventory(args.out_dir, catalog, curves, names)
    print(f"Prepared {len(names)} galaxies / {sum(len(curves[n]) for n in names)} points in {args.out_dir}")


if __name__ == "__main__":
    main()

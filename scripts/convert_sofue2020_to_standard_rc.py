#!/usr/bin/env python3

from pathlib import Path
import csv

src = Path("data/milky_way/sofue_2020_unified_rc.csv")
dst = Path("data/milky_way/sofue_2020_standard_rc.csv")

dst.parent.mkdir(parents=True, exist_ok=True)

with src.open() as f_in, dst.open("w", newline="") as f_out:
    reader = csv.DictReader(f_in)
    writer = csv.DictWriter(
        f_out,
        fieldnames=[
            "R_kpc",
            "Vobs_kms",
            "eV_kms",
            "Vgas_kms",
            "Vdisk_kms",
            "Vbul_kms",
            "source",
        ],
    )

    writer.writeheader()

    for row in reader:
        writer.writerow({
            "R_kpc": row["R_kpc"],
            "Vobs_kms": row["v_kms"],
            "eV_kms": row["v_err_kms"],
            "Vgas_kms": "",
            "Vdisk_kms": "",
            "Vbul_kms": "",
            "source": "Sofue2020",
        })

print(f"wrote {dst}")

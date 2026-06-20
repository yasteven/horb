#!/usr/bin/env python3

from pathlib import Path
import numpy as np

out = Path("data/milky_way/mw_starter_total_curve.csv")
out.parent.mkdir(parents=True, exist_ok=True)

# Starter total circular velocity curve.
#
# This is NOT final raw data.
# It is a clean first-pass target:
#   ~235 km/s at R=8.2 kpc
#   gently declining outer curve
#   usable immediately for pipeline testing.
#
# Later replace with Eilers/Sofue tabulated values.
r = np.array([
    5.0, 6.0, 7.0, 8.2, 10.0,
    12.0, 15.0, 18.0, 20.0, 22.0, 25.0,
    30.0, 40.0, 50.0, 60.0, 80.0, 100.0
])

v = np.array([
    240.0, 238.0, 236.0, 235.0, 232.0,
    229.0, 224.0, 219.0, 216.0, 212.0, 207.0,
    200.0, 185.0, 170.0, 155.0, 130.0, 110.0
])

# Loose placeholder uncertainties for visual overlay.
err = np.array([
    8.0, 7.0, 6.0, 5.0, 5.0,
    5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    12.0, 15.0, 18.0, 20.0, 25.0, 30.0
])

with out.open("w") as f:
    f.write("R_kpc,v_kms,v_err_kms,source\n")
    for rr, vv, ee in zip(r, v, err):
        f.write(f"{rr:.6f},{vv:.6f},{ee:.6f},starter_mw_total\n")

print(f"wrote {out}")

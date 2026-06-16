#!/usr/bin/env python3

from pathlib import Path
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

paths = sorted(Path("plots").glob("dz2_total_a1.5_m*_vs_sofue2020.png"))

if not paths:
    raise SystemExit("No plots found matching plots/dz2_total_a1.5_m*_vs_sofue2020.png")

cols = 2
rows = (len(paths) + cols - 1) // cols

fig, axes = plt.subplots(rows, cols, figsize=(cols * 7, rows * 5))

if rows == 1:
    axes = [axes]

flat = []
for row in axes:
    try:
        flat.extend(row)
    except TypeError:
        flat.append(row)

for ax, path in zip(flat, paths):
    ax.imshow(mpimg.imread(path))
    ax.set_title(path.stem)
    ax.axis("off")

for ax in flat[len(paths):]:
    ax.axis("off")

fig.tight_layout()
out = Path("plots/contact_sheet_horb_mass_scan_vs_sofue2020.png")
fig.savefig(out, dpi=160)
plt.close(fig)

print(f"wrote {out}")

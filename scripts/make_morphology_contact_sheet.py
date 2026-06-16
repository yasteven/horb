#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Make HORB morphology contact sheet.")
    parser.add_argument(
        "--mode",
        choices=["log", "linear", "threshold015"],
        default="threshold015",
    )
    parser.add_argument("--mass", default="1e11")
    parser.add_argument(
        "--a0",
        nargs="+",
        default=["1.0", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "2.0"],
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=None,
    )
    args = parser.parse_args()

    images = []
    labels = []

    for a0 in args.a0:
        path = Path("plots") / f"dz2_xz_a{a0}_m{args.mass}_{args.mode}.png"
        if not path.exists():
            print(f"missing: {path}")
            continue

        images.append(mpimg.imread(path))
        labels.append(f"a0={a0} kpc")

    if not images:
        raise SystemExit("no images found")

    n = len(images)
    cols = 3
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))

    if rows == 1:
        axes = [axes]

    flat_axes = []
    for row in axes:
        try:
            flat_axes.extend(row)
        except TypeError:
            flat_axes.append(row)

    for ax, img, label in zip(flat_axes, images, labels):
        ax.imshow(img)
        ax.set_title(label)
        ax.axis("off")

    for ax in flat_axes[len(images):]:
        ax.axis("off")

    fig.tight_layout()

    if args.out is None:
        args.out = Path("plots") / f"morphology_contact_sheet_{args.mode}_m{args.mass}.png"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=160)
    plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

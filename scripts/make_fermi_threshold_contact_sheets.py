#!/usr/bin/env python3

from pathlib import Path
import argparse
import matplotlib.image as mpimg
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mass", default="2.5e11")
    parser.add_argument("--thresholds", nargs="+", default=["005", "010", "015", "020", "025", "030", "040", "050"])
    parser.add_argument("--a0", nargs="+", default=["1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "2.0"])
    args = parser.parse_args()

    for th in args.thresholds:
        paths = []
        labels = []

        for a0 in args.a0:
            path = Path("plots") / f"dz2_xz_a{a0}_m{args.mass}_threshold{th}.png"
            if path.exists():
                paths.append(path)
                labels.append(f"a0={a0}")

        if not paths:
            print(f"no images found for threshold {th}")
            continue

        cols = 4
        rows = (len(paths) + cols - 1) // cols

        fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))

        if rows == 1:
            axes = [axes]

        flat = []
        for row in axes:
            try:
                flat.extend(row)
            except TypeError:
                flat.append(row)

        for ax, path, label in zip(flat, paths, labels):
            ax.imshow(mpimg.imread(path))
            ax.set_title(label)
            ax.axis("off")

        for ax in flat[len(paths):]:
            ax.axis("off")

        fig.tight_layout()

        out = Path("plots") / f"fermi_morphology_contact_threshold{th}_m{args.mass}.png"
        fig.savefig(out, dpi=160)
        plt.close(fig)

        print(f"wrote {out}")


if __name__ == "__main__":
    main()

"""Cumulative Simpson quadrature for the one-dimensional SWELLS radial grids.

The interval integrals use the three-point quadratic formula documented at
https://docs.scipy.org/doc/scipy-1.12.0/reference/generated/scipy.integrate.cumulative_simpson.html
(Cartwright 2017, equation 8). Implementing this small NumPy-only operation
keeps the same quadrature on Jetson installations with older SciPy.
This is deliberately a 1D API, not a replacement for every SciPy option.
"""
import numpy as np


def cumulative_simpson_1d(y, *, x, initial=0.0):
    """Integrate real samples on an increasing, possibly nonuniform grid.

Each pair of intervals uses its three-point quadratic. If a final interval
is unpaired, use the final three samples. Two samples use a trapezoid.
With ``initial=None`` the first grid point is omitted from the result;
otherwise prepend the scalar initial integral and add it to every value.
Negative samples and cumulative integrals retain their signs.
"""
    y, x = np.asarray(y, dtype=float), np.asarray(x, dtype=float)
    if y.ndim != 1 or x.ndim != 1 or y.shape != x.shape or not len(x):
        raise ValueError("quadrature needs nonempty 1D arrays of equal length")
    widths = np.diff(x)
    if not np.isfinite(y).all() or not np.isfinite(x).all() or np.any(widths <= 0):
        raise ValueError("quadrature needs finite samples and strictly increasing x")
    if initial is not None and (np.ndim(initial) != 0 or not np.isfinite(initial)):
        raise ValueError("initial integral must be a finite scalar or None")
    intervals = np.empty(len(x)-1)
    if len(x) == 2:
        intervals[0] = widths[0]*(y[0]+y[1])/2
    elif len(x) >= 3:
        a, b = widths[:-1], widths[1:]
        total = a+b
        # Integrate the interpolating quadratic over its left/right interval.
        left_ratio, right_ratio = a*a/(b*total), b*b/(a*total)
        left = a/6*((3-a/total)*y[:-2] + (3+a/total+left_ratio)*y[1:-1] - left_ratio*y[2:])
        right = b/6*((3-b/total)*y[2:] + (3+b/total+right_ratio)*y[1:-1] - right_ratio*y[:-2])
        intervals[:-1:2] = left[::2]
        intervals[1::2] = right[::2]
        intervals[-1] = right[-1]
    result = np.cumsum(intervals)
    if initial is None:
        return result
    return np.r_[float(initial), result+initial]

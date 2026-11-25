"""Independent integral checks and comparison with optional modern SciPy."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from horb_swells.quadrature import cumulative_simpson_1d


class SwellsQuadratureTests(unittest.TestCase):
    def test_nonuniform_signed_quadratic_and_uniform_cubic(self):
        x = np.array([-3., -2.7, -1., -.2, .05, .8, 1., 1.4])
        primitive = lambda r: r**3/3-r**2-4*r
        for n in (3, 4, 7, 8):
            grid = x[:n]
            actual = cumulative_simpson_1d(grid**2-2*grid-4, x=grid, initial=2.)
            np.testing.assert_allclose(actual, primitive(grid)-primitive(grid[0])+2, atol=2e-14)
        x = np.linspace(-2., 3., 21)
        actual = cumulative_simpson_1d(x**3, x=x)
        # Full two-interval Simpson panels are exact for cubics.
        np.testing.assert_allclose(actual[::2], (x[::2]**4-x[0]**4)/4, atol=1e-13)

    def test_short_inputs_and_initial_offset(self):
        np.testing.assert_equal(cumulative_simpson_1d([7.], x=[1.]), [0.])
        np.testing.assert_equal(cumulative_simpson_1d([7.], x=[1.], initial=None), [])
        np.testing.assert_equal(cumulative_simpson_1d([-2., 6.], x=[1., 4.]), [0., 6.])
        np.testing.assert_equal(cumulative_simpson_1d([-2., 6.], x=[1., 4.], initial=5.), [5., 11.])

    def test_rejects_invalid_grids_and_preserves_signed_linearity(self):
        for y, x in (([], []), ([1., 2.], [0.]), ([[1., 2.]], [[0., 1.]]),
                     ([1., 2.], [1., 1.]), ([1., 2.], [1., 0.]),
                     ([1., np.nan], [0., 1.]), ([1., 2.], [0., np.inf])):
            with self.assertRaises(ValueError):
                cumulative_simpson_1d(y, x=x)
        x = np.linspace(-8., 8., 101)
        y = np.sin(x)*np.exp(-x*x/4)
        np.testing.assert_allclose(cumulative_simpson_1d(-3*y, x=x),
                                   -3*cumulative_simpson_1d(y, x=x), atol=2e-15)

    def test_agrees_with_modern_scipy_when_available(self):
        try:
            from scipy.integrate import cumulative_simpson
        except ImportError:
            self.skipTest("native cumulative_simpson is optional; analytic tests still run")
        rng = np.random.default_rng(341)
        for n in (1, 2, 3, 4, 5, 6, 11, 12, 3199, 3200):
            for uniform in (False, True):
                x = np.arange(n, dtype=float) if uniform else np.cumsum(rng.uniform(.1, 2., n))
                y = rng.normal(size=n)
                for initial in (None, 0., 2.5):
                    np.testing.assert_allclose(cumulative_simpson_1d(y, x=x, initial=initial),
                                               cumulative_simpson(y, x=x, initial=initial),
                                               rtol=3e-13, atol=3e-13)


if __name__ == "__main__":
    unittest.main()

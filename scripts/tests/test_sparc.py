import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
import fit_sparc
from horb_sparc import (DATA, G, load_data, baryonic_components, prepare_galaxy,
                        select_galaxies, digest, save_json)
from horb_signed import baryons_function
from fit_signed_density import dataset_metadata


class SparcTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog, cls.curves, cls.source = load_data()

    def options(self, *argv):
        args = fit_sparc.parser().parse_args(list(argv))
        names, config = fit_sparc.validate(args, self.catalog, self.curves)
        return args, names, config

    def test_complete_published_sample_including_low_quality_and_negative_forces(self):
        self.assertEqual(len(self.catalog), 175)
        self.assertEqual(sum(map(len, self.curves.values())), 3391)
        self.assertEqual([sum(c['quality'] == q for c in self.catalog.values()) for q in (1, 2, 3)], [99, 64, 12])
        self.assertEqual(sum(np.any(c[:, 3:6] < 0) for c in self.curves.values()), 48)
        np.testing.assert_array_equal(self.curves['NGC2403'][0], [.16, 24.5, 2.83, 0, 23.21, 0, 1105.79, 0])
        self.assertEqual(self.catalog['CamB']['disk_scale_kpc'], .47)
        self.assertNotIn('M31', self.catalog)

    def test_all_normalized_observations_and_baryons_preserve_published_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, curve in self.curves.items():
                dest = Path(tmp)/name
                prepare_galaxy(dest, self.catalog[name], curve, self.source)
                np.testing.assert_array_equal(np.loadtxt(dest/'rotation_curve.csv', delimiter=',', skiprows=1), curve[:, :3])
                expected = curve[:, 3]*abs(curve[:, 3]) + .5*curve[:, 4]*abs(curve[:, 4]) + .7*curve[:, 5]*abs(curve[:, 5])
                np.testing.assert_allclose(baryons_function(dest/'baryons.csv')(curve[:, 0]), expected, rtol=2e-15)
                dataset_metadata(dest/'dataset.json', dest/'rotation_curve.csv', dest/'baryons.csv')

    def test_signed_baryonic_formula_scales_force_once_and_never_adds_helium_again(self):
        curve = np.array([[1, 20, 2, -2, -3, -4, 0, 0], [2, 30, 3, 2, 3, 4, 0, 0]], dtype=float)
        parts = baryonic_components(curve, .5, .7)
        np.testing.assert_allclose(parts, [[-4, -4.5, -11.2], [4, 4.5, 11.2]])
        for disk, bulge in [(float('nan'), .7), (-1, .7), (.5, float('inf'))]:
            with self.assertRaises(ValueError):
                baryonic_components(curve, disk, bulge)

    def test_source_corruption_fails_before_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)/'source'
            shutil.copytree(DATA, dest)
            file = dest/'MassModels_Lelli2016c.mrt'
            file.write_bytes(file.read_bytes()+b'\n')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                load_data(dest)

    def test_relative_scales_mass_normalization_and_small_sample_component_cap(self):
        _, _, config = self.options()
        for name in self.catalog:
            row, curve = self.catalog[name], self.curves[name]
            p = fit_sparc.galaxy_parameters(config, row, curve)
            np.testing.assert_allclose(p['scales_kpc'], np.array([.125, .25, .5, 1, 2])*row['disk_scale_kpc'])
            self.assertEqual(p['max_components_effective'], min(4, len(curve)-1))
            self.assertAlmostEqual(p['ref_mass_msun']*G/curve[-1, 0], max(curve[:, 1])**2)
        _, _, fixed = self.options('--scales', '0.5,2', '--ref-mass', '1e10')
        p = fit_sparc.galaxy_parameters(fixed, self.catalog['D512-2'], self.curves['D512-2'])
        self.assertEqual(p['scales_kpc'], [.5, 2])
        self.assertEqual(p['ref_mass_msun'], 1e10)
        self.assertEqual(p['max_components_effective'], 3)

    def test_selection_validation_and_exhaustive_budget(self):
        self.assertEqual(len(select_galaxies('all', self.catalog)), 175)
        for argv in [('--galaxies', 'M31'), ('--galaxies', 'NGC2403,NGC2403'), ('--name', '../escape'),
                     ('--scale-factors', '1,1'), ('--scales', 'nan'), ('--jobs', '0'),
                     ('--max-components', '0'), ('--states', '1s,2p_x'), ('--search', 'exhaustive'),
                     ('--budget', '1')]:
            with self.assertRaises(ValueError):
                self.options(*argv)
        self.options('--search', 'exhaustive', '--budget', '31931')

    def test_resume_verifies_every_artifact_and_settings(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(fit_sparc, 'ROOT', Path(tmp)):
            root = Path(tmp)
            (root/'fit.json').write_text('a valid result')
            status = root/'status.json'
            save_json(status, {'status': 'complete', 'cache_key': 'key', 'outputs': {'fit.json': digest(root/'fit.json')}})
            self.assertIsNotNone(fit_sparc.completed(status, 'key'))
            self.assertIsNone(fit_sparc.completed(status, 'changed-code-or-inputs'))
            (root/'fit.json').write_text('truncated')
            self.assertIsNone(fit_sparc.completed(status, 'key'))
            (root/'fit.json').unlink()
            self.assertIsNone(fit_sparc.completed(status, 'key'))
            status.write_text('{partial')
            self.assertIsNone(fit_sparc.completed(status, 'key'))
            _, _, config = self.options()
            run = root/'run.json'
            save_json(run, {'schema_version': 1, 'config': config})
            fit_sparc.check_run_config(run, config)
            with self.assertRaisesRegex(ValueError, 'new --name'):
                fit_sparc.check_run_config(run, {**config, 'ml_disk': .6})

    def test_per_run_lock_prevents_concurrent_writers(self):
        with tempfile.TemporaryDirectory() as tmp:
            with fit_sparc.run_lock(Path(tmp)):
                with self.assertRaisesRegex(ValueError, 'already active'):
                    with fit_sparc.run_lock(Path(tmp)):
                        pass
            with fit_sparc.run_lock(Path(tmp)):
                pass

    def test_child_failure_is_recorded_and_batch_continues_then_returns_nonzero(self):
        # Two real failed child processes. The second galaxy must still run,
        # both errors must survive in status/summary, and exit status must fail.
        with tempfile.TemporaryDirectory() as tmp, patch.object(fit_sparc, 'ROOT', Path(tmp)):
            root = Path(tmp)
            binary = root/'failing curve'
            binary.write_text('#!/bin/sh\nexit 71\n'); binary.chmod(0o755)
            argv = ['fit_sparc', '--name', 'failure', '--galaxies', 'D512-2,DDO064', '--plots', 'none', '--curve-bin', str(binary)]
            with patch.object(sys, 'argv', argv), patch.object(fit_sparc, 'software_fingerprint', return_value={'test': 'failure'}):
                self.assertEqual(fit_sparc.main(), 1)
            for name in ['D512-2', 'DDO064']:
                status = json.loads((root/'reports/sparc/failure'/name/'status.json').read_text())
                self.assertEqual(status['status'], 'failed')
                self.assertIn('71', status['error'])
                self.assertNotIn('outputs', status)
            summary = json.loads((root/'reports/sparc/failure/summary.json').read_text())
            self.assertEqual(summary['counts'], {'failed': 2, 'not_selected': 173})


if __name__ == '__main__':
    unittest.main()

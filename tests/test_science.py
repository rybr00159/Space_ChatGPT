"""Small integrity/time-index tests; these do not establish scientific validity."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from archive import verify_fits_checksums
from search import flatten, mask_known, read_curve
from validate import local_depth


class ScientificUtilities(unittest.TestCase):
    def test_checksum_without_datasum_and_detects_corruption(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "raw.fits"
            hdu = fits.PrimaryHDU(np.arange(20, dtype=np.int16))
            hdu.add_checksum(override_datasum=True)
            hdu.writeto(path)
            with fits.open(path) as hdul:
                self.assertNotIn("DATASUM", hdul[0].header)
                self.assertEqual(verify_fits_checksums(path, hdul)[0]["status"], "verified_raw_ones_complement")
                offset = hdul[0].fileinfo()["datLoc"]
            raw = bytearray(path.read_bytes())
            raw[offset + 1] ^= 1
            path.write_bytes(raw)
            with fits.open(path) as hdul:
                with self.assertRaises(ValueError):
                    verify_fits_checksums(path, hdul)

    def test_read_curve_filters_and_sorts_arrays_together(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lc.fits"
            columns = [fits.Column(name="TIME", format="D", array=[3., 1., 2., 4., np.nan]),
                       fits.Column(name="PDCSAP_FLUX", format="E", array=[30., 10., 20., 40., 50.]),
                       fits.Column(name="PDCSAP_FLUX_ERR", format="E", array=[3., 1., 2., 4., 5.]),
                       fits.Column(name="QUALITY", format="J", array=[0, 0, 1, 0, 0])]
            fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns)]).writeto(path)
            t, f, e, _ = read_curve(path)
            np.testing.assert_equal(t, [1., 3., 4.])
            np.testing.assert_equal(f, [10., 30., 40.])
            np.testing.assert_equal(e, [1., 3., 4.])
            self.assertTrue(f.dtype.isnative)

    def test_bjd_epoch_and_nonchronological_masks(self):
        t = np.array([12., 11., 10., 10.03, 9.9])
        planet = [{"period": 2., "epoch_bjd": 2457010., "duration_hours": 1.}]
        np.testing.assert_equal(mask_known(t, planet), [False, True, False, False, True])

    def test_flatten_does_not_bridge_large_gap(self):
        t = np.r_[np.arange(40) / 1440, 10 + np.arange(40) / 1440]
        f = np.r_[np.ones(40) * 100, np.ones(40) * 1000]
        e = f * .001
        tb, fb, eb = flatten(t, f, e, width=1.)
        np.testing.assert_allclose(fb, 1.)
        self.assertEqual(len(tb), len(fb))
        self.assertTrue(np.all(eb > 0))
        self.assertFalse(np.any((tb > .1) & (tb < 9.9)))

    def test_local_baseline_preserves_box_on_linear_trend(self):
        t = np.arange(0, 5, .002)
        epoch, period, duration, depth = .3, 1.3, .08, .0012
        inside = np.abs((t - epoch + period/2) % period - period/2) < duration/2
        f = 1 + .0002 * t - depth * inside
        result = local_depth(t, f, np.full(len(t), .0001), period, epoch, duration)
        self.assertEqual(result["n_events"], 4)
        self.assertAlmostEqual(result["depth"], depth, places=10)


if __name__ == "__main__":
    unittest.main()

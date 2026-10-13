# Known-host pilot: independent residual vetting

Analysis date: 2026-10-13T16:49:29.786054+00:00

The strongest residuals on TOI-1064 (TIC 79748331) and TOI-4342 (TIC 354944123) were re-examined in raw SAP and PDCSAP, both median-detrending widths, individual usable segments, and other observed sectors. The prominent long-period screen peaks are dominated by different quality-gap or sector-boundary excursions. These residual periods are not supported as new planets.

Time units were verified directly from every FITS: `TIMESYS=TDB`, `BJDREFI=2457000`, `BJDREFF=0`. Catalogue BJD_TDB epochs were converted to BTJD once. `CROWDSAP` was read from the light-curve extension; high values here are aperture contamination estimates, not blend exclusion.

| TIC | Sector | Original top period (d) | Original diagonal S/N | Top S/N after 0.5d gap/edge trim, width 1d | CROWDSAP |
| --- | --- | --- | --- | --- | --- |
| 79748331 | 94 | 12.030505 | 20.95 | 7.07 | 0.9612 |
| 79748331 | 104 | 12.421142 | 11.27 | 6.91 | 0.9603 |
| 79748331 | 105 | 5.693564 | 9.32 | 7.23 | 0.9392 |
| 354944123 | 101 | 11.540718 | 10.46 | 5.03 | 0.9896 |
| 354944123 | 102 | 11.150664 | 10.03 | 8.53 | 0.9800 |
| 354944123 | 103 | 6.157538 | 5.80 | 4.44 | 0.9828 |

## Event-level evidence

- TIC 79748331, sector 94: deepest predicted event at BTJD 3881.045291, local depth 3735 ppm, nearest usable-segment boundary 0.0415 d. Fixed original-ephemeris S/N changes from 20.94873599738573 to None after trimming. [Diagnostic plot](tic79748331_s0094_boundary_diagnostics.png).
- TIC 79748331, sector 104: deepest predicted event at BTJD 4197.931256, local depth 1639 ppm, nearest usable-segment boundary 0.0616 d. Fixed original-ephemeris S/N changes from 11.270306868472233 to 1.6650183410445039 after trimming. [Diagnostic plot](tic79748331_s0104_boundary_diagnostics.png).
- TIC 79748331, sector 105: deepest predicted event at BTJD 4225.769072, local depth 1072 ppm, nearest usable-segment boundary 0.0968 d. Fixed original-ephemeris S/N changes from 9.323699106325916 to 6.093659041967962 after trimming. [Diagnostic plot](tic79748331_s0105_boundary_diagnostics.png).
- TIC 354944123, sector 101: deepest predicted event at BTJD 4101.477856, local depth 3219 ppm, nearest usable-segment boundary 0.0343 d. Fixed original-ephemeris S/N changes from 9.59705025333912 to None after trimming. [Diagnostic plot](tic354944123_s0101_boundary_diagnostics.png).
- TIC 354944123, sector 102: deepest predicted event at BTJD 4127.543371, local depth 2889 ppm, nearest usable-segment boundary 0.0680 d. Fixed original-ephemeris S/N changes from 10.030721546815698 to None after trimming. [Diagnostic plot](tic354944123_s0102_boundary_diagnostics.png).
- TIC 354944123, sector 103: deepest predicted event at BTJD 4177.499014, local depth 2267 ppm, nearest usable-segment boundary 0.0698 d. Fixed original-ephemeris S/N changes from 5.797616355078229 to None after trimming. [Diagnostic plot](tic354944123_s0103_boundary_diagnostics.png).

## Scope and limitations

Each original period/epoch is propagated without retuning to the other sectors in `known_host_vetting.json`. Individual predicted events, masks, local depths and uncertainties are retained, including missing predictions. This is a retrospective consistency check; the other sectors were already searched and do not constitute independent prospective validation. Per-sector periods do not justify an ephemeris uncertainty small enough for a distant-season null result to exclude every possible period alias.

The 0.5d exclusion applies to every contiguous good-quality segment, including 0.217d gaps. It tests the hypothesis that instrumental edge excursions dominate the BLS score. True transits can occur near gaps; trimming sensitivity is reduced, and this diagnostic rejects the promoted residual interpretation rather than proving either star lacks additional planets.

Odd/even statistics with only one supported event are undefined in scientific practice even when the screen reports zero. The TOI-4342 sector-102 ~9.379d interpretation also alternates deep events with an essentially absent middle event; inspect individual events rather than counting an available time window as a transit. Exact figures for both detrending widths are preserved in JSON.

No formal transit fit, localization, or false-positive probability was justified for these rejected artifacts. Known ephemerides have also been audited with fixed-phase and freely shifted-phase depths; a free-phase detection is not automatically timing variation, and must be consistent across sectors. TOI-4342 has a previously reported ~47.5d RV candidate according to the parallel literature review; a future long-period transit search must treat that as prior recognition.

Reproduce with `python src/vet_known_host_leads.py` after downloading the same FITS products listed by SHA256 in each sector JSON. The script writes only this folder and does not modify the archive manifest or screen implementation. Raw data are retained in the project cache.

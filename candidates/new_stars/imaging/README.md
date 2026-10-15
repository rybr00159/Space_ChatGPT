# Image-level search: no validated new star

Twelve pilot fields and 64 expanded fields were analyzed using original,
model and residual Legacy Survey g/r/z images. A total 1035 recurring peaks
were screened, including 495 Gaia-missing but deeper-catalogued detections.
Peaks are not equivalent to stars.

The final expanded-sample compact peak at RA 202.550980,Dec−8.576843 lies
inside native g/r/z saturation masks (value 18462). [Its plot](expanded64_survivor.png)
shows the bright-star subtraction structure. It cannot support a reliable
new-star detection; this does not prove no object exists behind that mask.

One faint pilot peak at RA 216.428278,Dec+23.612155 remains unverified and
broader than the nominal PSF. [The field](pilot_peak.png),
[independent image comparison](independent_peak_check.png) and
[stationary/common-motion checks](independent_motion_check.png) retain the
evidence. PS1 and SDSS checks did not establish another detection. It could
be an extended/variable source or artifact; it is not promoted as a star.

All 218 intermediate plots (186265828 bytes) were moved intact to the cache;
[their SHA manifest](../../../data/new_stars/imaging/plot_preservation.json)
records exact identities. Four selected plots remain in Git. Numerical
tables, snapshots of earlier methods, rejected peaks, catalogue matches,
injections and failed queries are in `data/new_stars/imaging/`.

Reproduce bulk plots outside Git:

```sh
python scripts/run_uncatalogued_imaging.py --expanded --run-name expanded64_v4
```

Use the archived method snapshot for an exact earlier version. Injection
recovery measures only the detection stage, not full completeness. Close
faint companions remain poorly accessible; see [the detailed image report](../../../reports/new_stars/imaging.md).

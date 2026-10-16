# Multi-sector recognition and targeted context

This directory records prior-report recognition separately from the blind
light-curve search. A missing catalog match is not a discovery or an exhaustive
literature result. See `data/multisector/catalogs/README.md` for snapshot scopes
and the reproducible API.

`signals.csv.gz` is an annotated snapshot of the growing survey queue. Its
adjacent manifest records the exact input hash, retained compressed original
input, row count, output hash, method hashes and catalog provenance. It can
lag the running survey; compare its manifest to the survey's final export and
rerun recognition at a checkpoint. `signals_initial_4080.csv` preserves the
first, explicitly preliminary annotation. `prioritized_unknown.csv` is a
numerical follow-up queue, not a credible-candidate list. Different detrending
widths use the same photons and are not independent detections.

`positive_controls.json` records 15 checks covering actual confirmed/RV planets,
CTOI/RAVEN recognition, reported-manuscript signals, harmonics, actual-epoch
masking, nearby contamination, known-host additional-period eligibility and
service/schema errors. Run them with:

```sh
python src/multisector_recognition.py --controls
```

`live_context_summary.json` derives finite public ExoFOP, VSX, SIMBAD and exact
TIC-text arXiv checks. Failed arXiv requests remain failures. The immutable
source responses and manifests are in timestamped
`data/multisector/catalogs/live/` directories. `deep_context_summary.json` adds
actual Gaia properties/classifications/NSS rows, the published ASAS-SN 2021
catalog, and a 60-arcsecond SIMBAD alias/bibliography search for selected leads.
References retain each source's name and coordinates: nearby objects are not
silently assigned to the target. Valid retained responses are reused by URL
and SHA-256 on later deep checks.

## Useful falsification context

The following are prior-context findings, not new planet or stellar discoveries:

| TIC | Actual public context | Consequence for vetting |
|---|---|---|
| 366489851 | A known RRAB source, Gaia DR3 2250689002001074688, is 21.34 arcsec away with VSX period 0.492685 days. | The searched approximately 0.98339-day peak is close to twice a known nearby variable's period; test contamination. |
| 441739038 | Same-position ASASSN-V J171724.73+720637.3, classified ROT, published period 8.1887 days. | Stellar variability is already reported; it does not alone identify every searched peak. |
| 261083172 | Gaia DR3 1956652654701605504, classified DSCT/GDOR/SXPHE, lies 16.32 arcsec away. | Nearby pulsation is an alternative source of aperture variability. |
| 395316864 | Same-position Gaia DR3 2772242578858228480 is VSX class RS. | The target's variability is already catalogued. |
| 468714627 | A DSCT/GDOR/SXPHE source 23.12 arcsec away has period 0.36818 days. | The approximately 1.843-day searched period is close to five times that period. |
| 329083521 | SIMBAD identifies emission-line star G 262-38, with 21 references. | Activity is a relevant alternative to isolated dips near positive excursions. |

TIC **368248545** has an especially useful concrete contaminant: the known
EW eclipsing binary **ZTF J043919.48+501122.7 / ATO J069.8311+50.1896**, at
24.2524 arcsec, has published VSX period **0.321179 days**. Gaia DR3 source
260154557196225536 independently has class ECL. The strongest TESS modulation
measured by follow-up, approximately 0.1605405 days, is within 0.031% of half
that binary period. The approximately 1.124-day BLS peak can therefore be a
multiple of the neighboring binary's half-period. This does not replace pixel
localization, but it makes “pulsation on the target” an unnecessarily specific
interpretation of the modulation. Actual SIMBAD bibliography rows point to
ATLAS 2018 (`2018AJ....156..241H`), ZTF 2020 (`2020ApJS..249...18C`), ZTF
anomalies 2022 (`2022ApJ...932..118C`) and the Gaia 2023 EB catalog
(`2023A&A...674A..16M`), rather than relying on one catalog's absence.

TIC **302550769**, Gaia DR3 521607117926177408, has a **published NSS astrometric
Orbital solution** with period **161.673704 +/- 0.333587 days**, eccentricity
0.14239 +/- 0.05391 and RUWE 2.284. The source is already in an astrometric
non-single-star table. Its raw NSS record also has `flags=64`, `efficiency=0`,
and `goodness_of_fit=6.34`; these are retained and have not been decoded or
independently validated. No companion mass or confirmed stellar nature follows
from this check, and the 161.67-day solution does not identify the separate
approximately 0.22458-day photometric modulation. The final deeper snapshot
successfully checks Gaia, ASAS-SN and SIMBAD; earlier rejected Gaia requests
are preserved as errors, never interpreted as source absence.

The live checks do not exhaust all aliases, all literature, all TCEs, current
ground-based catalogs or source-localization tests. Physical signal validation
and cross-sector ephemeris tests belong to `candidates/multisector/`.

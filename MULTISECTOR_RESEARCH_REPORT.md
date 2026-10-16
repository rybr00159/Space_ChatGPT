# Broad multi-sector exoplanet investigation

Started 15 October 2026 (America/Halifax; machine manifests use UTC 16 October).
**The new search is in progress. No credible new planet candidate is established.**
This phase responds to the user's authorization to expand across sectors after
the TIC 4206066 audit exposed the earlier survey's restricted coverage.

## Fixed first sample and actual execution

The first tranche selects 1,024 sector-target observations in each of QLP
Sectors **83, 86, 95 and 98**: **4,096 observations of 3,613 distinct TICs**.
These are selected counts, not completed-analysis counts. There are 321 TICs
in common with the older Sector 1 period-searched sample; 3,292 selected stars
were outside that older sample. Counts from different phases must not be added
without exact-ID deduplication.

The two northern and two southern fields were chosen using archive availability,
timing warnings, storage and repeat coverage before seeing their flux. For each
complete official inventory, 8,192 TIC IDs were ordered by a fixed SHA-256 salt.
Their public TIC metadata were retrieved. Eligibility requires finite Tmag<=12.5,
Teff 3000–6500 K and radius 0.15–2.5 solar radii; the first 1,024 eligible hash-ordered
rows in each sector form the immutable sample. These catalogue estimates need
survivor-specific verification. The sample is not a census of any whole sector.

See [archive selection and version evidence](reports/multisector/archive_probe.md),
[fixed selection](data/multisector/selection.csv) and
[selection manifest](data/multisector/selection.json).
Current actual completion counts live in [the survey summary](results/multisector/summary.json).
Raw download counts, insufficient inputs, signal peaks and unique stars remain
separate. Raw products remain cached; none is deleted to fit more targets.

## Changes from the old search

The original Sector 1 source and results remain frozen. The new wrapper calls
the same verified core with two gap-separated median widths (0.75 and 1.5 days),
eight retained peaks, refinement around 32 coarse seeds and more duration trials.
It retains recurring peaks down to local diagonal S/N 5.5 for further evaluation.
This increases follow-up sensitivity and false alarms; it does not increase the
significance of a given peak. The two detrenders share the same photons.

The nominal single-sector period cap is 20 days, but the underlying method also
requires P<=baseline/2. A typical 27-day observation therefore searches only to
approximately 13.5 days. A separate **blind pooled search of every repeated star**
addresses signals below individual-sector thresholds and extends the nominal
period range to 40 days where the combined baseline permits. It is not restricted
to individually retained peaks. Per-sector tests at the pooled fitted ephemeris
are conditional consistency checks, not independent held-out confirmations.

Known transits are also masked for pooled residual searches when usable
ephemerides exist. Alternative detrending, individual events, raw SAP flux,
aperture changes, odd/even and secondary features, variability, catalogues,
external sectors and pixel localization test stronger leads.

## Measured sensitivity and false-alarm limitations

[Recovery controls](candidates/multisector_controls/README.md) contain 288 distinct
actual-noise injections, each evaluated by both old and new methods. In the
SPOC background, retention improves from 10/144 to 33/144; fixing the same S/N>=8
criterion changes 10 to 12, showing that most improvement comes from looser
selection. No 300-ppm injections survive these tested backgrounds. These two
backgrounds and finite seeded phases do not estimate population completeness.

Three genuine shallow planets recover at the expected period and epoch in all
nine control observations. The reported TIC 4206066 signal is retained in Sectors
32 and 98 by the expanded screen. Its known ephemeris across three sectors
provides a separate conditional accumulation control, not a discovery.

[Thirty-two block-surrogate searches](candidates/multisector_controls/NULLS.md)
from eight preselected stars produce seven retained cases, concentrated in two
stars. Astrophysical fragments, block joins and nonstationarity prevent a formal
null/FAP interpretation. No tiny p-value, global false-alarm rate or planetary
probability follows from this exercise.

## Data-product limitations and recognition

Hundreds of Sector 95 files contain good-quality finite flux but no positive
reported DET_FLUX_ERR. They are explicitly insufficient under the original
error-based screen. A **separate versioned alternative** estimates constant
per-observation uncertainty from robust adjacent flux differences. Original
FITS errors and flags stay unchanged; assumption labels and the failed original
search remain recorded. These alternative searches must be counted separately
and deduplicated when computing searched stars.

Fresh NASA confirmed/TOI and ExoFOP CTOI snapshots, published candidates,
RAVEN release tables and the supplied TIC 4206066 paper are used for recognition.
Controls distinguish known periods, harmonics, known hosts at different periods,
nearby contaminants and finite catalogue nonmatches. A nonmatch is not novelty.
Live source/alias checks supplement the tables for leading signals.

## Preservation and resumption

Raw products, prepared curves and SQLite checkpoints reside under
`/workspace/tess-research-data/multisector/`. External cloud-storage durability
is unverified. Exact raw query snapshots and metadata source hashes are tracked;
a compressed consistent SQLite backup restores completed searches without
rerunning flux analysis. The initial backup restored exactly and verified every
associated raw observation. Original snapshots are preserved when methods change.

Use the prepared Python environment and keep OMP/BLAS threads at one per process.
The active command is `python src/multisector_survey.py run --workers 3 --downloads 12`.
Completed rows are skipped; do not launch a second writer. Scientific method
hashes prohibit mixing configurations. The pooled and error-recovery stages have
separate checkpoints. See [PROGRESS.md](PROGRESS.md) for the current exact commands.
No paid models, external messages, submissions or observing requests were used.

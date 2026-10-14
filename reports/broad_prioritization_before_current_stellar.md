# Broad-survey prioritization

Snapshot: 2026-10-13T21:41:32.046020+00:00. The SQLite source was opened in read-only/query-only mode
with a consistent read transaction; the running survey was not changed.
This report and `results/broad/shortlist.csv` can be regenerated as the survey
grows. These are screening leads, not credible candidates or discoveries.

Source: `/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite`. Actual checkpoint statuses: `{'insufficient_data': 180, 'screened': 58789}`.
Examined 117,752 stored peaks. Before 1% within-star alias deduplication:

- unrecognized_lead: 479
- stellar_size_or_depth_warning: 776
- known_eb_control: 70
- wide_duty_warning: 438
- odd_even_warning: 102
- known_period_control: 39

Saved at most 200 rows per category, with exact JSON evidence
references, raw product URLs/paths/hashes, and source method hashes. Category
counts above are before deduplication/capping; CSV ranks are after these steps.

| TIC | Period d | Depth ppm | Local S/N | Events | Rp proxy Earth radii | Flags |
|---|---:|---:|---:|---:|---:|---|
| 389350449 | 1.271931 | 1547 | 15.3 | 12 | 3.53 |  |
| 358458645 | 0.840003 | 2082 | 17.5 | 20 | 8.80 |  |
| 281703963 | 1.229090 | 287 | 9.7 | 15 | 7.76 |  |
| 231705838 | 0.779829 | 1379 | 17.1 | 22 | 9.00 |  |
| 273863946 | 1.269869 | 2361 | 23.7 | 13 | 8.35 |  |
| 29478453 | 0.910703 | 3834 | 17.8 | 13 | 4.94 |  |
| 229810652 | 0.893387 | 1444 | 11.6 | 23 | 4.93 |  |
| 358462105 | 0.788781 | 2614 | 14.8 | 21 | 9.03 |  |
| 355765534 | 1.040316 | 2832 | 12.8 | 18 | unknown | unknown stellar radius |
| 278831310 | 1.685940 | 3243 | 11.5 | 10 | 6.25 |  |


Required support: initial screening guard pass; >=3 sampled events; >=3
positive-depth events and >=70% positive event fraction; local diagonal S/N
>=8.0; no event supplying >65% of positive event squared S/N. Event
sampling/sideband rules come from `fast_screen.event_stats`. Odd/even
differences >3.5 formal diagonal sigma are retained in a warning category.
Durations occupying >8% of a period are retained in a wide-duty warning
category, because broad variable-star troughs dominated the initial rankings.
This heuristic can demote real short-period planets and is not a veto.
This does not account for red noise, chosen-period effects, or global trials.

Crossmatching is recomputed against all prepared support catalogs. Exact or
harmonic matches are recovered controls. EBs are a separate control group;
period-matched nearby sources are contamination warnings. Known hosts at
different periods remain eligible for follow-up. Coarse 1% period matching is
a screening approximation, not ephemeris confirmation.

Rplanet proxy = Rstar * sqrt(max(BLS depth, local depth)), without a dilution
correction or transit fit. Negative sentinel radius/temperature/magnitude
values are invalid, never physical values. Proxy radii >22 Earth radii or depth
>2% are warnings, except deep cases around measured Rstar<=0.8 solar radii
with proxy radius<=22 Earth radii and depth<=30%, which are retained. These
limits prioritize follow-up and do not establish a physical planet boundary.
They can miss diluted EBs or planets around badly characterized stars.

If a measured stellar mass is available, compare duration with the circular
central-transit maximum and flag durations >1.5 times that prediction. QLP
survey metadata currently generally omit mass, so missing-mass cases are
explicitly flagged; no inferred stellar density is silently substituted.
Eccentricity, impact parameter, exposure integration, and stellar parameter
errors prevent this proxy from being a validation test.

The priority score saturates S/N at 25 and rewards recurring positive events
with distributed power, then demotes physical/odd-even/contamination warning
cases. It also reduces the rank of duty fractions >4% and depths >3,000 ppm
without discarding these cases; this addresses amplitude/variability-dominated
follow-up queues. The score and diagonal S/N have no false-alarm-probability meaning.

Next: load each leading raw QLP product; inspect raw/alternative apertures,
refine periods with an independent detrending choice and event-level residuals,
test variability and instrumental/common-mode timing, compare all sectors and
independent SPOC/TGLC data, then perform pixel/nearby-star and literature
checks. A strong folded BLS peak alone is insufficient.

Reproduce: `python src/broad_prioritize.py --db /workspace/tess-research-data/broad/qlp_s0001/survey.sqlite`.
Rejection reason counts (nonexclusive): `{'local_diagonal_snr_below_threshold': 109553, 'initial_screen_guard_failed': 33663, 'fewer_than_three_sampled_events': 37535, 'insufficient_positive_event_support': 40862, 'one_event_dominates_or_unavailable': 29420}`.

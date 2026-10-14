# Broad-survey prioritization

Snapshot: 2026-10-14T05:29:34.575709+00:00. The SQLite source was opened in read-only/query-only mode
with a consistent read transaction; the running survey was not changed.
This report and `results/broad/shortlist.csv` can be regenerated as the survey
grows. These are screening leads, not credible candidates or discoveries.

Current stellar table: `data/broad_catalogs/current_stars_20261014/current_stellar.csv`; validated SHA-256 `d0810cee943f505838f729e93688cc82ebf82024c704cacd6ad5c61017032533`. Exact current-TIC selection contains 3,410 targets. Missing current values do not fall back to old QLP headers.
Prioritizer method SHA-256: `7889c35b8d7601ce1dc201ad3a6e9f2553240c0708a27d2676ba950f0996f061`.

Source: `/workspace/tess-research-data/broad/qlp_s0001/survey.sqlite`. Actual checkpoint statuses: `{'insufficient_data': 495, 'screened': 159505}`.
Examined 319,496 stored peaks. Before 1% within-star alias deduplication:

- stellar_size_or_depth_warning: 3,209
- wide_duty_warning: 738
- known_period_control: 115
- unrecognized_lead: 427
- evolved_star_warning: 142
- known_eb_control: 193
- odd_even_warning: 140
- nearby_known_period_warning: 141

Saved at most 1000 rows per category, with exact JSON evidence
references, raw product URLs/paths/hashes, and source method hashes. Category
counts above are before deduplication/capping; CSV ranks are after these steps.

| TIC | Period d | Depth ppm | Local S/N | Events | Rp proxy Earth radii | Flags |
|---|---:|---:|---:|---:|---:|---|
| 369707296 | 1.309597 | 3471 | 22.6 | 12 | 15.50 |  |
| 260820871 | 1.345744 | 475 | 13.8 | 13 | 5.97 |  |
| 92232994 | 1.109926 | 3187 | 15.0 | 14 | 8.02 |  |
| 79449567 | 1.108486 | 1240 | 13.1 | 15 | 5.49 |  |
| 237906138 | 1.055504 | 1563 | 13.3 | 15 | 4.53 |  |
| 358458645 | 0.840003 | 2082 | 17.5 | 20 | 7.69 |  |
| 263077994 | 1.489532 | 3041 | 15.1 | 11 | 13.83 |  |
| 231705838 | 0.779829 | 1379 | 17.1 | 22 | 8.84 |  |
| 278779647 | 1.000668 | 1765 | 12.0 | 17 | 8.92 |  |
| 404972821 | 1.334864 | 2019 | 9.7 | 13 | 13.92 |  |


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
The current nearby-object radius is 120.0 arcsec (host
association still requires exact TIC or <=2 arcsec). This captures known
nearby periods beyond the former 30-arcsec threshold without assigning the
signal to a neighboring source or discarding the target.

Rplanet proxy = Rstar * sqrt(max(BLS depth, local depth)), without a dilution
correction or transit fit. Negative sentinel radius/temperature/magnitude
values are invalid, never physical values. Proxy radii >22 Earth radii or depth
>2% are warnings, except deep cases around measured Rstar<=0.8 solar radii
with proxy radius<=22 Earth radii and depth<=30%, which are retained. These
limits prioritize follow-up and do not establish a physical planet boundary.
They can miss diluted EBs or planets around badly characterized stars.
Current TIC giant classification or radius>3 solar radii receives an evolved
star warning rather than the isolated-small-star interpretation. Planetary
signals can exist in evolved systems; these warnings preserve the evidence
but require a physically appropriate model and stellar variability testing.
Old QLP radius/temperature/magnitude columns remain in the table for comparison.
When current mass and radius imply semimajor axis <= stellar radius, the
catalog-host orbital model receives an explicit
`semimajor_axis_inside_catalog_photosphere` warning and physical demotion.
Catalog mass/radius uncertainties are retained; the proxy does not identify
the source of a blend or establish a universal exclusion independent of
stellar-parameter errors. Missing current mass never invents this test.

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

Reproduce: `python src/broad_prioritize.py --db /workspace/tess-research-data/broad/qlp_s0001/survey.sqlite --catalog-radius-arcsec 120.0 --stellar-table data/broad_catalogs/current_stars_20261014/current_stellar.csv`.
Rejection reason counts (nonexclusive): `{'local_diagonal_snr_below_threshold': 297320, 'initial_screen_guard_failed': 90904, 'fewer_than_three_sampled_events': 101932, 'insufficient_positive_event_support': 110914, 'one_event_dominates_or_unavailable': 79745}`.

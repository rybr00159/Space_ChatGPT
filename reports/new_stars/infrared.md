# Infrared motion reconnaissance and failed stellar-novelty leads

The fixed experiment queried twelve disjoint 5-degree-radius fields at Galactic
latitude±60 degrees. Ten succeeded, returning 1946 quality/motion-selected CatWISE
sources across approximately 785 square degrees. Two timed out and are not counted
as searched. These are already catalogue sources; only their physical nature
could be newly characterized.

Gaia checks cover all 1946:1696 have nearby entries, 250 do not. The 250 unmatched
cases plus ten controls were checked against 2MASS at actual observation epochs,
using CatWISE motion positions at their documented MJD 57170 reference.27 unmatched
cases have reasonably consistent independent old positions, but **all 27 already
have TIC STAR classifications**. Several are explicitly known brown dwarfs.
No newly recognized star follows from absence in Gaia or SIMBAD.

Twenty-eight no-2MASS optical follow-ups gave one apparent proper-motion match,
CatWISE J 124349.41+005812.0. Two additional SDSS epochs falsify the proposed
optical confirmation: the counterpart is classified a galaxy and shifts only
0.126 arcsec from MJD 51258.33 to 53766.49, over 6.87 years. CatWISE motion predicts
2.64 arcsec. Thus a near-coincidental PS1 position must not be advertised as
independent motion confirmation. This does not exclude a different faint mover
within the blended WISE beam. Legacy catalogue requests returned 502 and are
retained as access failures, not zero matches.

Methods and raw query parameters are in `src/infrared_star_scout.py`. Stages
`fetch`, `gaia`, `counterparts`, `novelty`, and `optical` reuse hash-verified
responses. Data tables and failure manifests are in `data/new_stars/infrared/`;
the SDSS falsification is `data/new_stars/infrared_followup/optical_motion_rejection.json`.
All stages were rerun from the retained cache after adding explicit integrity
checks. The published CatWISE motion fit alone is not stellar identification,
mass classification or first recognition.

# Image-level source search and falsification

The fixed 12-field bright pilot produced 84 recurring peaks;32 were missing
from Gaia but in deeper catalogues. Apparent residual sources were affected
by saturation. A fixed 64-field expansion selected 16 primaries in each of
G 10–11, 11–12, 12–13.5, 13.5–15 by SHA 256(source_id), within the documented
sky/proper-motion/parallax cuts. All 64 were processed successfully.

The expansion retained 951 original/residual two-band peaks:463 deeper-catalogue
matches, 83 Gaia matches, 206 nonstellar morphologies, 110 halo-noise failures,
61 saturation trails, 26 incoherent/broad profiles, one nonrecurring residual
and one initially compact candidate. Native maskbits show that final peak
is saturated in every g/r/z band; its measured morphology cannot validate a star.
This final mask check is preserved separately rather than silently rewriting
the earlier classifier table.

The remaining pilot peak has no established independent detection at a
stationary or the primary's common-motion position. Exact Legacy exposure
epochs 2017.24–2018.20 and PS1 means 2011.78/2012.66/2012.07 were used. A
stationary-only comparison would have been inadequate for a primary moving
1.37 arcsec/year. The peak is also off the primary's own earlier/later track.
Its broadened profile and absence of reliable independent confirmation leave
extension, variability and artifacts unresolved.

3832 static flat-AB Gaussian injections test only matched-filter recovery.
At 4 arcsec, rates for magnitudes 17/19/21/23 are 61.4/19.1/2.85/0.41%; at 16 arcsec
they are 99.6/99.6/98.7/90.0%. These omit true PSFs, source motion, color,
subtraction and artifact classification. The compact-band gate accepts 200/277
detected catalogue-PSF controls. There is no calibrated full-search completeness
or false-alarm probability.

Raw products and URLs are retained outside Git. Earlier method files and all
numerical tables remain tracked; all bulk diagnostic images are preserved via
their hash manifest. [Curated evidence](../../candidates/new_stars/imaging/README.md)
links the final tested cases. Further useful work requires independent-night
or higher-resolution imagery, not treating Gaia absence as source novelty.

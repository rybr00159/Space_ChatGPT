# TIC 214311908: resolved ground-based transit brief and cost options

Prepared 17 October 2026. This is an observer-facing planning draft; no request, application, booking or purchase has been sent. Telescope suitability and the actual observing window require an observer's assessment.

## Scientific objective

Measure the target and nearby potential contaminating stars separately through a predicted eclipse. TESS pixels are about 21 arcseconds wide, blending the nominated dwarf with several sources. Ground images of adequate quality can distinguish the approximately 3.15-arcsecond neighbor. Here, **resolved** means spatially separating their photometry, not just taking more densely sampled measurements of the same unresolved blend.

The immediate question is: does TIC 214311908 itself undergo the approximately 1–2% loss of light, or does another star in its TESS aperture eclipse? One good complete event may settle this particular source question. Repeat observations strengthen it. A ground transit alone does not measure companion mass or exclude every unresolved binary/compact companion.

This is the kind of observation described by [TFOP SG1](https://tess.mit.edu/followup/): seeing-limited photometry to identify nearby eclipsing binaries contaminating TESS planet-candidate images. Describing a suitable SG1-type observation does not establish TFOP acceptance or entitlement to telescope time.

## Target and provisional timing

| Item | Planning value |
|---|---|
| Target | TIC 214311908, Gaia DR3 5925373272850784768 |
| Reference coordinates | Approximately RA 260.969699°, Dec −51.323833°; propagate proper motion for the observation epoch |
| Brightness | TIC T magnitude 13.01; red dwarf, with band-dependent optical brightness |
| Period | 5.691506 ± 0.000011 days, provisional adopted external-report value |
| Reference epoch | 2460840.9150 ± 0.0013 BJD_TDB |
| Expected event | Approximately 1.7 hours; approximately 1.6% of target-reference TESS flux if the nominated dwarf is the host |
| Principal close alternative | TIC 1515709327 / Gaia 5925373272827260544, approximately 3.15 arcseconds away in catalogue context |
| Other field sources | TIC 214311910 wide-companion pair at approximately 11.7 arcseconds; bright variable Mira TIC 214311906 at approximately 17.8 arcseconds |

Separations change with epoch. Do not schedule from these catalogue coordinates or convert BTJD/full BJD by guesswork. Verify an updated ephemeris and uncertainty, propagate covariance/systematic timing allowance, convert BJD_TDB to site-appropriate UTC, and check visibility, darkness, airmass, Moon and weather before requesting time. The declination favors southern observing sites such as Chile, South Africa and Australia. No exact future observing window is independently certified in this brief.

## Suggested observation

- A continuous series of calibrated images covering the approximately 1.7-hour event, plus roughly 1–1.5 hours before and after and appropriate timing margins: usually **about five hours**, with a longer reservation if needed.
- Aim for approximately 1–2-arcsecond image quality and sufficiently fine sampling, roughly 0.3–0.7 arcseconds/pixel as a planning target. These are not guarantees: the observer should test whether the target/close-neighbor contrast can actually be separated, using suitable apertures or PSF-fitting photometry. Deliberate defocusing would worsen this field's source-confusion problem.
- A stable red filter such as i, z or a suitable R/I band is a sensible first localization observation. Actual exposure times/cadence depend on telescope, throughput, saturation, and precision; aim for roughly minute-scale sampling with precision comfortably below the expected 1–2% target dip after comparison-star correction. An experienced observer should set the exposure plan.
- Measure light curves for the nominated target, the close faint neighbor and the wide companion pair, together with constant comparison stars. The bright Mira requires explicit attention to saturation, bleeding and PSF wings.
- Simultaneous multiband measurements are preferable for additional blend/colour tests, but a high-quality single-band resolved light curve can already provide important host identification. A two-band result must be modelled with limb darkening and dilution; equal raw depth in every band is not a universal planet criterion.
- Retain raw science/calibration frames, precise timestamps/time standard, WCS and seeing estimates, comparison-star identities, extraction settings, fluxes/errors/flags, and plots of the event and relevant controls.

An advanced amateur setup might succeed in good conditions, but aperture alone does not guarantee the necessary resolution and precision. A suitable approximately 0.5–1-metre facility is a sensible request; a larger telescope is not automatically necessary for the initial host-localization goal. A 2-metre simultaneous multiband instrument would provide stronger characterization.

## Costs and practical access

**It need not cost the project owner money.** A willing experienced observer, a research collaborator with allocated time, or an eligible funded allocation may contribute the observations without a commercial observing charge to the project owner. Availability, eligibility, observer effort and scientific priority must be agreed; this is not a promised free observation. It does not require owning a telescope or an existing ExoFOP candidate designation to approach possible collaborators.

The [AAVSO research-opportunities page](https://aavso.org/research-opportunities/) describes its role in connecting research requests with observers. This is a possible collaboration route, not an already approved campaign. A scientifically substantive collaborator may merit acknowledgment or coauthorship according to their contribution and an agreed publication plan; authorship is not something to buy.

Institutional access is programme-specific. For example, [LCO's current allocation rules](https://lco.global/observatory/proposal/process/) distinguish partner standard proposals, large Key Projects and discretionary time. For 2026B–2028B, the stated worldwide open-access support covers Key Projects and Director's Discretionary Time, not an unrestricted standard-proposal route. Discretionary time needs its stated urgency/science justification. A collaborator's appropriate existing allocation may be more practical than assuming a one-target application will be awarded. LCO also [sells network time](https://lco.global/observatory/sales/); no project-specific price was obtained.

Commercial remote access is an optional paid route. As an illustrative checked example, [iTelescope's standard-plan page](https://www.itelescope.net/pricing-standard) currently advertises approximately **US$19–118 per hour**, depending on telescope and plan, plus **US$39.95 or US$89.95 per 28-day membership**. Five chargeable hours therefore imply a cost in the hundreds of US dollars for many configurations. This is not a quote for a selected suitable telescope, and charged exposure/session time may differ from wall-clock reservation time. Confirm all of the following before spending: southern-site access, field resolution/sampling, filter/cadence control, continuous reservation, calibration data, billing basis, weather/repeat policy and any reduction/support charges. A low monthly membership price does not itself purchase a complete suitable transit observation.

The recommended first route is a scientific collaboration using this brief and the existing evidence dossier. Paid rental should follow an actual suitability check and quote. No observing fee or other payment has been incurred.

Supporting sources and exact snapshots are preserved under `data/ground_followup_options/20261017/`. The scientific assessment and uncertainties remain in the [comparison report](../candidates/tic214311908_external_comparison/README.md).

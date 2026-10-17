# TIC 214311908: publication and ExoFOP submission path

Checked 17 October 2026 against current official ExoFOP pages. This is guidance and a readiness assessment, not an external submission or a prediction of journal acceptance.

## Current rule

The [19 August 2026 announcement](https://exofop.ipac.caltech.edu/tess/news.php) says new community candidates must first be published in a peer-reviewed journal with online access, and approval to upload candidates is required. An ordinary ExoFOP account alone does not provide this approval.

The more detailed [candidate guidelines](https://exofop.ipac.caltech.edu/tess/candidate_help.php) add an exception:

> Candidates can also be presented in a Research Note of the American Astronomical Society (RNAAS) if the methodology to detect and vet the candidate has been previously published in a peer-reviewed journal that is available online and is cited by the Research Note.

RNAAS itself is [not peer reviewed](https://journals.aas.org/research-notes/). Its eligibility here is this explicit exception. Citing a standard BLS/library paper should not be assumed to establish that the full custom deblending and vetting methodology meets the condition. ExoFOP support can clarify a proposed RNAAS case; no query has been sent. GitHub or an arXiv preprint alone does not meet the stated journal-publication requirement.

## What the current evidence justifies

The evidence justifies preparing a carefully qualified **candidate paper**, seeking expert scientific review and approaching observers with a concrete follow-up case. Planet confirmation or a mass measurement is not universally required to publish a candidate analysis; journal editors and referees determine whether the contribution and evidence meet their standard.

The coherent six-sector approximately 5.6915-day signal, planet-compatible shape/size, improved PRF localization, source-specific alternative tests and reproducibility material provide a substantial starting point. They support the current comparative assessment, **strong, unvalidated planet candidate**. They do not establish a confirmed planet, a calibrated exact 99.9% probability or discovery priority.

I would complete these steps before journal submission:

1. **Make the central inference fully reproducible.** Archive the final deblended light curves with time system, flux scale, errors, flags, extraction settings, PRF/WCS/colour calibration and key injection outputs. The external report says some large intermediary products are absent from its repository; our saved raw FITS do not by themselves reproduce every external extraction. Supply a runnable end-to-end example with fixed versions and a durable data release.
2. **Present defensible error and false-positive accounting.** Show PRF/sky/noise sensitivity, close-source alternatives and their physical assumptions. Retain eccentric/compact/unresolved-host caveats where tested coverage is incomplete. Separate model-conditional FPP, finite simulation tails and empirical measurements. A candidate paper can state that validation is pending rather than rely on an exact small FPP.
3. **Complete the paper-level novelty and reference review.** Recheck candidate catalogues, aliases, literature/full-text and supplementary search reports near submission. Have human authors verify the central results, citations and claims; an exoplanet observer or transit-analysis collaborator would materially improve the review. Follow the selected journal's authorship and AI-assistance rules.

Resolved multiband photometry during a predicted event would be the most valuable new datum: determine whether the nominated dwarf or a nearby star dims. High-resolution companion imaging and reconnaissance spectroscopy/RVs address other remaining scenarios. These are strongly recommended scientific improvements; the public ExoFOP rule does not itself require all of them or an RV mass before a candidate can be published.

## Process once there is a qualifying publication

1. Log into the existing ExoFOP account and complete the [Published Candidate Upload Request](https://exofop.ipac.caltech.edu/tess/pub_candidate_upload_request.php), giving the published paper's URL and requested information. The form is login-gated, so its private fields were not inspected. Wait for approval before candidate creation.
2. Recheck the [TIC 214311908 page](https://exofop.ipac.caltech.edu/tess/target.php?id=214311908), current Community Candidates, TOIs and confirmed planets. Avoid duplicate or simple harmonic entries. No Community Candidate or TOI was shown on the public target page at this check; that must be refreshed later.
3. Use the approved candidate interface to choose the next available TIC-based candidate identifier, select discovery source **TESS**, and provide disposition **PC**, the paper URL, a descriptive initial note and a data tag. Let the form determine the available identifier; no new designation has been created here.
4. Enter period in days, full barycentric Julian transit epoch, depth in ppm and duration, with realistic uncertainties. Creation requires at least two positive entries among period/epoch/depth, but TESS TOI-working-group consideration requires **all four** of period, epoch, depth and duration to be positive. We recommend all four regardless. Initial creation parameters are fixed; later revisions are additional parameter sets, so check units/time standards carefully.
5. Upload supporting light curves, fitted-transit and false-positive diagnostics to the host TIC and associate them with the **same data tag** as the candidate. The [tag instructions](https://exofop.ipac.caltech.edu/tess/tag_help.php) explain this linkage. Candidate upload does not itself confer TOI status or planet confirmation.

The current data already contain the requested observables. The present obstacle to creating a new community candidate is the publication/approval requirement, not the absence of a measured ephemeris.

## Provisional values for preparation

These are examples from the supplied analysis's adopted solution, not final publication-approved values or an upload request. Reconcile and freeze the paper's own parameter set and errors before submission.

| Quantity | Provisional example |
|---|---|
| Period | 5.691506 ± 0.000011 days |
| Epoch | 2460840.9150 ± 0.0013 BJD_TDB; full JD, not BTJD |
| Depth | 15,800 ppm; quoted statistical uncertainty approximately 900 ppm, additional normalization/extraction systematic +1500/−1100 ppm |
| Duration | Approximately 1.70 hours; the quoted ±0.07-hour value is conditional on the fitted density/shape assumptions |
| Interpretation | Planet candidate; approximately five Earth radii if hosted by TIC 214311908 |

Official page snapshots and hashes are retained in `data/exofop_policy/20261017/`. Original science and the separate comparative evidence remain unchanged. No account credentials were accessed; no upload, approval request, email or observing request was sent.

## Journal shortlist checked 17 October 2026

**First choice to investigate: The Open Journal of Astrophysics (OJAp).** Its [official author guidelines](https://astro.theoj.org/for-authors) explicitly include extrasolar planets under astro-ph.EP, describe conventional editorial/referee review, and make submission, review and publication free to authors and readers. Its scientific-quality, originality, relevance and comprehensibility criteria still require a sound paper; zero fees are not evidence of easier acceptance. It uses an arXiv-overlay workflow. Posting first to arXiv is preferred, direct initial submission is possible, and the final version must be on arXiv for publication. The refereed journal publication, not the initial preprint, would support the ordinary ExoFOP route.

Its current policy explicitly permits generative AI assistance with disclosure in the acknowledgments. If AI is used for calculation, analysis or data visualization, the paper must describe the steps taken to test and validate it. This is relevant to this project's actual workflow; it is not permission to substitute an unverified model-generated report for a scientific manuscript. First-time arXiv submission may require endorsement.

**Most conventional specialist venue: The Astronomical Journal (AJ).** Its [scope](https://journals.aas.org/astronomical-journal/) emphasizes significant observational astronomical research, making an exoplanet/transit manuscript a natural subject fit. This does not establish that a single unvalidated candidate would be accepted. The [2026 publication charges](https://journals.aas.org/article-charges-and-copyright/) start at US$1,425 for at most 30 digital quanta, rise to US$3,162 for 31–50 and US$5,581 for 51–100, with additional long-article surcharges. These are publication charges, not payments that secure acceptance. The journal is fully open access. Check current charges and any assistance eligibility before committing to a route.

**Potential focused-observational alternative: JAAVSO.** The [AAVSO website](https://aavso.org/) identifies the Journal of the American Association of Variable Star Observers as its peer-reviewed journal. The amateur/professional observational community makes it worth an editorial inquiry about a compact, carefully qualified transit-candidate paper. However, its current application/author pages returned403 to this environment. Current charges and its acceptance of a single archival, unvalidated planet-candidate report were not verified. Do not label it the easiest or promise acceptance on that basis.

No reliable comparative acceptance-rate data were obtained. These are assessments of subject fit, cost and workflow, not measured acceptance probabilities. A legitimate peer-reviewed venue does not remove the need for reproducibility, careful false-positive treatment and a sufficiently useful scientific contribution.

Practical next step: prepare a concise abstract plus the folded signal, six-sector localization summary, principal alternative tests and limitations, then ask the selected editor whether an archival candidate paper without ground-based validation is suitable. No inquiry has been sent. A suggested question is:

> Would a research article reporting a coherent 5.6915-day, approximately five-Earth-radius TESS planet candidate around TIC 214311908 be within scope, based on six-sector photometry, calibrated source localization and explicit false-positive tests, while resolved ground-based validation remains pending? The manuscript would retain candidate status, disclose its AI-assisted workflow and provide reproducible analysis products.

Full manuscripts should be submitted to one journal at a time. RNAAS remains a separate, non-peer-reviewed possibility only under ExoFOP's explicit previously-published-methodology exception; it is not a general substitute for the refereed venues above. Source snapshots and access limitations are in `data/journal_options/20261017/`.

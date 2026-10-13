# TOI-4342 final joint gap-trimmed peak

P=19.454330 d, epoch BTJD 4107.418112, duration 5.51 h; screening S/N 9.37, supported events 2. Selected using all three sectors, so no independent heldout claim applies.

|Reduction|Width (d)|Local baseline|Depth (ppm)|Diagonal S/N|
|---|---|---|---|---|
|PDCSAP_FLUX|1|linear|1380|8.09|
|PDCSAP_FLUX|1|quadratic|1479|6.75|
|PDCSAP_FLUX|2|linear|1600|9.01|
|PDCSAP_FLUX|2|quadratic|1486|6.54|
|SAP_FLUX|1|linear|-340|-2.01|
|SAP_FLUX|1|quadratic|-72|-0.33|
|SAP_FLUX|2|linear|-338|-1.89|
|SAP_FLUX|2|quadratic|50|0.22|
|Native SAP|—|linear|-218|-1.37|
|Native SAP|—|quadratic|127|0.62|

Assessment: the corrected-flux dips lack a positive native-SAP counterpart (linear S/N about -1.37; quadratic about 0.62). No credible candidate is promoted. This diagnoses reduction dependence; it does not prove that the star lacks further planets.
These are local boxes with separate nuisance baselines, not calibrated false-alarm probabilities. Selection spans the full joint data and multiple grids. SAP/PDC share photons and local fits can still depend on baseline window. No pixel localization or physical transit inference is justified by weak or reduction-dependent support. Original boundary/alias diagnostics remain preserved separately.

Reproduce: `python src/vet_toi4342_joint.py`.

The two supported cycles have only 2 versus 31 ten-minute in-box samples and the same parity; odd/even agreement cannot be evaluated. Partial coverage is not evidence of two well-observed transits.

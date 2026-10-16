# Six-sector signal investigations

16 October 2026 UTC. **No credible new planet or previously unknown star is established.** [TIC229146531](tic229146531/README.md) is the strongest unresolved recurrent transit/eclipsing signal. It occurs in an already known astrometric binary, and grazing planet/star/blend alternatives remain unresolved.

The completed search covers5,058observations /4,391distinct stars in parts of Sectors83/86/95/96/98/105. Read the [campaign report](../../MULTISECTOR_RESEARCH_REPORT.md), [ranked assessments](../../results/multisector/ranked_assessments.csv), [deduplicated observation table](../../results/multisector/all_observations.csv) and [final batch report](../../reports/multisector/batch02_followup.md).

`followup_summary.csv` contains166diagnostic variants for66TICs /153distinct target-sector pairs, with overlapping periods, controls and assumptions. The separate final pooled queue adds14selected TICs /30observation tests; its table and full quantitative records are in `results/multisector_final_pool_followup/`. These counts overlap the primary survey and are not new stars or independent detections.

Individual `vetting.json` files preserve raw provenance, error assumptions, alternate baselines, individual event statistics, secondaries, parity, centroids and limitations. Their generic warning/status strings are screening outputs, not automatic astrophysical verdicts. Prefer the current dossier/ranked assessment over an early generated sector README. Frozen helper-source snapshots are in `method_snapshots/`.

TIC368248545is strongly contaminated by a known neighboring eclipsing binary. TIC182517223matches a published GaiaSB1. TIC353933809/183373155/71274048weaken under native-pixel/background tests. A weak Gaussian-only offset for TIC293268396was invalid and is not source attribution. Gap events, ill-conditioned Fourier baselines and a preserved zero-event helper exception do not count as physical nondetections.

All diagnostics were retained. BulkPNG/NPZfiles were moved intact outside ordinary Git with original/retained paths, byte counts and SHA256in [the preservation manifest](../../data/multisector/diagnostic_preservation.json). Important arrays and figures remain tracked, including the full TIC229146531dossier. Its six original QLPlight curves plus native pixel cube also have an exact, tested offline archive in `data/multisector/candidate_checkpoint/`. Public retrieval manifests preserve other originals; cloud filesystem persistence is unverified.

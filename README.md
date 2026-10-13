# Autonomous TESS candidate research

A skeptical search for overlooked transit signals in public TESS observations. A detected dip is not a confirmed planet or evidence of novelty. This project preserves negative findings and candidate falsification alongside promising results.

## Installation

Python 3.12; create a virtual environment and install `requirements.txt`. Current task uses `/workspace/openrouter-check-venv/bin/python`. Source material: `Reddit Post.zip`, retained unchanged.

## Initial investigation

Compare additional-planet searches in small known systems with less-studied stars. Start with a manageable nearby cool-star sample and positive controls, mask known ephemerides, and inspect residual signals. Catalog and data access checks and sample selection are in progress; no new candidate has yet been established.

## Preservation and resumption

Raw observations and intermediate arrays are cached outside ordinary Git under `/workspace/tess-research-data`. This is an overlay filesystem, not a separately verified durable volume. Codex cloud publication can snapshot prepared filesystem state; persistence across tasks is not assumed without publication. Recorded archive product IDs, URLs, processing metadata, and SHA-256 hashes support re-download. Small tables, diagnostic figures, code, and reports belong in Git. See `PROGRESS.md` for the exact checkpoint.

## Integrity and spending

Use catalog and literature checks before suggesting novelty. Transit S/N is not a survey false-alarm probability. Nearby eclipsing binaries, systematics, and detrending artifacts require investigation. No paid model requests are required by the pipeline; OpenRouter use during this research phase is currently $0 against the authorized $10 cap.

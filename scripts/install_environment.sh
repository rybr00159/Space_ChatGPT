#!/usr/bin/env bash
set -euo pipefail
cd /workspace/Space_ChatGPT
TESS_VENV=/workspace/openrouter-check-venv
if [[ ! -x "$TESS_VENV/bin/python" ]]; then
  python3 -m venv "$TESS_VENV"
fi
"$TESS_VENV/bin/python" -m pip install --cache-dir /workspace/.pip-cache -r requirements.txt
mkdir -p /workspace/tess-research-data/lightcurves /workspace/tess-research-data/intermediate /workspace/tess-research-data/matplotlib
export MPLCONFIGDIR=/workspace/tess-research-data/matplotlib
"$TESS_VENV/bin/python" -m pip check
"$TESS_VENV/bin/python" -c 'import astropy, astroquery, lightkurve, numpy, scipy, pandas; print("Astronomy environment import check passed")'

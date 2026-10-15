#!/usr/bin/env python3
"""Run the retained image method with bulk plots outside ordinary Git."""
from pathlib import Path
import os
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import uncatalogued_imaging as analysis

analysis.PLOTS=Path(os.environ.get('NEW_STAR_PLOT_DIR','/workspace/tess-research-data/new_stars/imaging/plot_archive'))
if __name__=='__main__':analysis.main()

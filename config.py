"""Central paths for every Python script. Edit here (or set the environment variables) - nowhere else.
Run every script from the repository root, e.g.  python media_PyDESeq2/media_K1_pydeseq2.py"""
import os
DATA_DIR = os.environ.get("CHO_DATA_DIR", "data")      # input files (counts, metadata, gene lists)
OUT_DIR = os.environ.get("CHO_OUT_DIR", "results")     # everything the scripts write
os.makedirs(OUT_DIR, exist_ok=True)

#!/usr/bin/env python3
"""Run the shared GENCODE v44 TSS ±2 kb promoter-BED generator."""

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
COMMON_GENERATOR = SCRIPT_DIR.parent / "common" / "make_promoter_bed_2kb_2kb.py"
GTF = "/share/LCZX_Data/ref/gencode.v44.basic.annotation.gtf"
OUTPUT = "/share/home/rzli/METHSCAN/Meth_diff/common_resources/promoter_annotation/gencode_v44_basic_promoter_2kb_up_2kb_down.bed"

subprocess.run(
    [sys.executable, str(COMMON_GENERATOR), "--gtf", GTF, "--output", OUTPUT],
    check=True,
)

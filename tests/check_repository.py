#!/usr/bin/env python3
"""Lightweight checks that do not require UK Biobank data or MATLAB."""

from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    "00_run_matlab_pipeline.m",
    "01_extract_ukb_tabular_data.R",
    "02_compute_sfc_dsfc_glasser.m",
    "03_prepare_sfc_analysis_tables.m",
    "04_run_demographic_age_associations.m",
    "05_run_phenotype_lifestyle_associations.m",
    "06_run_genetic_associations.m",
    "07_prepare_gwas_inputs.m",
    "07_run_gwas.sh",
    "08_run_mediation_analysis.R",
    "09_run_robustness_analyses.m",
    "config_example.m",
]

errors = []
for relative_path in REQUIRED:
    if not (ROOT / relative_path).is_file():
        errors.append(f"Missing required file: {relative_path}")

code_files = [
    path
    for suffix in ("*.m", "*.R", "*.sh")
    for path in ROOT.rglob(suffix)
    if ".git" not in path.parts
]
han_pattern = re.compile(r"[\u3400-\u9fff]")
private_path_pattern = re.compile(r"/public/home|/home1/|/Users/")
for path in code_files:
    text = path.read_text(encoding="utf-8")
    if han_pattern.search(text):
        errors.append(f"Non-English text found in code: {path.relative_to(ROOT)}")
    if private_path_pattern.search(text):
        errors.append(f"Private absolute path found in code: {path.relative_to(ROOT)}")

subprocess.run(["bash", "-n", str(ROOT / "07_run_gwas.sh")], check=True)

for path in sorted(ROOT.glob("*.m")) + sorted((ROOT / "utilities").glob("*.m")):
    stack = []
    for line_number, raw_line in enumerate(path.read_text().splitlines(), 1):
        line = raw_line.split("%", 1)[0].strip()
        opener = re.match(
            r"^(if|for|parfor|while|switch|try|function|classdef|methods|properties)\b",
            line,
        )
        if opener:
            stack.append((opener.group(1), line_number))
        if re.match(r"^end\b", line):
            if stack:
                stack.pop()
            else:
                errors.append(f"Unmatched end in {path.name}:{line_number}")
    if stack:
        errors.append(f"Unclosed MATLAB blocks in {path.name}: {stack}")

if errors:
    print("\n".join(errors), file=sys.stderr)
    raise SystemExit(1)
print("Repository checks passed.")

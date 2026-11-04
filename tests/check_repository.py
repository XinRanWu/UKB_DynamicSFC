#!/usr/bin/env python3
from pathlib import Path
import re, subprocess, sys

root = Path(__file__).resolve().parents[1]
required = [
    "README.md", "config/config.example.m", "config/config.example.R",
    "matlab/run_primary_pipeline.m", "matlab/01_compute_primary_coupling.m",
    "matlab/02_prepare_analysis_tables.m", "matlab/03_descriptive_and_age_associations.m",
    "matlab/04_phenotype_associations.m", "matlab/05_lifestyle_associations.m",
    "matlab/06_genetic_associations.m", "matlab/07_prepare_gwas_inputs.m",
    "R/03_age_gamm_visualization.R",
    "R/08_lifestyle_indirect_effects.R", "R/09_age_adprs_structural_equation_models.R",
    "R/10_socioeconomic_association_sensitivity.R", "R/11_formal_sfc_dsfc_age_comparison.R",
    "matlab/10_robustness_and_sensitivity.m", "matlab/10b_intrinsic_bold_adjustment.m", "shell/07_run_gwas.sh",
    "docs/ANALYSIS_MAP.md", "docs/DATA_CONTRACT.md"
]
errors = [f"Missing required file: {p}" for p in required if not (root/p).is_file()]
code = [p for pattern in ("*.m","*.R","*.sh") for p in root.rglob(pattern) if ".git" not in p.parts]
private = re.compile(r"/public/home|/home1/|/Users/xinranwu")
han = re.compile(r"[\u3400-\u9fff]")
for path in code:
    text = path.read_text(encoding="utf-8")
    if private.search(text): errors.append(f"Private absolute path: {path.relative_to(root)}")
    if han.search(text): errors.append(f"Non-English code text: {path.relative_to(root)}")

# Lightweight MATLAB structure check. This is not a replacement for MATLAB,
# but catches accidental missing end statements without protected UKB data.
openers = re.compile(r"^(if|for|parfor|while|switch|try|function|arguments)\b")
for path in root.rglob("*.m"):
    stack = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if line.startswith("%"):
            continue
        match = openers.match(line)
        if match: stack.append((match.group(1), number))
        closes = re.match(r"^end\b", line) or (match and re.search(r"[;,]\s*end\s*;?$", line))
        if closes:
            if stack: stack.pop()
            else: errors.append(f"Unmatched end: {path.relative_to(root)}:{number}")
    if stack: errors.append(f"Unclosed MATLAB blocks: {path.relative_to(root)} {stack}")
subprocess.run(["bash","-n",str(root/"shell/07_run_gwas.sh")],check=True)
if errors:
    print("\n".join(errors),file=sys.stderr); raise SystemExit(1)
print(f"Repository checks passed ({len(code)} code files).")

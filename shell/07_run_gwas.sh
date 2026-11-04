#!/usr/bin/env bash
set -euo pipefail

: "${PLINK2_BIN:?Set PLINK2_BIN to the PLINK 2 executable}"
: "${BFILE_PATTERN:?Set BFILE_PATTERN, for example /data/ukb_imp_chr%d}"
: "${PHENO_FILE:?Set PHENO_FILE to coupling_phenotypes.tsv}"
: "${COVAR_FILE:?Set COVAR_FILE to the analysis-ready genetic covariate file}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR for GWAS results}"

mkdir -p "$OUTPUT_DIR"
phenotypes=(SFC_VN SFC_SMN SFC_DAN SFC_VAN SFC_LN SFC_FPN SFC_DMN \
            DSFC_VN DSFC_SMN DSFC_DAN DSFC_VAN DSFC_LN DSFC_FPN DSFC_DMN)
covariates="Sex_0_0,AgeAttend_2_0,BMI_2_0,HeadMotion_2_0,SNR_2_0,Vol_WB_TIV_2_0,Centre_0_0,PC1-PC10"

for chr in $(seq 1 22); do
  bfile=$(printf "$BFILE_PATTERN" "$chr")
  for phenotype in "${phenotypes[@]}"; do
    "$PLINK2_BIN" \
      --bfile "$bfile" \
      --pheno "$PHENO_FILE" \
      --pheno-name "$phenotype" \
      --covar "$COVAR_FILE" \
      --covar-name "$covariates" \
      --glm hide-covar cols=+a1freq,+nobs \
      --out "$OUTPUT_DIR/chr${chr}_${phenotype}"
  done
done

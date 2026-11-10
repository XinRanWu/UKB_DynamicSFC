#!/usr/bin/env bash
set -euo pipefail

: "${PLINK2_BIN:?Set PLINK2_BIN to the PLINK 2 executable}"
: "${BFILE_PATTERN:?Set BFILE_PATTERN, for example /data/ukb_imp_chr%d}"
: "${PHENO_FILE:?Set PHENO_FILE to coupling_phenotypes.tsv}"
: "${COVAR_FILE:?Set COVAR_FILE to gwas_covariates.tsv}"
: "${KEEP_FILE:?Set KEEP_FILE to unrelated_imaging_keep.tsv}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR for GWAS results}"

mkdir -p "$OUTPUT_DIR"
phenotypes=(SFC_VN SFC_SMN SFC_DAN SFC_VAN SFC_LN SFC_FPN SFC_DMN \
            DSFC_VN DSFC_SMN DSFC_DAN DSFC_VAN DSFC_LN DSFC_FPN DSFC_DMN)
phenotype_names=$(IFS=,; echo "${phenotypes[*]}")
covariates="Sex_0_0,AgeAttend_2_0,AgeAttend_2_0_c2,BMI_2_0,HeadMotion_2_0,SNR_2_0,Vol_WB_TIV_2_0,ImagingCentre,PC1-PC20"
threads=${PLINK_THREADS:-10}

for chr in $(seq 1 22); do
  bfile=$(printf "$BFILE_PATTERN" "$chr")
  "$PLINK2_BIN" \
    --bfile "$bfile" \
    --keep "$KEEP_FILE" \
    --maf 0.01 \
    --geno 0.05 \
    --hwe 1e-6 \
    --pheno "$PHENO_FILE" \
    --pheno-name "$phenotype_names" \
    --covar "$COVAR_FILE" \
    --covar-name "$covariates" \
    --covar-variance-standardize \
    --glm hide-covar cols=+a1freq,+nobs \
    --threads "$threads" \
    --out "$OUTPUT_DIR/gwas_chr${chr}_7n"
done

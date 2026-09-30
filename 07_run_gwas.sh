#!/usr/bin/env bash
set -euo pipefail

# Run chromosome-wise PLINK 2 GWAS after 07_prepare_gwas_inputs.m.
# Required environment variables:
#   PLINK2_BIN, BFILE_PATTERN, PHENO_FILE, COVAR_FILE, OUTPUT_DIR
# BFILE_PATTERN must contain one printf placeholder, for example:
#   /approved/ukb/qc/ukb_imp_chr%d

: "${PLINK2_BIN:?Set PLINK2_BIN to the PLINK 2 executable}"
: "${BFILE_PATTERN:?Set BFILE_PATTERN with a chromosome placeholder}"
: "${PHENO_FILE:?Set PHENO_FILE to the phenotype TSV}"
: "${COVAR_FILE:?Set COVAR_FILE to the genetic covariate file}"
: "${OUTPUT_DIR:?Set OUTPUT_DIR for GWAS results}"

mkdir -p "${OUTPUT_DIR}"

PHENO_NAMES="${PHENO_NAMES:-SFC_Visual,SFC_Somatomotor,SFC_DorsalAttention,SFC_VentralAttention,SFC_Limbic,SFC_Frontoparietal,SFC_Default,SFC_Global,DSFC_Visual,DSFC_Somatomotor,DSFC_DorsalAttention,DSFC_VentralAttention,DSFC_Limbic,DSFC_Frontoparietal,DSFC_Default,DSFC_Global}"
COVAR_NAMES="${COVAR_NAMES:-Sex_0_0,AgeAttend_2_0,BMI_2_0,HeadMotion_2_0,Vol_WB_TIV_2_0,Centre_0_0,PC1-PC20}"

for chromosome in $(seq 1 22); do
  printf -v bfile "${BFILE_PATTERN}" "${chromosome}"
  "${PLINK2_BIN}" \
    --bfile "${bfile}" \
    --glm hide-covar \
    --pheno "${PHENO_FILE}" \
    --pheno-name "${PHENO_NAMES}" \
    --covar "${COVAR_FILE}" \
    --covar-name "${COVAR_NAMES}" \
    --covar-variance-standardize \
    --out "${OUTPUT_DIR}/gwas_chr${chromosome}"
done

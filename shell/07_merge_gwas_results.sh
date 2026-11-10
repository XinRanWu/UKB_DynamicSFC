#!/usr/bin/env bash
set -euo pipefail

: "${GWAS_CHR_DIR:?Set GWAS_CHR_DIR to the chromosome-level PLINK output directory}"
: "${GWAS_MERGED_DIR:?Set GWAS_MERGED_DIR for the merged phenotype files}"

mkdir -p "$GWAS_MERGED_DIR"
phenotypes=(SFC_VN SFC_SMN SFC_DAN SFC_VAN SFC_LN SFC_FPN SFC_DMN \
            DSFC_VN DSFC_SMN DSFC_DAN DSFC_VAN DSFC_LN DSFC_FPN DSFC_DMN)

for phenotype in "${phenotypes[@]}"; do
  merged="$GWAS_MERGED_DIR/gwas.${phenotype}.glm.linear"
  first="$GWAS_CHR_DIR/gwas_chr1_7n.${phenotype}.glm.linear"
  [[ -s "$first" ]] || { echo "Missing or empty file: $first" >&2; exit 1; }
  head -n 1 "$first" > "$merged"

  for chr in $(seq 1 22); do
    part="$GWAS_CHR_DIR/gwas_chr${chr}_7n.${phenotype}.glm.linear"
    [[ -s "$part" ]] || { echo "Missing or empty file: $part" >&2; exit 1; }
    tail -n +2 "$part" >> "$merged"
  done

  gzip -c "$merged" > "$merged.gz"
done

echo "Merged 22 chromosomes for ${#phenotypes[@]} phenotypes into $GWAS_MERGED_DIR"

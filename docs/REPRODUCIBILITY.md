# Reproducibility and provenance

The released workflow is aligned to the manuscript dated 4 November 2026. It replaces exploratory scripts that mixed data preparation, modeling, and plotting and that contained private cluster paths. Those scripts remain recoverable from Git history but are not part of the executable release.

## Fixed analytical decisions

- Primary FC is Pearson correlation within 30-TR windows advanced by 5 TRs.
- Regional window-wise SFC is Spearman correlation between SC and FC profiles after removing the self-edge.
- The primary whole-profile definition retains zero-valued SC entries and signed FC values.
- Static SFC is the temporal mean and DSFC is the temporal standard deviation of window-wise regional SFC.
- Participants with mean framewise displacement above 0.20 mm are excluded from association analyses.
- Continuous outcomes and continuous predictors are standardized within each model sample. Binary predictors remain 0/1.
- Assessment centre is a random intercept in the primary MATLAB LMMs.
- Statistical indirect effects are estimated in lavaan with covariates explicitly included in the mediator and outcome equations. No mediator residualization is used.

## Sensitivity definitions

`matlab/10_robustness_and_sensitivity.m` recomputes the metric for 40- and 50-TR windows, mean-FA structural weights, GSR-preprocessed fMRI, positive-only FC, non-zero SC, and both edge restrictions. It also calculates parcel-level AC1 and BOLD variability. Within-participant map correspondence uses Spearman correlation; corresponding network values across participants use Pearson correlation; group maps use Spearman correlation. Spin-test inference requires the same surface projection and 10,000 rotation indices used for the manuscript and should be archived with the release metadata.

## Software

The manuscript reports MATLAB R2025a, R 4.2.3, and PLINK 2.0. R package requirements are `data.table`, `lavaan`, `lme4`, `lmerTest`, `emmeans`, `dplyr`, `tidyr`, and `mgcv`. Run `sessionInfo()` after each final R analysis and archive the text with the corresponding outputs.

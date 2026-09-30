# Reproducibility decisions requiring author confirmation

These items were found while aligning the public code with the manuscript dated 30 September 2026. They are scientific or reporting decisions, so the cleanup does not silently guess the intended answer.

## 1. Streamline normalization

The current Methods state that no additional streamline weighting or filtering was applied. The original MATLAB calculation in Git history divides every streamline-count edge by the product of the two parcel surface sizes. These procedures are not equivalent and can change SFC estimates.

The cleaned `config_example.m` defaults to `NORMALIZE_STREAMLINES_BY_PARCEL_SIZE = false` to follow the manuscript. If the reported results were produced with parcel-size normalization, update the Methods and set the flag to true with a documented `roi_size` file.

## 2. MATLAB version

The computation notes in the repository identify MATLAB R2018b, while the manuscript Code availability section identifies MATLAB 2025a. Confirm the version used for the reported results and state it consistently. A second version may be listed as a tested compatibility version after the pipeline is rerun.

## 3. Number of PRS measures

The Results text refers to 87 PRSs, whereas the Methods describe 90 Standard and Enhanced PRSs. The committed `UKB_PRS_info.xlsx` and the final model output should be used to report the exact number tested after exclusions.

## 4. Multiple-comparison correction

The mixed-model scripts return raw P values. The manuscript describes family-wise correction over different test families, including 7 x phenotype and 360-region families. The exact family definitions and corrected columns should be added to the final result-export step rather than inferred from plotting code.

## 5. Mediation covariate residualization

The manuscript says that fixed and random effects were removed before mediation. The original script residualizes the mediator with a mixed model and enters age and APOE dosage directly into the structural-equation model. The cleaned script preserves that executable behavior and labels the analysis as a statistical indirect-effect model. Confirm whether the outcome was also residualized in the reported analysis.

## 6. Schaefer-to-Glasser spatial comparison

The original robustness script in Git history uses local surface-projection and spin-permutation resources that are not present in the repository. Window-length and FA comparisons are reproducible from configured MAT files. To reproduce the cross-atlas spatial statistic, add the exact atlas-to-surface resources, source/version citations, and spin-permutation file or replace the analysis with a fully documented public package.

## 7. Upstream imaging preparation

The repository starts from participant-level cell arrays derived from UKB Connectome releases. The exact release identifiers, download dates, matrix orientation, and any reorganization commands should be recorded in a release manifest. This is necessary because participant-level UKB inputs cannot be redistributed publicly.

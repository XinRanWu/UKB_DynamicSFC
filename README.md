# UKB_DynamicSFC
This repository contains the data analysis code used for the paper entitled 
“Stable brain architecture but flexible dynamics: distinct structure–function 
coupling signatures link aging and genetic Alzheimer's risk to cognitive decline 
in 34,067 adults.”


## Data provenance

The analyses in this study were conducted using UK Biobank data obtained
under the legacy UKB data-access framework, before the transition to the
current UK Biobank Research Analysis Platform (RAP).

The imaging inputs were obtained from data releases of the UKB Connectome
project. Resting-state fMRI regional time series and diffusion MRI structural
connectivity matrices were extracted from the corresponding bulk releases
and reorganized into MATLAB-compatible participant-level data structures.

For each modality, participant-level imaging data were stored as an `N × 1`
MATLAB cell array, where each cell contained the imaging matrix or regional
time series for one participant. A corresponding vector of UK Biobank
participant identifiers (`eid`) was stored alongside each cell array.

For example, the principal Glasser-360 inputs used by the dynamic
structure-function coupling analysis were:

    fMRI_Glasser_TS_2_0
        N × 1 cell array containing participant-level regional
        resting-state fMRI time series.

    fMRI_Glasser_eID_2_0
        N × 1 vector containing the corresponding UK Biobank participant IDs.

    connectome_streamline_count_10M_2_0
        M × 1 cell array containing participant-level structural
        connectivity matrices based on streamline counts.

    dMRI_Glasser_S1_eID_2_0
        M × 1 vector containing the corresponding UK Biobank participant IDs.

Participants present in both imaging modalities were matched using `eid`
before calculation of structure-function coupling.

Phenotypic, demographic, cognitive, lifestyle, health, and imaging-quality
variables were extracted from authorized UK Biobank bulk tabular files.
The script `01_extract_ukb_tabular_data.R` documents the study-specific UKB
field IDs and reproduces the corresponding tabular extraction step.

Participant-level UK Biobank data are not redistributed with this repository.
Researchers wishing to reproduce the analysis must obtain the corresponding
data through their own approved UK Biobank project.

Because the original study was performed before migration to the UK Biobank
RAP, the supplied code reflects the original local/HPC-based data structure
rather than a RAP-native data-access workflow. Researchers working within the
current RAP environment should retrieve the corresponding UKB variables and
UKB Connectome-derived imaging data within their approved project and map
them to the documented input structures.

## Dynamic SFC computation

Participant-level dynamic functional connectivity was calculated from the
Glasser-360 resting-state fMRI time series using a sliding-window procedure.

The primary analysis used a 30-TR window with a step size of 5 TR. Alternative
window lengths of 40 and 50 TR were additionally evaluated as sensitivity
analyses.

Structural connectivity matrices were restricted to the 360 cortical Glasser
parcels and streamline counts were normalized according to parcel-pair surface
size.

Dynamic structure-function coupling was calculated separately for each
participant and sliding window. The temporal mean of structure-function
coupling was used as the average SFC measure, whereas the temporal standard
deviation was used as the DSFC measure.

The relevant script is:

    02_compute_sfc_dsfc_glasser.m

An example SLURM submission script is provided as:

    submit_sfc_dsfc.slurm

The original calculations were performed using MATLAB R2018b on a Linux HPC
cluster. MATLAB Parallel Computing Toolbox was used for participant-level
parallelization.

The supplied SLURM file documents the original computational workflow.
Cluster-specific settings such as partition name, memory allocation, module
names, and filesystem locations should be adapted to the user's computing
environment.

Code for

Weed L, Jamgochian A, St. Hilaire MA, Cheng P, Kochenderfer MJ, Zeitzer JM. Circadian Phase Estimation From Ambulatory Wearables With Particle Filtering: Accuracy Depends on Initialization, Recording Duration, and Light Exposure. *Journal of Biological Rhythms*. 2025;0(0). doi:10.1177/07487304251392289
  

### Overview
This repository contains scripts to preprocess wearable data, compute light diet characteristics, run circadian light-response models, and evaluate estimation performance. The workflow is primarily MATLAB, with supporting R code for FPCA.

### Repository layout
- `a_OrganizeData/`: raw-to-tidy organization
  - `organizeData.m`
- `b_PreprocessData/`: handling missingness and basic preprocessing
  - `findMissing.m`, `medianImputation.m`, `preprocess.m`
- `c_SubsetData/`: inclusion/exclusion and cohort subsetting
  - `includedSWsubjects.m`, `subsetData.m`
- `d_LightDietCharacteristics/`: light-features and FPCA
  - `calcIVIS.m`, `LightDietTable_betweenDatasets.m`, `fpca_code.r`, `func_fpca.r`, `fpcaPlots.m`
- `e_modeling/`: circadian model integrations and batch runs
  - `processL_stHilaire2007.m`, `processPODE.m`, `runLightModel_Wearable_ODE15_AMP.m`, `submitAll.sh`
- `f_Performance/`: accuracy metrics, distributions, correlations
  - `calcIVIS.m`, `DLMODistributions.m`, `estimateVonMisesParamsWithEntropy.m`, `examineInitDays.m`, `examineInits_14Days.m`, `extractTinitDays.m`, `lightDietCorrelations.m`, `LinsBias.m`

### Citation
```bibtex

@article{Weed2026_Circadian,
   author = {Weed,  Lara and Jamgochian,  Arec and St. Hilaire,  Melissa A. and Cheng,  Philip and Kochenderfer,  Mykel J. and Zeitzer,  Jamie M.},
   title = {Circadian Phase Estimation From Ambulatory Wearables With Particle Filtering: Accuracy Depends on Initialization, Recording Duration, and Light Exposure},
   journal = {Journal of Biological Rhythms},
   volume = {41},
   number = {1},
   pages = {42-52},
   ISSN = {0748-7304},
   DOI = {10.1177/07487304251392289},
   url = {https://doi.org/10.1177/07487304251392289},
   year = {2026},
   type = {Journal Article}
}



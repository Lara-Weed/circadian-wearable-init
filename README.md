Code for

Weed L, Jamgochian A, St. Hilaire MA, Cheng P, Kochenderfer MJ, Zeitzer JM. Circadian phase estimation from ambulatory wearables with particle filtering: accuracy depends on initialization, recording duration, and light exposure. *Journal of Biological Rhythms*. In press, 2025.

### Overview
This repository contains scripts to preprocess wearable data, subset cohorts, compute light diet characteristics, run circadian light-response models, and evaluate estimation performance. The workflow is primarily MATLAB, with supporting R code for FPCA.

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
If you use this code, please cite the paper above. Placeholder BibTeX:
```bibtex
@article{Weed2025Circadian,
  title = {Circadian phase estimation from ambulatory wearables with particle filtering: accuracy depends on initialization, recording duration, and light exposure},
  author = {Weed, Lara and Jamgochian, Arec and St. Hilaire, Melissa A. and Cheng, Philip and Kochenderfer, Mykel J. and Zeitzer, Jamie M.},
  year = {2025},
  journal = {Journal of Biological Rhythms},
  notes = {In Press as of 7 Oct 2025}
}


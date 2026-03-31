# Repeated hold-out validation versus k-fold cross validation

## Title ##
Repeated hold-out validation splits can be more truthful than k-fold cross-validation splits in supervised machine learning

## Description ##
This repository contains R scripts that execute experiments to compare the results of repeated hold-out validation and the results of k-fold cross validation.
For five different datasets, the script here reads the target feature and applies repeated hold-out validation and k-fold cross-validation several times.
For each dataset, the script compares the ground truths of the test sets of the repeated hold-out with the global dataset targets, and does the same with the k-fold cross validation test set targets, through DTW.

## Datasets availability ##

The EHR datasets employed in this study are publicly available in the supplementary materials of the corresponding original publications under the *Creative Commons Attribution 4.0 International (CC BY 4.0)* license. Direct access URLs are provided below:

- **Cardiac Arrest**:
  https://figshare.com/articles/dataset/Mortality_after_out-of-hospital_cardiac_arrest_in_a_Spanish_Region/4876247?file=8166893

- **ColonRectal Cancer**:
  https://figshare.com/articles/dataset/The_effect_of_epidural_analgesia_on_cancer_progression_in_patients_with_stage_IV_colorectal_cancer_after_primary_tumor_resection_A_retrospective_cohort_study/6846365?file=12464069

- **Type 1 Diabetes**:
  https://figshare.com/articles/dataset/Circulating_osteocalcin_as_a_bone-derived_hormone_is_inversely_correlated_with_body_fat_in_patients_with_type_1_diabetes/8079389?file=15057092

- **Heart Failure & Depression**:
  https://figshare.com/articles/dataset/Comorbid_Depression_and_Heart_Failure_A_Community_Cohort_Study/3916224?file=6130425

- **Neuroblastoma**:
  https://doi.org/10.7717/peerj.5665/supp-5

## Datasets citations ##
The real datasets are derived from electronic health records and include:

- **Neuroblastoma**
  Data derived from the cohort described in:
  Ma Y., Zheng J., Feng J., Chen L., Dong K., Xiao X.
  *Neuroblastomas in eastern China: a retrospective series study of 275 cases in a regional center.*
  PeerJ, 6:e5665, 2018.
  DOI: https://doi.org/10.7717/peerj.5665

- **Type one Diabetes**
  Data informed by the following studies:
  Takashi Y., Ishizu M., Mori H., Miyashita K., Sakamoto F., Katakami N., et al.
  *Circulating osteocalcin as a bone-derived hormone is inversely correlated with body fat in patients with type 1 diabetes.*
  PLOS One, 14(5):e0216416, 2019.
  DOI: https://doi.org/10.1371/journal.pone.0216416

- **Sepsis & Systemic Inflammatory Response Syndrome (SIRS)**
  Data derived from the following sources:
  Gucyetmez B., Atalan H.K.
  *C-reactive protein and hemogram parameters for the non-sepsis systemic inflammatory response syndrome and sepsis: what do they mean?*
  PLOS One, 11(2):e0148699, 2016.
  DOI: https://doi.org/10.1371/journal.pone.0148699

- **Heart Failure & Depression**
  Data based on the cohort described in:
  Jani B.D., Mair F.S., Roger V.L., Weston S.A., Jiang R., Chamberlain A.M.
  *Comorbid depression and heart failure: a community cohort study.*
  PLOS One, 11(6):e0158570, 2016.
  DOI: https://doi.org/10.1371/journal.pone.0158570

- **Cardiac Arrest**
  Data derived from:
  Requena-Morales R., Palazón-Bru A., Rizo-Baeza M.M., Adsuar-Quesada J.M., Gil-Guillén V.F., Cortés-Castell E.
  *Mortality after out-of-hospital cardiac arrest in a Spanish region.*
  PLOS One, 12(4):e0175818, 2017.
  DOI: https://doi.org/10.1371/journal.pone.0175818


## Code information, installation, and requirements ##

To run the tests here, you need to have the following programs and packages installed in your computer:

* R (version > 4)
* R packages `pacman, dtw, randomForest, MLmetrics, crayon`

You can install them through these commands

    options(repos = list(CRAN="http://cran.rstudio.com/"))

    list.of.packages <- c("pacman")ì
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

    library("pacman")
    p_load("randomForest", "caret", "MLmetrics", "dtw", "crayon")

## Methodology ##

No preprocessing phase is needed on the raw data to run the scripts provided in this repository.
The target features of the datasets used are numeric or ordinal, and not categories.

## Usage instructions ##

To run the tests, you just need to execute the bash script this way:

    ./script_launching_repeated_holdout_and_cross_validation.sh

## License ##
This software package is released under the GNU General Public License v2.0.

## Article ##
Additional information about this project will be available in the following article:

> Davide Chicco, Giuseppe Jurman, "Repeated hold-out validation split can be more truthful than k-fold cross-validation split in supervised machine learning", 2026, in preparation.

## Contributions ##

This repository was developed by [Davide Chicco](https://www.DavideChicco.it). Questions and possible contributions should be addressed to davidechicco(AT)davidechicco.it

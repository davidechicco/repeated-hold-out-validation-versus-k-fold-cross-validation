# Repeated hold-out validation versus k-fold cross validation

Experiments to compare the results of repeated hold-out validation and the results of k-fold cross validation.

## Installation ##

To run the tests here, you need to have the following programs and packages installed in your computer:

* R (version > 4)
* R packages `pacman, dtw, randomForest, MLmetrics, crayon`

You can install them through these commands

    options(repos = list(CRAN="http://cran.rstudio.com/"))

    list.of.packages <- c("pacman")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)

    library("pacman")
    p_load("randomForest", "caret", "MLmetrics", "dtw", "crayon")

## Execution instructions ##

To run the tests, you just need to execute the bash script this way:

    ./script_launching_repeated_holdout_and_cross_validation.sh

## Article
Additional information about this project will be available in the following article:

> Davide Chicco, Giuseppe Jurman, "Repeated hold-out validation split can be more truthful than k-fold cross-validation split in supervised machine learning", 2026, in preparation.

## Contacts ##

This repository was developed by [Davide Chicco](https://www.DavideChicco.it). Questions should be addressed to davidechicco(AT)davidechicco.it

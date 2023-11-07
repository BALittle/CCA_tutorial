# A tutorial for running CCA for brain-behaviour associations in R

The script CCA_tutorial.R is intended to run through the process of using canonical correlation analysis (CCA) for investigating brain-behaviour associations.

It is written to accommodate cortical thickness data as exported from FreeSurfer in the default parcellation atlas (Desikan-Killiany).

## Data

The script requires two data sets saved in csv format:

-   'brain' data (e.g. cortical thickness for each brain region)

-   'behaviour' data (e.g., scores onseveral cognitive tests)

These datasets should be formatted such that each row represents a participant and each column represents a variable. All variables should be continuous, except the first column, which should contain subject IDs.

### Tutorial data

The tutorial uses a publicly available dataset of healthy adults (available at <https://openneuro.org/datasets/ds004215/versions/1.0.2>). This is used as 'dummy' data and has not been cleaned or pre-processed.

[**Please note**]{.underline}: If you use this code for your own data, the script assumes that the data has been cleaned (i.e., any skewed distributions and outliers have been dealt with, and missing data have been removed or imputed). There are some suggestions of how to check the data in the script. Any other pre-processing of the data should also be done before running the CCA_tutorial.R script, for example, accounting for confounding variables and z-scoring.

## Results

Results are saved in csv, txt and pdf files in the folder "Results".

------------------------------------------------------------------------

| Bethany Little [bethany.little\@ncl.ac.uk](mailto:bethany.little@ncl.ac.uk){.email}
| www.cnnp-lab.com
| 2023

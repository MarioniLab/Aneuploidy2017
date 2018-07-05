# Aneuploidy2017

## Contents

### Paper Supplementary files

This repo contains the latest version of the supplementary report for the paper *Mosaic autosomal aneuploidies are detectable from single-cell RNA-seq data*. This report is called `PaperSupplemental.html`

There are three R scripts that are used to generate this report:

* `PaperSupplemental.Rmd` is the code used to create figures and perform analyses.

* `aneu_functions.R` is a set of functions that were used in the analyses.

* `data_proc.R` is a script that converts raw data into a processed set for the analyses.

One further shell script, `get_data.sh` will download the raw data to this folder. You must have the folder containing this repo set as your working directory to download the data.

### *scploid* R package

Files for the *scploid* R package can be found in the `package` folder. However, it is easiest to install the package from R by installing the `devtools` package and using these commands:

`library(devtools)`

`devtools::install_github("MarioniLab/Aneuploidy2017", subdir = "package")`

A usage guide for the package is present, called `scploid_usage.html`. This file walks a user through the steps they would want to take in order to make an aneuploidy assessment on their own data.

The code used to generate it is also present, called `scploid_usage.Rmd`.

## Rendering the report yourself

Rendering of the report is computationally intensive, using ~60GB of RAM and 12 processing threads for around 3 hours of total CPU time. This is largely down to the simulations that are repeated to provide mean values. We therefore recommend running this report on a compute cluster if you wish to reproduce it in its entirity (some example submission scripts are in the `cluster` folder. However, it is likely sufficient to review the report if you seek clarity on any of our analyses (as code can be expanded in the .html file), or otherwise to extract and use particular functions from the `aneu_functions.R` file for your own analyses.

If you wish to compile the report, you will need to:

* Ensure the packages are installed that are loaded at the top of the `aneu_functions.R` file (best achieved using the `BiocInstaller` package).

* Run in R the command `rmarkdown::render(aneu_report.Rmd)` (recommended on a compute cluster; example scripts are present in the `cluster` folder)

The processed data (mostly *scploid* R objects) are available in this repo. Raw data may be downloaded using the `get_data.sh` script. By default, the `data_proc.R` script is disabled from running in the Rmd - you will need to un-comment the `source()` command in the first code chunk, and comment-out the data loading step if you would like to run this from the Rmd script.

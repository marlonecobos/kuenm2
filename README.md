kuenm2: Detailed Development of Ecological Niche Models
================
Marlon E. Cobos, Weverton Trindade, Luis F. Arias-Giraldo, Luis
Osorio-Olvera, and A. Townsend Peterson

- [Package description](#package-description)
- [Installing the package](#installing-the-package)
- [Workflow description](#workflow-description)
  - [Data preparation](#data-preparation)
  - [Model calibration](#model-calibration)
  - [Model explorations](#model-explorations)
  - [Model projections](#model-projections)
- [Variability and uncertainty](#variability-and-uncertainty)

<hr>

<br>

## Package description

The **kuenm2** R package implements multiple tools to help with the
development of detailed ecological niche models using distinct
algorithms. Pre-modeling analyses and explorations can be done to
prepare data. Model calibration (model selection) is done by training
and testing several candidate models. Handy options for producing final
models with transfers are included. Other tools to assess extrapolation
risks and variability in model transfers are also available.

<br>

## Installing the package

Note: Internet connection is required to install the package.

To install the latest release of **kuenm2** use the following line of
code:

<br>

The development version of **kuenm2** can be installed using the code
below.

``` r
# Installing and loading packages
if(!require(remotes)){
  install.packages("remotes")
}

# To install the package use
remotes::install_github("marlonecobos/kuenm2")

# To install the package and its vignettes use (if needed use: force = TRUE)  
#remotes::install_github("marlonecobos/kuenm2", build_vignettes = TRUE)  # in the process
```

<br>

*Having problems?*

If you have any problems during installation of the development version
from GitHub, restart R session, close other RStudio sessions you may
have open, and try again. If during the installation you are asked to
update packages, do so if you don’t need a specific version of one or
more of the packages to be installed. If any of the packages give an
error when updating, please install it alone using `install.packages()`,
then try installing **kuenm2** again.

<br>

To load the package use:

``` r
# Load kuenm2
library(kuenm2)
```

<br>

## Workflow description

### Data preparation

As shown in Fig. 1, to use **kuenm2** …

<br>

### Model calibration

After preparing data,

<br>

### Model explorations

### Model projections

After the selection of best parameter settings and the fitting of such
best performing models, projections to multiple scenarios can be done.
To facilitate projections to simple or complex combinations of
scenarios, use the **kuenm2** function `prepare_projection()`.

<br>

## Variability and uncertainty

Variability in results can be summarized and represented using multiple
options in **kuenm2**.

To represent uncertainty, the Mobility Oriented-Parity (MOP) metric is
implemented in **kuenm2**.

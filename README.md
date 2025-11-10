# Empirical Variational Bayes Factorization for Longitudinal Matrix Data (EVBLM)

R code to accompany "Empirical Variational Bayes Factorization for Longitudinal Matrix Data"
## Overview

This repository contains an implementation of an empirical variational Bayes (EVB) approach for factorizing high-dimensional longitudinal data. The method accommodates flexible covariance structures over time and automatically determines hyperparameters, including the underlying rank (number of factors).

## Installation

### Dependencies

The following R packages are required:

```r
install.packages(c("mvnfast", "MASS", "matrixcalc", "gplots", 
                   "Matrix", "flashier", "softImpute", "optimx"))
```

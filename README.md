# Empirical Variational Bayes Factorization for Longitudinal Matrix Data (EVBLM)

R code to accompany "Empirical Variational Bayes Factorization for Longitudinal Matrix Data"
## Overview

This repository contains an implementation of an empirical variational Bayes (EVB) approach for factorizing high-dimensional longitudinal data. The method accommodates flexible covariance structures over time and automatically determines hyperparameters, including the underlying rank (number of factors).

**Key features:**
- Flexible modeling of temporal correlation in longitudinal data
- Automatic rank selection
- Handles both aligned and irregularly-spaced timepoints across subjects
- Efficient missing data imputation (both entrywise and blockwise)
- Multiple covariance structures: Exchangeable (EC), Radial Basis Function (RBF), and Free covariance

## Installation

### Dependencies

The following R packages are required:

```r
install.packages(c("mvnfast", "MASS", "matrixcalc", "gplots", 
                   "Matrix", "flashier", "softImpute", "optimx"))
```

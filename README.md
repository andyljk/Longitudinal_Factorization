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

## Files

- **`EVBLM_functions.R`**: Core functions for the EVB longitudinal matrix factorization method with aligned timepoints
- **`EVBLM_irreg.R`**: Functions for handling irregularly-spaced timepoints across subjects
- **`EVBLM_example.R`**: Example usage demonstrating simulation, estimation, and imputation

## Quick Start

### Basic Usage (Aligned Timepoints)

```r
source("EVBLM_functions.R")

# Set parameters
set.seed(714)
n <- 50   # number of subjects
p <- 300  # number of variables
D <- 1:6  # timepoints
rank <- 3

# Generate simulated data with RBF covariance
# pars: c(noise_variance, eta, alpha)
dat_list <- sim_data(n, p, D, rk = rank, pars = c(1, 0.5, 3), fn = "RBF")

# Estimate model
est_rbf <- full_est(dat_list$X, D, fn = "RBF")

# Calculate relative squared error
sum((est_rbf$S - dat_list$S)^2) / sum(dat_list$S^2)
```

### Missing Data Imputation

```r
# Create 10% random missing data
X_miss <- dat_list$X
M <- sample(1:prod(dim(X_miss)), size = round(0.1 * prod(dim(X_miss))), replace = FALSE)
X_miss[M] <- NA

# Impute missing data
est_miss <- full_impute(X_miss, D, fn = "RBF", method = "greedy+backfit")

# Evaluate imputation accuracy
sum((est_miss$S[M] - dat_list$S[M])^2) / sum(dat_list$S[M]^2)
```

### Irregular Timepoints

For irregularly-spaced timepoints, use functions from `EVBLM_irreg.R`:

```r
source("EVBLM_irreg.R")

# X should be a list of length n, where each element is a p × d_i matrix
# D should be a list of length n, where each element is a vector of timepoints for that subject

# Example structure:
# X <- list(X_1, X_2, ..., X_n)  # each X_i is p × d_i
# D <- list(D_1, D_2, ..., D_n)  # each D_i is a vector of length d_i
```

## Output Structure

All estimation and imputation functions return a list with:

- **`S`**: Estimated low-rank signal (array)
- **`u`**: Loading factor estimates and posterior moments
- **`v`**: Score factor estimates and posterior moments
- **`elbo`**: Evidence Lower Bound values across iterations
- **`noise`**: Estimated noise variance(s)
- **`n_factors`**: Number of factors (rank) selected


## Citation

If you use this code, please cite:

```
[Citation will be added upon publication]
```

## Contact

For questions or issues, please open an issue on GitHub or contact the authors.

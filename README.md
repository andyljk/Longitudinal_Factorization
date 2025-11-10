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

#### Entrywise Missingness

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

#### Blockwise Missingness (Entire Timepoints Missing)

```r
# Create blockwise missingness (10% of columns missing)
X_block <- dat_list$X
X_unfold <- do.call(cbind, lapply(D, function(d) X_block[,,d]))
M <- sample(1:(p * length(D)), size = round(0.1 * p * length(D)), replace = FALSE)
X_unfold[,M] <- NA
X_block <- array(X_unfold, dim = c(n, p, length(D)))
M <- which(is.na(X_block))

# Impute blockwise missing data
est_block <- full_impute(X_block, D, fn = "RBF", method = "greedy+backfit")

# Evaluate imputation accuracy
sum((est_block$S[M] - dat_list$S[M])^2) / sum(dat_list$S[M]^2)
```

## Main Functions

### Data Simulation

- **`sim_data(n, p, D, rk, pars, fn)`**: Generate simulated longitudinal data
  - `n`: Number of subjects
  - `p`: Number of variables
  - `D`: Vector of timepoints
  - `rk`: True rank
  - `pars`: Parameters for covariance structure (varies by `fn`)
  - `fn`: Covariance function ("EC", "RBF", "Free", or "Block")

### Model Estimation

- **`full_est(X, D, fn, thres, verbose)`**: Complete estimation using greedy initialization + backfitting
  - `X`: Data array (n × p × length(D))
  - `D`: Vector of timepoints
  - `fn`: Covariance structure ("EC", "RBF", or "Free")
  - `thres`: Convergence threshold (default: 1e-8)
  - `verbose`: Print progress (default: FALSE)

- **`backfit(X, D, R, fn, u_mt, l_mt, thres, max_iter)`**: Coordinate ascent algorithm for rank R model
- **`greedy_fac(X, D, fn, thres)`**: Greedy sequential rank-1 fitting

### Missing Data Imputation

- **`full_impute(X, D, R, fn, method, thres, max_iter, verbose)`**: Impute missing data
  - `X`: Data array with missing values (NA)
  - `D`: Vector of timepoints
  - `R`: Rank (if NULL, determined automatically with method="greedy+backfit")
  - `fn`: Covariance structure
  - `method`: "greedy+backfit" (recommended) or "backfit"
  - `thres`: Convergence threshold (default: 1e-5)
  - `max_iter`: Maximum iterations (default: 100)

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

## Covariance Structures

The method supports several covariance functions for modeling temporal correlation:

1. **Exchangeable Covariance (EC)**: `K(s,t) = ρ·σ²_v` for s ≠ t
   - Parameters: `(sigma_v, rho)`
   - Good for: Consistent correlation across all timepoint pairs

2. **Radial Basis Function (RBF)**: `K(s,t) = σ²_v · exp(-||s-t||²/(2α²))`
   - Parameters: `(sigma_v, alpha)`
   - Good for: Smooth decay of correlation with time distance

3. **Free Covariance**: Unstructured covariance matrix
   - Parameters: Full covariance matrix
   - Good for: Maximum flexibility when structure is unclear

## Output Structure

All estimation and imputation functions return a list with:

- **`S`**: Estimated low-rank signal (array)
- **`u`**: Loading factor estimates and posterior moments
- **`l` or `v`**: Score factor estimates and posterior moments
- **`elbo`**: Evidence Lower Bound values across iterations
- **`noise`**: Estimated noise variance(s)
- **`n_factors`**: Number of factors (rank) selected

## Simulation Studies

The code reproduces the simulation experiments from the paper:

- **Complete data signal recovery** (Section 3.1)
- **Entrywise missingness imputation** (Section 3.2)
- **Blockwise missingness imputation** (Section 3.2)

See `EVBLM_example.R` for implementation details.

## Citation

If you use this code, please cite:

```
[Citation will be added upon publication]
```

## Contact

For questions or issues, please open an issue on GitHub or contact the authors.

## License

[License information to be added]

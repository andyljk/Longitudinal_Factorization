setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("EVBLM_functions.R")


set.seed(714)
n = 50
p = 300
D = 1:6
rank = 3
# returns two arrays true underlying data (S), and noisy data (X) in a list
dat_list = sim_data(n,p,D,rk=rank,pars = c(1,0.5,3), fn="RBF") # for RBF, the pars represent (sigma^2, eta, alpha)

est_rbf = full_est(dat_list$X, D, fn="RBF")

# show the relative mean squared error
sum((est_rbf$S - dat_list$S)^2)/sum(dat_list$S^2)


# imputation
# entrywise missingness
X_miss = dat_list$X
M = sample(1:prod(dim(X_miss)), size= round(0.1*prod(dim(X_miss))), replace=F)
X_miss[M] = NA

est_miss = full_impute(X_miss, D, fn="RBF", method="greedy+backfit")
sum((est_miss$S[M]-dat_list$S[M])^2)/sum(dat_list$S[M]^2)

# blockwise missingness
X. = dat_list$X
X_unfold = do.call(cbind, lapply(D, function(d) X.[,,d])) # cbind the frontal slices
M = sample(1:(p*length(D)), size= round(0.1*p*length(D)), replace=F)
X_unfold[,M] = NA
X. = array(X_unfold,dim=c(n,p,length(D)))
M = which(is.na(X.))
est_block = full_impute(X., D, fn="RBF", method="greedy+backfit")
sum((est_block$S[M]-dat_list$S[M])^2)/sum(dat_list$S[M]^2)

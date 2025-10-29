library(mvnfast)
library(MASS)
library(matrixcalc)
library(gplots)
library(Matrix)
library(flashier)
#source("/Users/andyl/Desktop/Bayesian_BIDIFAC/EBMF.R")

avg_loglik = function(par,X,S,D){ # avg log likelihood for GP with the RBF kernel
  p = nrow(X)
  # D = ncol(X)
  eta_hat = exp(par[1])
  alp_hat = exp(par[2])
  pairwise_diffs = outer(D, D, function(i, j) (i-j)^2)
  K <- eta_hat^2 * exp(-pairwise_diffs / (2*alp_hat^2))
  #print(c(eta_hat,alp_hat))
  sum(dmvn(X, mu=rep(0,length(D)), sigma=K+S, log=T))/p/1000000
}
grad = function(par,X,S,D){
  p = nrow(X)
  #D = ncol(X)
  t1_hat = par[1] # t1 = log(eta)
  t2_hat = par[2] # t2 = log(alpha)
  pairwise_diffs <- outer(D, D, function(i, j) (i-j)^2)
  # Apply the RBF covariance function
  K_y <- exp(2*t1_hat) * exp(-pairwise_diffs / (2*exp(2*t2_hat))) + S
  K_y_inv = solve(K_y) # pre compute the inverse of K_y
  # element-wise deriv of K_y wrt t1 & t2
  dK_dt1 = 2*exp(2*t1_hat) * exp(-pairwise_diffs / (2*exp(2*t2_hat))) 
  dK_dt2 = exp(2*t1_hat) * exp(-pairwise_diffs / (2*exp(2*t2_hat)))*pairwise_diffs/exp(2*t2_hat)
  
  gam = K_y_inv%*%t(X) # intermediate term to speed up computation
  
  #partial derivs of log_lik wrt t1 & t2
  dl_dt1 = 1/2*sum(diag(gam%*%t(gam)%*%dK_dt1))/p - 1/2*sum(diag(K_y_inv%*%dK_dt1))
  dl_dt2 = 1/2*sum(diag(gam%*%t(gam)%*%dK_dt2))/p - 1/2*sum(diag(K_y_inv%*%dK_dt2))
  
  c(dl_dt1, dl_dt2)/1000000
}

GridSearchRBF = function(X,S,D,fine=0.1,range_max=10){
  p = nrow(X)
  #D = ncol(X)
  pairwise_diffs <- outer(D, D, function(i, j) (i-j)^2)
  range. = seq(fine,range_max,fine)
  N = length(range.)
  log_lik = matrix(,nrow=N,ncol=N)
  for (a in 1:N){
    for (b in 1:N){ # iterate through all value pairs of eta and alpha
      log_lik[a,b] = avg_loglik(log(c(range.[a],range.[b])),X,S,D)
    }
  }
  ind = which(log_lik==max(log_lik),arr.ind=T)
  return(fine*c(ind[1],ind[2]))
}


mvEBNM = function(X,s,D,fn = c("EC","Free","RBF")){
  results = list() # a list to save everything: prior, par, loglik, negKL (for obj func calc)
  p = dim(X)[1] # number of variables in model (samples in mvEBNM)
  d = length(D)#D = dim(X)[2] # number of longitudinal timepoints in model (dimension in mvEBNM)
  S = diag(s^2,nrow=length(s)) # diag matx of known variances in mvEBNM
  results$prior = fn # prior structure
  
  # Empirical estimation of prior via max marginal likelihood
  if (fn == "EC"){ # equicorrelation model (compound symmetry)
    sample_cov = t(X)%*%X/p
    eta_hat2 = max(mean(diag(sample_cov)-s^2),1e-50)
    rho_hat = min(0.95,max((sum(sample_cov) - sum(diag(sample_cov)))/(d*(d-1))/(eta_hat2),1e-50))
    if (eta_hat2 < 1e-49){rho_hat = 0}
    V = (matrix(rho_hat,nrow=d,ncol=d) + diag(1-rho_hat, nrow=d))%*%diag(eta_hat2,nrow=d)
    results$par = c(eta_hat2^0.5, rho_hat)
    names(results$par) = c("eta","rho")
    results$loglik = sum(dmvn(X, mu=rep(0,d), sigma = V+S, log=T))
  }else if (fn=="Free"){ # unparametrized covariance matrix
    sample_cov = t(X)%*%X/p
    V = sample_cov - diag(s^2,nrow=length(s))
    if (any(diag(V)<=0)){
      V=diag(1e-3,nrow=d)
    }
    results$par = V
    names(results$par) = "Cov_Mat"
    results$loglik = sum(dmvn(X, mu=rep(0,d), sigma = V+S, log=T))
  }else if (fn=="RBF"){ # Gaussian Process with RBF kernel
    pairwise_diffs = outer(D, D, function(i, j) (i-j)^2)
    par_init = log(GridSearchRBF(X,S,D,fine=0.2,range_max = 10))
    if (exp(par_init[1]) == 0.1){par_init[1] = log(GridSearchRBF(X,S,fine=0.01,range_max = 0.2)[1])}
    res = optim(par=par_init,fn = avg_loglik, grad, method="BFGS",control=list(fnscale=-1),X=X, S=S,D=D)
    eta_hat = max(exp(res$par[1]),1e-50)
    alp_hat = exp(res$par[2])
    V <- eta_hat^2 * exp(-pairwise_diffs / (2*alp_hat^2))+1e-5
    results$loglik = p*res$value*1000000
    results$par = c(eta_hat, alp_hat)
    names(results$par) = c("eta","alpha")
  }
  
  # Calculate the 1st and 2nd posterior moments
  S_inv = solve(diag(s^2,nrow=length(s)))
  V_inv = svd.inverse(V)#V_inv = solve(V)
  results$pos = array(dim = c(p,2,d))
  results$pos[,1,] = t(svd.inverse(S_inv+V_inv)%*%(S_inv%*%t(X)))
  results$pos[,2,] = t(apply(results$pos[,1,], 1, function(row) diag(row %*% t(row) + svd.inverse(S_inv + V_inv))))
  
  # Calculate the negative KL divergence that will be used in objective function calculation
  term = 0
  for (j in 1:p){
    term = term + d/2*log(2*pi) + 1/2*sum(diag(log(S))) + 1/2*(t(X[j,])%*%S_inv%*%X[j,] - 2*t(X[j,])%*%S_inv%*%results$pos[j,1,] + sum(1/s^2*results$pos[j,2,]))
  }
  results$neg_KL = results$loglik + term
  return(results)
}

# univariate EBNM function
ebnm_gaussian = function(x,s){
  n = length(x)
  results = list() # a list to save everything: prior, par, loglik, negKL (for obj func calc)
  results$prior = "Normal"
  mu_hat = 0
  eta_hat2 = max((mean(x^2)-s[1]^2),1e-50)
  results$par = c(0, eta_hat2^0.5)
  names(results$par) = c("mu","eta")
  results$loglik = sum(dnorm(x,mu_hat,sqrt(s^2+eta_hat2),log=TRUE))
  results$pos = array(dim = c(n,2))
  results$pos[,1] = (eta_hat2*x+mu_hat*s[1]^2)/(eta_hat2+s[1]^2)
  results$pos[,2] = ((eta_hat2*x+mu_hat*s[1]^2)/(eta_hat2+s[1]^2))^2+eta_hat2*(s[1]^2)/(eta_hat2+s[1]^2)
  results$neg_KL = results$loglik + 0.5*sum(log(2*pi*s^2) + 1/s^2*(x^2+results$pos[,2]-2*x*results$pos[,1]))
  return(results)
}



fnorm = function(X){
  return(sum(X^2)/prod(dim(X)))
}

get_Rt = function(X_t, u, u_2, l, l_2){
  u = as.matrix(u)
  u_2 = as.matrix(u_2)
  l = as.matrix(l)
  l_2 = as.matrix(l_2)
  n=dim(X_t)[1]
  p=dim(X_t)[2]
  R = matrix(NA,nrow=n,ncol=p)
  for (i in 1:n){
    for (j in 1:p){
      R[i,j] = (X_t[i,j]-sum(u[i,]*l[j,]))^2 - sum(u[i,]^2*l[j,]^2) + sum(u_2[i,]*l_2[j,])
    }
  }
  return(R)
}

x_u = function(X, l_mt, tau){ # X_1-3 are longitudinal data matrices, l's should be size p
  vec = c()
  n = dim(X)[1]
  for (i in 1:n){
    vec[i] = sum((l_mt[,1,]*X[i,,])%*%diag(tau,nrow=length(tau)))/sum(l_mt[,2,]%*%diag(tau,nrow=length(tau)))
  }
  return(vec)
}

x_l = function(X, u_mt, tau){ # X is a n*p*d array (mode 3 tensor), dim(u_mt)=p,2 ,tau are precisions length d
  p = dim(X)[2] # this should equal to p
  u = u_mt[,1]
  u_2 = u_mt[,2]
  S = diag(1/tau, nrow=length(tau))
  mat = matrix(,nrow=p,ncol=length(tau))
  for (j in 1:p){
    mat[j,] = (1/sum(u_2)) %*% colSums(u*(X[,j,]))
  }
  return(mat)
}

s_u = function(l_mt,tau){
  return(max(sum(l_mt[,2,]%*%diag(tau,nrow=length(tau))),1e-10)^(-1/2))
}

s_l = function(u_mt, tau){
  S = diag(1/tau, nrow=length(tau))
  return(1/(sum(u_mt[,2])*tau)^(0.5))
}

get_F = function(X, D, u_mt, l_mt, tau){ # u_mt and l_mt have $pos, $loglik, $par
  n = dim(X)[1]
  p = dim(X)[2]
  d = dim(X)[3]#D = dim(X)[3]
  loglik = 0
  for (t in 1:d){
    R_t = get_Rt(X[,,t], u_mt$pos[,,1], u_mt$pos[,,2], l_mt$pos[,,1,t], l_mt$pos[,,2,t])
    loglik = -n*p/2*log(2*pi/tau[t]) - tau[t]/2*sum(R_t) + loglik
  }
  Free_energy = loglik + sum(u_mt$neg_KL) + sum(l_mt$neg_KL)
  return(Free_energy)
}

init_factors = function(X,D,R){
  n = dim(X)[1]
  p = dim(X)[2]
  #D = dim(X)[3]
  X_unfold = do.call(cbind, lapply(1:length(D), function(d) X[,,d])) # cbind the frontal slices
  u_mt = l_mt = list()      
  svd_init = svd(X_unfold)
  s_vals = sqrt(svd_init$d[1:R])
  u = svd_init$u[,1:R]
  u_mt$pos = array(c(u%*%diag(s_vals,nrow=R),u^2%*%diag(s_vals^2,nrow=R)), dim=c(n,R,2)) # posteriors as n*rk*2 array for u_mt
  l_mt$pos = array(unlist(lapply(1:length(D), function(d) { # posteriors as n*rk*2*length(D) array for l_mt
    l_d <- svd_init$v[((d-1)*p + 1):(d*p), 1:R]
    c(l_d%*%diag(s_vals,nrow=R), l_d^2%*%diag(s_vals^2,nrow=R))
  })), dim = c(p, R, 2, length(D)))
  u_mt$par = l_mt$par = list() # $par as length rk list (since each column may have > 1 prior parameters)
  return(list(u_mt,l_mt))
}

update_u_factor = function(u_mt, hold_u, k){
  u_mt$par[[k]] = hold_u$par
  u_mt$pos[,k,] = hold_u$pos
  u_mt$loglik[k] = hold_u$loglik
  u_mt$neg_KL[k] = hold_u$neg_KL
  return(u_mt)
}

update_l_factor = function(l_mt, hold_l, k){
  l_mt$par[[k]] = hold_l$par
  l_mt$pos[,k,,] = hold_l$pos
  l_mt$loglik[k] = hold_l$loglik
  l_mt$neg_KL[k] = hold_l$neg_KL
  return(l_mt)
}

single_factor = function(X, D, fn=c("EC","Free","RBF"), thres=1e-8, max_iter = 50){
  n = dim(X)[1]
  p = dim(X)[2]
  d = length(D) # D = dim(X)[3]
  R = 1
  init = init_factors(X,D,R)
  u_mt = init[[1]]
  l_mt = init[[2]]
  S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
  S_old = array(Inf, dim=c(n,p,d))
  tau = rep(1,d)
  F_vec=c()
  count = 0
  while(fnorm(S_new - S_old)>thres & count <= max_iter){
    S_old = S_new
    for (m in 1:d){
      R_m = get_Rt(X[,,m],u_mt$pos[,,1],u_mt$pos[,,2],l_mt$pos[,,1,m],l_mt$pos[,,2,m])
      tau[m] = n*p/sum(R_m)
    }
    for (k in 1:R){
      R_k = array(sapply(1:d, function(m) X[,,m] - u_mt$pos[,-k,1] %*% t(l_mt$pos[,-k,1,m])), dim = c(n, p, d))
      hold_u = ebnm_gaussian(x = x_u(R_k, l_mt$pos[,k,,], tau), s = s_u(l_mt$pos[,k,,],tau))
      u_mt = update_u_factor(u_mt,hold_u,k)
      hold_l = mvEBNM(X = x_l(R_k, u_mt$pos[,k,], tau), s = s_l(u_mt$pos[,k,], tau), D=D, fn=fn)
      l_mt = update_l_factor(l_mt,hold_l,k)
    }
    F_vec = c(F_vec,get_F(X,D,u_mt,l_mt,tau))
    S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
    count=count+1
    #print(fnorm(S_new - S_old))
  }
  results = list(u = u_mt, l = l_mt, elbo = F_vec, S = S_new, noise = tau, n_iter = count)
  return(results)
}

backfit = function(X, D, R, fn = c("EC","Free","RBF"), 
                   u_mt = NULL, l_mt = NULL, thres=1e-8, 
                   verbose = F, max_iter = 50, null_check=T){
  n = dim(X)[1];p = dim(X)[2];d = length(D) # D = dim(X)[3]
  if (is.null(u_mt) | is.null(l_mt)){
    init = init_factors(X,D,R=R)
    u_mt = init[[1]];l_mt = init[[2]]
  }
  S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
  S_old = array(Inf, dim=c(n,p,d))
  
  tau = rep(1,d)
  
  F_vec=c()
  count = 0
  while(fnorm(S_new - S_old)>thres & count <= max_iter){
    S_old = S_new
    for (m in 1:d){
      R_m = get_Rt(X[,,m],u_mt$pos[,,1],u_mt$pos[,,2],l_mt$pos[,,1,m],l_mt$pos[,,2,m])
      tau[m] = n*p/sum(R_m)
    }
    for (k in 1:R){
      R_k = array(sapply(1:d, function(m) X[,,m] - u_mt$pos[,-k,1] %*% t(l_mt$pos[,-k,1,m])), dim = c(n, p, d))
      hold_u = ebnm_gaussian(x = x_u(R_k, l_mt$pos[,k,,], tau), s = s_u(l_mt$pos[,k,,],tau))
      u_mt = update_u_factor(u_mt,hold_u,k)
      hold_l = mvEBNM(X = x_l(R_k, u_mt$pos[,k,], tau), s = s_l(u_mt$pos[,k,], tau), D=D, fn=fn)
      l_mt = update_l_factor(l_mt,hold_l,k)
    }
    F_vec = c(F_vec,get_F(X,D,u_mt,l_mt,tau))
    S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
    count=count+1
    if (verbose){print(fnorm(S_new - S_old))}
  }
  results = list(u = u_mt, l = l_mt, elbo = F_vec, S = S_new, noise = tau, n_factors = R)
  if (null_check){results = nullcheck(results,D)}
  return(results)
}

nullcheck = function(results,D){
  n = dim(results$S)[1]
  p = dim(results$S)[2]
  d = length(D) # D = dim(results$S)[3]
  r = dim(results$u$pos)[2]
  if (any(colSums(abs(matrix(results$u$pos[,,1],ncol=r)))<1e-8) & any(colSums(abs(matrix(results$l$pos[,,1,],ncol=r*d)))<1e-8)){
    ind = which(colSums(abs(matrix(results$u$pos[,,1],ncol=r)))<1e-8)
    results$u$pos = results$u$pos[,-ind,]
    results$l$pos = results$l$pos[,-ind,,]
    results$u$par = results$u$par[-ind]
    results$l$par = results$l$par[-ind]
    results$u$loglik = results$u$loglik[-ind]
    results$l$loglik = results$l$loglik[-ind]
    results$u$neg_KL = results$u$neg_KL[-ind]
    results$l$neg_KL = results$l$neg_KL[-ind]
    results$n_factors = results$n_factors - length(ind)
  }
  return(results)
}

greedy_fac = function(X, D, fn=c("EC","Free","RBF"), thres=1e-8){
  n = dim(X)[1];p = dim(X)[2];d = length(D) # D = dim(X)[3]
  u_mt = l_mt = list()
  u_mt$pos = array(0,dim=c(n,min(n,p),2))
  l_mt$pos = array(0,dim=c(p,min(n,p),2,d))
  u_mt$par = l_mt$par = list()
  S_new = array(0,dim=c(n,p,d))
  k=0
  done=F
  F_vec = c()
  while (!done){
    R_k = X - S_new
    new = single_factor(R_k, D=D, fn=fn, thres=thres)
    if (abs(sum(new$u$pos[,,1])) < 1e-10 & abs(sum(new$l$pos[,,1,])) < 1e-10){break}
    k = k+1
    u_mt$pos[,k,] = new$u$pos
    u_mt$par[[k]] = new$u$par[[1]]
    u_mt$loglik[k] = new$u$loglik
    u_mt$neg_KL[k] = new$u$neg_KL
    l_mt$pos[,k,,] = new$l$pos
    l_mt$par[[k]] = new$l$par[[1]]
    l_mt$loglik[k] = new$l$loglik
    l_mt$neg_KL[k] = new$l$neg_KL
    tau = new$noise
    if (k>1){
      if (get_F(X, D, u_mt, l_mt, tau) <= F_vec[k-1]){
        k = k-1
        break
      }
    }
    F_vec[k] = get_F(X, D, u_mt, l_mt, tau)
    S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
  }
  u_mt$pos = array(u_mt$pos[,1:k,], dim=c(n,k,2))
  l_mt$pos = array(l_mt$pos[,1:k,,], dim=c(p,k,2,d))
  return(list(u = u_mt, l = l_mt, elbo = F_vec, S = S_new, noise = tau, n_factors = k))
}

full_est = function(X, D, fn=c("EC","Free","RBF"), thres=1e-8, verbose=F){
  print("Initializing via Greedy estimation")
  g_init = greedy_fac(X,D,fn=fn,thres=thres)
  paste0("Estimated factors:", g_init$n_factors)
  print("Refining via full loop algorithm")
  results = backfit(X,D,R=g_init$n_factors,fn=fn,u_mt=g_init$u,l_mt=g_init$l,thres=thres,verbose=verbose)
  print("done")
  return(results)
}

sim_data = function(n,p,D,rk,pars,fn=c("EC","Free","RBF")){
  tau. = rep(1,rk)
  d = length(D)
  if (fn=="RBF"){
    if (length(pars)!=3){return("# pars not correct (RBF)")}
    eta = pars[2]
    alp = pars[3]
    pairwise_diffs <- outer(D, D, function(i, j) (i-j)^2)
    K. <- eta^2 * exp(-pairwise_diffs / (2*alp^2)) # True GP prior covariance
  }else if (fn=="EC"){
    if (length(pars)!=3){return("# pars not correct (EC)")}
    eta = pars[2]
    rho = pars[3]
    K. = (array(rho,dim=c(d,d)) + diag(1-rho,nrow=d))*eta^2
  }else if (fn=="Free"){
    if (length(pars)!=2){return("# pars not correct (Free)")}
    A = array(rnorm(d^2,0,pars[2]),dim=c(d,d))
    K. = t(A)%*%A/d
  }else if (fn=="Block"){
    if (length(pars)!=4){return("# pars not correct (Block)")}
    eta = pars[2]
    rho = pars[3]
    cut = pars[4]
    K_1 = (array(rho,dim=c(cut,cut)) + diag(1-rho,nrow=cut))*eta^2
    K_2 = (array(rho,dim=c(d-cut,d-cut)) + diag(1-rho,nrow=d-cut))*eta^2
    K. = as.matrix(bdiag(K_1,K_2))
  }
  
  sig = rep(pars[1],d)
  U. = matrix(rnorm(n*rk, 0, 1),nrow = n) %*% diag(tau.,nrow=length(tau.))
  L. = array(,dim = c(p,rk,d))
  for (k in 1:rk){
    L.[,k,] = rmvn(n=p, mu=rep(0,d), sigma = K.)
  }
  S. = X. = array(,dim = c(n,p,d))
  for (m in 1:d){
    S.[,,m] = U.%*%t(L.[,,m])
    X.[,,m] = S.[,,m] + matrix(rnorm(p*n,0,sig[m]^0.5), nrow = n)
  }
  return(list(S = S., X = X.))
}

show.image = function(Image,sd=NULL,ylab=''){
  if (is.null(sd)){s=sd(Image)}
  else{s=sd}
  
  lower = mean(Image)-3*s
  upper = mean(Image)+3*s
  Image[Image<lower] = lower
  Image[Image>upper] = upper
  image(x=1:dim(Image)[2], y=1:dim(Image)[1], z=t(Image), zlim = c(lower,upper),axes=FALSE,col=bluered(100),xlab="",ylab=ylab)
} 

full_impute = function(X,D,R=NULL,fn=c("EC","RBF","Free"), 
                       method=c("greedy+backfit","backfit"),
                       thres=1e-5,max_iter=100,verbose=F){
  n = dim(X)[1]; p=dim(X)[2]; d = length(D) # D=dim(X)[3]
  M = which(is.na(X))
  X[M] = 0
  count = 0
  if (method=="greedy+backfit"){
    g_init = greedy_fac(X,D,fn=fn,thres=thres)
    R=g_init$n_factors; u_mt=g_init$u; l_mt=g_init$l
  }else {
    if (is.null(R)){return("n factors required for backfit only algorithm")}
    init = init_factors(X,D,R=R)
    u_mt = init[[1]];l_mt = init[[2]]
  }
  S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
  S_old = array(Inf, dim=c(n,p,d))
  tau = rep(1,d)
  
  F_vec=c(); count = 0
  while(fnorm(S_new[M] - S_old[M])>thres & count <= max_iter){
    S_old = S_new
    for (m in 1:d){
      R_m = get_Rt(X[,,m],u_mt$pos[,,1],u_mt$pos[,,2],l_mt$pos[,,1,m],l_mt$pos[,,2,m])
      tau[m] = n*p/sum(R_m)
    }
    for (k in 1:R){
      R_k = array(sapply(1:d, function(m) X[,,m] - u_mt$pos[,-k,1] %*% t(l_mt$pos[,-k,1,m])), dim = c(n, p, d))
      hold_u = ebnm_gaussian(x = x_u(R_k, l_mt$pos[,k,,], tau), s = s_u(l_mt$pos[,k,,],tau))
      u_mt = update_u_factor(u_mt,hold_u,k)
      hold_l = mvEBNM(X = x_l(R_k, u_mt$pos[,k,], tau), s = s_l(u_mt$pos[,k,], tau), D=D, fn=fn)
      l_mt = update_l_factor(l_mt,hold_l,k)
    }
    F_vec = c(F_vec,get_F(X,D,u_mt,l_mt,tau))
    S_new = array(sapply(1:d, function(m) u_mt$pos[,,1]%*%t(l_mt$pos[,,1,m])), dim = c(n,p,d))
    X[M] = S_new[M]
    count=count+1
    if (verbose){print(fnorm(S_new[M] - S_old[M]))}
  }
  results = list(u = u_mt, l = l_mt, elbo = F_vec, S = S_new, noise = tau, n_factors = R)
  results = nullcheck(results,D)
  return(results)
}


# full_impute = function(X,D,fn=c("EC","Free","RBF"),imp_thres=1e-5,thres=1e-5,max_iter=30){
#   n = dim(X)[1]; p=dim(X)[2]; d = length(D) # D=dim(X)[3]
#   M = which(is.na(X))
#   X[M] = 0
#   S_new = X
#   S_old = array(Inf, dim=c(n,p,d))
#   count = 0
#   while (fnorm(S_new[M] - S_old[M]) > imp_thres & count <= max_iter){
#     S_old = S_new
#     results = full_est(X, D, fn=fn, thres=thres) # backfit(X, D, R, fn=fn, thres=thres) #
#     S_new = results$S
#     X[M] = S_new[M]
#     print("impute iter:")
#     print(fnorm(S_new[M] - S_old[M]))
#     #print(sum((S_new[M] - S.[M])^2)/sum(S.[M]^2))
#     count = count + 1
#   }
#   return(results)
# }
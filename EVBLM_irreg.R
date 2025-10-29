library(mvnfast)
library(matrixcalc)
library(gplots)
library(softImpute)
library(optimx)

avg_loglik_rbf = function(par,X,s,D){ # avg log likelihood for GP with the RBF kernel
  n = length(D)
  eta_hat = exp(par[1])
  alp_hat = exp(par[2])
  loglik = mapply(
    function(Xi, Di) {
      S <- diag(s^2, nrow = length(Di))  # Create the diagonal matrix
      pairwise_diffs <- outer(Di, Di, function(a, b) (a - b)^2)  # Pairwise differences
      K <- eta_hat^2 * exp(-pairwise_diffs / (2 * alp_hat^2))  # Compute K
      dmvn(Xi, mu = rep(0, length(Di)), sigma = K+S, log = TRUE)  # Compute log-likelihood
    },
    X, D
    #SIMPLIFY = TRUE  # Ensure output is a vector
  )
  #print(c(eta_hat,alp_hat))
  value = -sum(loglik)/n
  # attr(value, 'gradient') <- grad(par,X,s,D)
  return(value)
}

avg_loglik_ec = function(par,X,s,D){ # avg log likelihood for GP with the EC kernel
  n = length(D)
  eta2_hat = par[1]
  rho_hat = par[2]
  loglik = mapply(
    function(Xi, Di) {
      S <- diag(s^2, nrow = length(Di))  # Create the diagonal matrix
      K = (array(rho_hat,dim=c(length(Di),length(Di))) + diag(1-rho_hat,nrow=length(Di)))*eta2_hat # compute K
      dmvn(Xi, mu = rep(0, length(Di)), sigma = K+S, log = TRUE)  # Compute log-likelihood
    },
    X, D
    #SIMPLIFY = TRUE  # Ensure output is a vector
  )
  #print(c(eta_hat,alp_hat))
  value = -sum(loglik)/n
  return(value)
}

grad = function(par,X,s,D){
  n = length(X)
  t1_hat = par[1] # t1 = log(eta)
  t2_hat = par[2] # t2 = log(alpha)
  dl_dt1 = dl_dt2 = 0
  # for (i in 1:n){
  #   pairwise_diffs <- outer(D[[i]], D[[i]], function(a, b) (a-b)^2)
  #   Si = diag(s^2, nrow=length(D[[i]]))
  #   K_i <- exp(2*t1_hat) * exp(-pairwise_diffs / (2*exp(2*t2_hat))) + Si
  #   K_i_inv = solve(K_i) # pre compute the inverse of K_y
  #   # element-wise deriv of K_y wrt t1 & t2
  #   dK_dt1 = 2*exp(2*t1_hat) * exp(-pairwise_diffs / (2*exp(2*t2_hat))) 
  #   dK_dt2 = exp(2*t1_hat) * exp(-pairwise_diffs / (2*exp(2*t2_hat)))*pairwise_diffs/exp(2*t2_hat)
  #   gam = K_i_inv%*%t(X[[i]])
  #   #partial derivs of log_lik wrt t1 & t2
  #   dl_dt1 = 0.5/n*sum(diag( (gam%*%t(gam)-K_i_inv)%*%dK_dt1 )) + dl_dt1
  #   dl_dt2 = 0.5/n*sum(diag( (gam%*%t(gam)-K_i_inv)%*%dK_dt2 )) + dl_dt2
  # }
  partials <- lapply(1:n, function(i) {
    D_i <- D[[i]]
    pairwise_diffs <- outer(D_i, D_i, function(a, b) (a - b)^2)
    Si <- diag(s^2, nrow = length(D_i))
    K_i <- exp(2 * t1_hat) * exp(-pairwise_diffs / (2 * exp(2 * t2_hat))) + Si
    K_i_inv <- solve(K_i)
    
    dK_dt1 <- 2 * exp(2 * t1_hat) * exp(-pairwise_diffs / (2 * exp(2 * t2_hat)))
    dK_dt2 <- exp(2 * t1_hat) * exp(-pairwise_diffs / (2 * exp(2 * t2_hat))) * pairwise_diffs / exp(2 * t2_hat)
    
    gam <- K_i_inv %*% t(X[[i]])
    common_term <- gam %*% t(gam) - K_i_inv
    
    dl_dt1_i <- 0.5 / n * sum(diag(common_term %*% dK_dt1))
    dl_dt2_i <- 0.5 / n * sum(diag(common_term %*% dK_dt2))
    
    c(dl_dt1_i, dl_dt2_i)
  })
  
  # Sum over all contributions from each i
  partials_mat <- do.call(rbind, partials)
  dl_dt1 <- sum(partials_mat[, 1])
  dl_dt2 <- sum(partials_mat[, 2])
  
  return(-c(dl_dt1, dl_dt2))
}

GridSearchRBF = function(X,s,D,fine=0.1,range_max=10){
  range. = seq(fine,range_max,fine)
  N = length(range.)
  log_lik = matrix(,nrow=N,ncol=N)
  for (a in 1:N){
    for (b in 1:N){
      log_lik[a,b] = avg_loglik_rbf(par=log(c(range.[a],range.[b])), X=X, s=s, D=D)
    }
  }
  ind = which(log_lik==min(log_lik),arr.ind=T)
  return(fine*c(ind[1],ind[2]))
}


get_R = function(X, D, u_mt, v_mt){
  n = length(X)
  p = dim(X[[1]])[1]
  R_2 = X
  for (i in 1:n){
    d_i = length(D[[i]])
    for (j in 1:p){
      for (d in 1:d_i){
        R_2[[i]][j,d] = (X[[i]][j,d] - sum(u_mt$pos[j,,1]*v_mt$pos[[i]][,d,1]))^2 - 
          sum(u_mt$pos[j,,1]^2*v_mt$pos[[i]][,d,1]^2) + sum(u_mt$pos[j,,2]*v_mt$pos[[i]][,d,2])
      }
    }
  }
  return(R_2)
}


x_u = function(X, v_mt, tau, k){ # should produce a vector of length n
  p = dim(X[[1]])[1]
  vec = c(); denom = sum(sapply(v_mt$pos, function(x) sum(x[k,,2])))*tau
  for (j in 1:p){ # X is a length n list of p*d_i matrices
    numer = sum(sapply(1:length(v_mt$pos), function(i) sum(v_mt$pos[[i]][k,,1] * X[[i]][j,])))*tau
    vec[j] = numer/denom
  }
  return(vec)
}

x_v = function(X, u_mt, tau){ # should produce a list of row matrices of corresponding vectors of irrEBNM
  n = length(X)
  x_list = vector(mode="list",length=n)
  for (i in 1:n){
    x_list[[i]] = matrix(colSums(u_mt[,1]*X[[i]])/sum(u_mt[,2]),nrow=1)
  }
  return(x_list)
}

s_u = function(v_mt, tau, k){ # v_mt is the moments of the k'th component, should be a length n list of 1*d_i*2 array
  t_sums <- sum(sapply(v_mt$pos, function(x) sum(x[k,,2])))
  return(1/(t_sums*tau)^0.5)
}

s_v = function(u_mt, tau){ # u_mt just input the k'th rank component
  return(1/(sum(u_mt[,2])*tau)^(0.5))
}

fnorm = function(X,Y=NULL){ # find the 2 norm of entries stored in two lists of matrices
  if (is.null(Y)){return(sum(unlist(X)^2)/length(unlist(X)))}
  return(sum((unlist(X)-unlist(Y))^2)/length(unlist(X)))
}

res = function(X, u_mt, v_mt, k){ # a function to remove the effects of non-k'th component
  n = length(X); R_k = vector(mode="list",length=n); R = dim(u_mt$pos)[2]
  if (dim(u_mt$pos)[2] == 1){return(X)}
  for (i in 1:n){
    R_k[[i]] = X[[i]] - u_mt$pos[,-k,1]%*%matrix(v_mt$pos[[i]][-k,,1], nrow=R-1)
  }
  return(R_k)
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

update_u_factor = function(u_mt, hold_u, k){
  u_mt$par[[k]] = hold_u$par
  u_mt$pos[,k,] = hold_u$pos
  u_mt$loglik[k] = hold_u$loglik
  u_mt$neg_KL[k] = hold_u$neg_KL
  return(u_mt)
}

update_v_factor = function(v_mt, hold_v, k){
  v_mt$par[[k]] = hold_v$par
  for (i in 1:length(v_mt$pos)){
    v_mt$pos[[i]][k,,] = t(hold_v$pos[[i]])
  }
  v_mt$loglik[k] = hold_v$loglik
  v_mt$neg_KL[k] = hold_v$neg_KL
  return(v_mt)
}

irrEBNM = function(X,s,D,fn = c("EC","RBF","Season")){
  # s is just one number now since noise variance is same across timepoints
  # X should be a list of vectors (potentially) of varying length
  # D is also a list of vectors, which contains the corresponding timepoints of the vectors in X
  n = length(X)
  results = list(); results$prior = fn
  if (fn == "EC"){ # equicorrelation model (compound symmetry)
    # Starting values
    start_params <- c(eta = 1, rho = 0.5)
    # Bounds
    lower <- c(1e-5, 1e-5); upper <- c(Inf, 1-1e-5)
    # Objective function (note: optimization usually minimizes)
    obj_fun <- function(par) {
      avg_loglik_ec(par = par, X = X, s = s, D = D)
    }
    # Optimize
    res. <- optimx(start_params, obj_fun, method = "L-BFGS-B", 
                     lower = lower, upper = upper)
    eta_hat2 = max(res.$eta,1e-50)
    rho_hat = min(0.99,max(res.$rho,1e-50))
    if (eta_hat2 < 1e-49){rho_hat = 0}
    results$loglik = -n*res.$value
    results$par = c(eta_hat2^0.5, rho_hat)
    names(results$par) = c("eta","rho")
    # calculate posterior point estimates
    pos = vector(mode="list", length=n)
    for (i in 1:n){
      di = length(D[[i]])
      S_inv = solve(diag(s^2, nrow=di))
      V = (matrix(rho_hat,nrow=di,ncol=di) + diag(1-rho_hat, nrow=di))%*%diag(eta_hat2,nrow=di)
      V_inv = svd.inverse(V)
      pos[[i]] = array(dim = c(2,di))
      pos[[i]][1,] = svd.inverse(S_inv+V_inv)%*%(S_inv%*%t(X[[i]]))
      pos[[i]][2,] = diag(pos[[i]][1,] %*% t(pos[[i]][1,]) + svd.inverse(S_inv + V_inv))
    }
    results$pos = pos
    # Calculate the negative KL divergence that will be used in objective function calculation
    term = 0
    for (i in 1:n){
      di = length(D[[i]])
      S = diag(s^2, nrow=di); S_inv = solve(S)
      term = term + di/2*log(2*pi) + 1/2*sum(diag(log(S))) + 1/2*(X[[i]]%*%S_inv%*%t(X[[i]]) - 2*X[[i]]%*%S_inv%*%pos[[i]][1,] + sum(1/s^2*pos[[i]][2,]))
    }
    results$neg_KL = term + results$loglik
  }else if (fn%in%c("RBF","Season")){ # Gaussian Process with RBF kernel
    par_init = log(GridSearchRBF(X=X,s=s,D, fine=0.1, range_max = 1)) # log(c(0.5,0.5)) # 
    #if (exp(par_init[1]) == 0.1){par_init[1] = log(GridSearchRBF(X,s,D,fine=0.01,range_max = 0.2)[1])}
    # Optimize
    
    res. = optimx(par=par_init,fn=avg_loglik_rbf,gr=grad,method="Nelder-Mead",X=X,s=s,D=D) # nlm(f=avg_loglik_rbf, p=par_init, X=X, s=s, D=D, iterlim=5) # minimize the negative loglik 
    eta_hat = max(exp(res.$p1),1e-50) # max(exp(res.$estimate[1]),1e-50) # 
    alp_hat = exp(res.$p2) # max(exp(res.$estimate[2]),0) # 
    #V <- eta_hat^2 * exp(-pairwise_diffs / (2*alp_hat^2))
    results$loglik = -n*res.$value # -n*res.$minimum # 
    results$par = c(eta_hat, alp_hat)
    names(results$par) = c("eta","alpha")
    # calculate posterior point estimates
    pos = vector(mode="list", length=n)
    for (i in 1:n){
      S_inv = solve(diag(s^2, nrow=length(D[[i]])))
      pairwise_diffs = outer(D[[i]], D[[i]], function(a, b) (a-b)^2)
      V = eta_hat^2 * exp(-pairwise_diffs / (2*alp_hat^2))
      V_inv = svd.inverse(V)
      pos[[i]] = array(dim = c(2,length(D[[i]])))
      pos[[i]][1,] = svd.inverse(S_inv+V_inv)%*%(S_inv%*%t(X[[i]]))
      pos[[i]][2,] = diag(pos[[i]][1,] %*% t(pos[[i]][1,]) + svd.inverse(S_inv + V_inv))
    }
    results$pos = pos
    # Calculate the negative KL divergence that will be used in objective function calculation
    term = 0
    for (i in 1:n){
      d = length(D[[i]])
      S = diag(s^2, nrow=d); S_inv = solve(S)
      term = term + d/2*log(2*pi) + 1/2*sum(diag(log(S))) + 1/2*(X[[i]]%*%S_inv%*%t(X[[i]]) - 2*X[[i]]%*%S_inv%*%pos[[i]][1,] + sum(1/s^2*pos[[i]][2,]))
    }
    results$neg_KL = term + results$loglik
  }
  return(results)
}

get_F = function(X, D, u_mt, v_mt, tau){
  p = nrow(X[[1]])
  R_2 = get_R(X, D, u_mt, v_mt)
  loglik = -length(unlist(D))*p/2*log(2*pi/tau) - tau/2*sum(unlist(R_2))
  elbo = loglik + sum(u_mt$neg_KL) + sum(v_mt$neg_KL)
  return(elbo)
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

init_factors = function(X,D,D_init=NULL,R){
  if (!is.null(D_init)){D=D_init}
  p = nrow(X[[1]]); n = length(X)
  X_unfold = do.call(cbind, X)
  
  u_mt = v_mt = list()
  svd_init = svd(X_unfold,nu=R,nv=R) # softImpute(X_unfold, rank.max=R) # 
  u = as.matrix(svd_init$u,ncol=R); s_vals = sqrt(svd_init$d); v = as.matrix(svd_init$v,ncol=R)
  u_mt$pos = array(c(u%*%diag(s_vals,nrow=R),u^2%*%diag(s_vals^2,nrow=R)), dim=c(p,R,2))
  pre_v_mt = array(c(v%*%diag(s_vals,nrow=R),v^2%*%diag(s_vals^2,nrow=R)), dim=c(nrow(v),R,2))
  
  v_mt$pos = vector(mode="list",length=n)
  ends <- cumsum(lengths(D))
  starts <- c(1, head(ends, -1) + 1)
  D_idx <- mapply(`:`, starts, ends, SIMPLIFY = FALSE)
  for (i in 1:n){
    to_perm = array(pre_v_mt[D_idx[[i]],,],dim=c(lengths(D)[i],R,2))
    v_mt$pos[[i]] = array(aperm(to_perm,c(2,1,3)),dim=c(R,length(D[[i]]),2)) # each horizontal slice is d_i*2, representing one rank component of first & second moments
  } # each element is a rk*d_i*2 array
  
  u_mt$par = v_mt$par = list() # $par as length rk list (since each column may have > 1 prior parameters)
  return(list(u_mt,v_mt))
}

single_factor = function(X, D, D_init = NULL, fn=c("EC","RBF"), thres=1e-8, max_iter = 50){
  n = length(X); p = nrow(X[[1]]); R=1
  init = init_factors(X,D,D_init=D_init,R=1)
  u_mt = init[[1]]
  v_mt = init[[2]]
  S_new = S_old = vector(mode="list",length=n)
  for (i in 1:n){
    S_new[[i]] = u_mt$pos[,,1] %*% matrix(v_mt$pos[[i]][,,1], nrow=1)
    S_old[[i]] = array(-Inf,dim=c(p,length(D[[i]])))
  }
  F_vec=c(); count = 0
  while(fnorm(S_new,S_old)>thres & count < max_iter){
    S_old = S_new
    R_2 = get_R(X,D,u_mt,v_mt)
    tau = p*length(unlist(D))/sum(unlist(R_2))
    for (k in 1:R){
      R_k = res(X, u_mt, v_mt, k)
      hold_u = ebnm_gaussian(x = x_u(R_k, v_mt, tau, k), s = s_u(v_mt,tau,k))
      u_mt = update_u_factor(u_mt,hold_u,k)
      hold_v = irrEBNM(X = x_v(R_k, u_mt$pos[,k,], tau), s = s_v(u_mt$pos[,k,], tau), D, fn=fn)
      v_mt = update_v_factor(v_mt,hold_v,k)
    }
    F_vec = c(F_vec,get_F(X, D, u_mt, v_mt, tau))
    S_new = lapply(v_mt$pos, function(v) u_mt$pos[,,1] %*% matrix(v[,,1], nrow = R))
    count=count+1
    print(fnorm(S_new,S_old))
  }
  results = list(u = u_mt, v = v_mt, elbo = F_vec, S = S_new, noise = tau, n_iter = count)
  return(results)
}

backfit = function(X, D, D_init = NULL, R, fn = c("EC","Free","RBF"), 
                   u_mt = NULL, v_mt = NULL, thres=1e-8, 
                   diagnose = F, max_iter = 50, null_check=T){
  n = length(X); p = nrow(X[[1]])
  if (is.null(u_mt) | is.null(v_mt)){
    init = init_factors(X,D,D_init=D_init,R=R)
    u_mt = init[[1]]; v_mt = init[[2]]
  }
  S_new = S_old = vector(mode="list",length=n)
  for (i in 1:n){
    S_new[[i]] = u_mt$pos[,,1] %*% matrix(v_mt$pos[[i]][,,1], nrow=R)
    S_old[[i]] = array(-Inf,dim=c(p,length(D[[i]])))
  }
  F_vec=c(); count = 0
  while(fnorm(S_new,S_old)>thres & count < max_iter){
    S_old = S_new
    R_2 = get_R(X,D,u_mt,v_mt)
    tau = p*length(unlist(D))/sum(unlist(R_2))
    for (k in 1:R){
      R_k = res(X, u_mt, v_mt, k)
      hold_u = ebnm_gaussian(x = x_u(R_k, v_mt, tau, k), s = s_u(v_mt,tau,k))
      u_mt = update_u_factor(u_mt,hold_u,k)
      hold_v = irrEBNM(X = x_v(R_k, u_mt$pos[,k,], tau), s = s_v(u_mt$pos[,k,], tau), D, fn=fn)
      v_mt = update_v_factor(v_mt,hold_v,k)
    }
    F_vec = c(F_vec,get_F(X, D, u_mt, v_mt, tau))
    S_new = lapply(v_mt$pos, function(v) u_mt$pos[,,1] %*% matrix(v[,,1], nrow = R))
    count=count+1
    if (diagnose){print(fnorm(S_new,S_old))}
  }
  results = list(u = u_mt, v = v_mt, elbo = F_vec, S = S_new, noise = tau, n_factors = R)
  if (null_check){results = nullcheck(results,D)}
  return(results)
}

nullcheck = function(results,D){
  n = length(results$S); p = nrow(results$S[[1]])
  r = dim(results$u$pos)[2]
  if (any(colSums(abs(matrix(results$u$pos[,,1],ncol=r)))<1e-8)){
    ind = which(colSums(abs(matrix(results$u$pos[,,1],ncol=r)))<1e-8)
    results$u$pos = array(results$u$pos[,-ind,], dim=c(p,r-length(ind),2))
    for (i in 1:n){
      results$v$pos[[i]] = array(results$v$pos[[i]][-ind,,],dim=c((r-length(ind)),length(D[[i]]),2))
    }
    results$u$par = results$u$par[-ind]
    results$v$par = results$v$par[-ind]
    results$u$loglik = results$u$loglik[-ind]
    results$v$loglik = results$v$loglik[-ind]
    results$u$neg_KL = results$u$neg_KL[-ind]
    results$v$neg_KL = results$v$neg_KL[-ind]
    results$n_factors = results$n_factors - length(ind)
  }
  return(results)
}

greedy_fac = function(X, D, D_init = NULL, fn=c("EC","RBF"), thres=1e-8){
  n = length(X); p = nrow(X[[1]])
  u_mt = v_mt = list()
  u_mt$pos = array(0,dim=c(p,min(n,p),2))
  v_mt$pos = lapply(D, function(d_i) array(0, dim = c(min(n,p),length(d_i),2)))
  u_mt$par = v_mt$par = list()
  S_new = lapply(v_mt$pos, function(v) u_mt$pos[,,1] %*% matrix(v[,,1], nrow = min(n,p)))
  k=0
  done=F
  F_vec = c()
  while (!done){
    R_k = lapply(1:n, function(i) X[[i]] - S_new[[i]])
    new = single_factor(R_k, D=D, D_init = D_init, fn=fn, thres=thres)
    if (sum(abs(new$u$pos[,,1])) < 1e-10){break}
    k = k+1
    u_mt$pos[,k,] = new$u$pos
    u_mt$par[[k]] = new$u$par[[1]]
    u_mt$loglik[k] = new$u$loglik
    u_mt$neg_KL[k] = new$u$neg_KL
    for (i in 1:n){
      v_mt$pos[[i]][k,,] = new$v$pos[[i]][1,,]
    }
    v_mt$par[[k]] = new$v$par[[1]]
    v_mt$loglik[k] = new$v$loglik
    v_mt$neg_KL[k] = new$v$neg_KL
    tau = new$noise
    if (k>1){
      if (get_F(X, D, u_mt, v_mt, tau) <= F_vec[k-1]){
        k = k-1
        u_mt$par = u_mt$par[1:k]; v_mt$par = v_mt$par[1:k]
        u_mt$loglik = u_mt$loglik[1:k]; v_mt$loglik = v_mt$loglik[1:k]
        u_mt$neg_KL = u_mt$neg_KL[1:k]; v_mt$neg_KL = v_mt$neg_KL[1:k]
        break
      }
    }
    F_vec[k] = get_F(X, D, u_mt, v_mt, tau)
    S_new = lapply(v_mt$pos, function(v) u_mt$pos[,,1] %*% matrix(v[,,1], nrow = min(n,p)))
    print(k)
  }
  u_mt$pos = array(u_mt$pos[,1:k,], dim=c(p,k,2))
  v_mt$pos = lapply(1:n, function(i) array(v_mt$pos[[i]][1:k,,], dim=c(k,length(D[[i]]),2)))
  return(list(u = u_mt, v = v_mt, elbo = F_vec, S = S_new, noise = tau, n_factors = k))
}


split_back <- function(big_mat, col_counts) {
  ends <- cumsum(col_counts)                # end positions
  starts <- c(1, head(ends, -1) + 1)        # start positions
  out <- mapply(function(s, e) big_mat[, s:e, drop = FALSE],
                starts, ends,
                SIMPLIFY = FALSE)
  out
}

full_impute = function(X,D,D_init=NULL,R=NULL,fn=c("EC","RBF","Free"), 
                       method=c("greedy+backfit","backfit"),
                       thres=1e-5,max_iter=100,diagnose=F,null_check=T){
  # find missing observations using the unfolded representation, fold it back to perform algorithm, then unfold again to impute
  n = length(X); p = nrow(X[[1]])
  X_unfold = do.call(cbind, X)
  M = which(is.na(X_unfold))
  count = 0
  X_unfold[M] = 0
  X_imp = split_back(X_unfold,lengths(D))
  if (method=="greedy+backfit"){
    g_init = greedy_fac(X_imp,D,D_init,fn=fn,thres=thres)
    R=g_init$n_factors; u_mt=g_init$u; v_mt=g_init$v
  }else {
    if (is.null(R)){return("n factors required for backfit only algorithm")}
    init = init_factors(X_imp,D,D_init,R=R)
    u_mt = init[[1]];v_mt = init[[2]]
  }
  S_new = S_old = vector(mode="list",length=n)
  for (i in 1:n){
    S_new[[i]] = u_mt$pos[,,1] %*% matrix(v_mt$pos[[i]][,,1], nrow=R)
    S_old[[i]] = array(-Inf,dim=c(p,length(D[[i]])))
  }
  F_vec=c(); count = 0
  while(fnorm(S_new,S_old)>thres & count < max_iter){
    S_old = S_new
    R_2 = get_R(X_imp,D,u_mt,v_mt)
    tau = p*length(unlist(D))/sum(unlist(R_2))
    for (k in 1:R){
      R_k = res(X_imp, u_mt, v_mt, k)
      hold_u = ebnm_gaussian(x = x_u(R_k, v_mt, tau, k), s = s_u(v_mt,tau,k))
      u_mt = update_u_factor(u_mt,hold_u,k)
      hold_v = irrEBNM(X = x_v(R_k, u_mt$pos[,k,], tau), s = s_v(u_mt$pos[,k,], tau), D, fn=fn)
      v_mt = update_v_factor(v_mt,hold_v,k)
    }
    F_vec = c(F_vec,get_F(X_imp, D, u_mt, v_mt, tau))
    S_new = lapply(v_mt$pos, function(v) u_mt$pos[,,1] %*% matrix(v[,,1], nrow = R))
    X_unfold[M] = do.call(cbind,S_new)[M] # replace missing entries with new imputed values
    X_imp = split_back(X_unfold,lengths(D)) # fold back data into list form
    count=count+1
    if (diagnose){print(fnorm(S_new,S_old))}
  }
  results = list(u = u_mt, v = v_mt, elbo = F_vec, S = S_new, noise = tau, n_factors = R)
  if (null_check){results = nullcheck(results,D)}
  return(results)
}
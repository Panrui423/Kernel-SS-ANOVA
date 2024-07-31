### AIC BIC
AICBIC = function(M, nlambda, s, Y.all, hat_f_m_p)
{
  N = length(Y.all)
  M <- M + nlambda * diag(N)
  M_inv = solve(M)
  Alam = diag(N) - (nlambda)*(M_inv - M_inv %*% s %*% solve(t(s) %*% M_inv %*% s) %*% t(s) %*% M_inv)
  
  res = c(N * log(mean((Y.all - hat_f_m_p)^2)) + 2  * sum(diag(Alam)), 
          N * log(mean((Y.all - hat_f_m_p)^2)) + log(N)  * sum(diag(Alam)), sum(diag(Alam)))
  
  res
}


## Selection of Representative Meta
select_rep <- function(N, Y, pro, uniq = TRUE) {   
  nq0 = floor(pro * N)+1  
  nk = 20  
  rang = range(Y)+c(-0.01,+0.01)
  interv = seq(rang[1], rang[2], length=nk+1)
  ninterv = NULL
  obsinterv = list()
  for(i in 1:nk){
    m0 = interv[i]
    m1 = interv[i+1]
    loc = which(Y >= m0 & Y < m1)
    ninterv = c(ninterv, length(loc))
    obsinterv[[i]] = loc
  }
  
  nnum = floor(ninterv/N*nq0)
  nnum[nnum == 0] = 1
  nq = sum(nnum)
  
  id.basis=NULL
  for(i in 1:nk){
    loc = obsinterv[[i]]
    id.basis = c(id.basis, loc[1:nnum[i]])
  }
  
  if(uniq) {
    id.basis = sort(unique(id.basis))
  } else{
    id.basis = sort(id.basis)
  }
  N.gss = length(id.basis)
  list_se = list()
  list_se[[1]] = id.basis
  list_se[[2]] = N.gss
  return(list_se)
}

## Tensor product function
Tensor_f <- function(K) {
  
  N = dim(K)[1]
  Nstar = dim(K)[2]
  pvar = dim(K)[3]
  
  Q <- array(0, dim = c(N, Nstar, pvar + pvar * (pvar - 1)/2 ) )
  Q[,,1:pvar] = K
  len = pvar
  
  for(i in 1:(pvar-1)) {
    for(j in (i+1):(pvar)){
      len = len + 1
      Q[, , len ] = K[,,i] * K[,,j]
    }
  }
  return(Q)
}


### bootstrap_sample
bootstrap_sample <- function(data, R) {
  samples <- matrix(NA, ncol = R, nrow = dim(data)[1])
  for (i in 1:R) {
    samples[, i] <- sample(data, replace = TRUE)
  }
  return(samples)
}




## known id.basis, estimate other parts
ker_est_all <- function(N, dist_train, id.basis, N.gss, pvar,h) {
  Q_rep_est <- array(0,dim = c(N.gss, N.gss, pvar)) 
  for (k in 1:pvar) {
    Q_rep_est[,,k] <- MetricKern(dist_train[id.basis, id.basis, k]/h[k])
  }
  
  dist_c_est <- array(0, dim = c(N, N.gss, pvar))
  for (i in 1:pvar) {
    dist_c_est[, , i] <- rbind(dist_train[id.basis, id.basis, i],
                               dist_train[-id.basis, id.basis, i])
  }
  
  ker_patch_est <- array(0, dim = c(N-N.gss, N.gss, pvar))
  ker_patch_r <- array(0, dim = c(N-N.gss, N-N.gss, pvar))
  for (k in 1:pvar) {
    ker_patch_est[,,k] <- MetricKtest(dist_c_est[,,k]/h[k], N.gss, N-N.gss)
    ker_patch_r[,,k] <- MetricKest_general(dist_c_est[,,k]/h[k], N.gss, N-N.gss)  # F test need 
  }
  
  Ker4_est <- array(0, dim = c(N, N, pvar))
  Ker4_est[id.basis, id.basis, ] <- Q_rep_est
  Ker4_est[-id.basis, id.basis, ] <- ker_patch_est
  for(i in 1:pvar)
    Ker4_est[id.basis, -id.basis, i] <- t(ker_patch_est[,,i])
  
  Ker4_est[-id.basis, -id.basis, ] <- ker_patch_r
  
  ## Tensor product
  Q_all_est = Tensor_f(Ker4_est)
  
  return(Q_all_est)
  
}


Ker_est_rep <- function(N, dist_train, rep.num, id_list, pvar, h){
  Q_rep <- array(0,dim = c(rep.num, rep.num, pvar)) 
  for (k in 1:pvar) {
    Q_rep[,,k] <- MetricKern(dist_train[id_list, id_list, k]/h[k])
  }
  
  dist_c <- array(0, dim = c(N, rep.num, pvar))
  for (i in 1:pvar) {
    dist_c[, , i] <- rbind(dist_train[id_list, id_list, i],
                           dist_train[-id_list, id_list, i])
  }
  
  ker_patch <- array(0, dim = c(N-rep.num, rep.num, pvar))
  for (k in 1:pvar) {
    ker_patch[,,k] <- MetricKtest(dist_c[,,k]/h[k], rep.num, N-rep.num)
  }
  
  Ker4 <- array(0, dim = c(N, N, pvar))
  Ker4[id_list, id_list, ] <- Q_rep
  Ker4[-id_list, id_list, ] <- ker_patch
  
  #for(i in 1:pvar)
  #  Ker4[id_list, -id_list, i] <- t(ker_patch[,,i])
  
  ## Tensor product
  Q_rep_est = Tensor_f(Ker4)
  
  return(Q_rep_est)
}

computeSigma <- function(theta, kernel_array, pall, N){
  h_K_tmp <- matrix(0, N, N)
  for (i in 1:pall) {
    h_K_tmp <- 10^theta[i] * kernel_array[, , i] + h_K_tmp
  }
  return(h_K_tmp)
}


kernalGaussian <- function(s,t,para){
  h <- length(s)*para
  kernalX <-exp(-sum((s-t)^2)/(2*h^2))
  return(kernalX)
}

kernalLap <- function(s, t, para)
{
  p <- length(s)
  h <- p*para
  kernalX <- exp(-sum((s - t)^2)^0.5/(2*h)) 
  
  return(kernalX)
}







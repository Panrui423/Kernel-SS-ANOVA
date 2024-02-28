library(gss)
library(stringr)
library(ggplot2)
library(doMC)
library(foreach)
library(ggpubr)
library(tidyverse)
library(SpatialPack)
library(Rmisc)

Rcpp::sourceCpp("Kernelmatrix.cpp")



####################################

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
    loc = which(Y >= m0 & Y <= m1)
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




######################################



load("ADNIdatatrain.Rdata")
load("ADNIdatatest.Rdata")

Y.train = log(data.train$Score + 2) # logrithm transformation 
N = length(Y.train)
xfun0.train = data.train$Xfun
s.train = data.train$Xscalar  # c(Intercept,  gender, Edu, age);  normalized age, Edu


Y.test = log(data.test$Score + 2)
xfun0.test = data.test$Xfun
s.test = data.test$Xscalar


Y.all = c(c(Y.train),c(Y.test))
s = rbind(s.train, s.test)
Xfun <- rbind(xfun0.train, xfun0.test)

m = length(Y.all)

label_left = read.csv('label-left-12.csv')
label_right = read.csv('label-right-12.csv')

# x$b[parasubiculum_203] <- 1
# x$b[presubiculum_204] <- 2
# x$b[subiculum_205] <- 3
# x$b[CA1_206] <- 4
# x$b[CA3_208] <- 5
# x$b[CA4_209] <- 6
# x$b[GC_DG_210] <- 7
# x$b[molecular_layer_HP_214] <- 8
# x$b[HATA_211] <- 9
# x$b[fimbria_212] <- 10
# x$b[HP_tail_226] <- 11
# x$b[hippocampal_fissure_215] <- 12


c_reg = c(7, 5, 4, 3) 
c_l = length(c_reg)

dist_matrix4 <- array(0,dim = c(m, m, c_l))
for (i in 1:c_l) {
  dist_matrix4[,,i] <- as.matrix(dist(Xfun[,15000 + which(label_right[,2] == c_reg[i])], 
                                      diag = T, upper = T, method = 'euclidean')) 
}

########### Parameter selection ##########
set.seed(1)
N = length(Y.all)
Y = Y.all
pho.r = 1/4
k.rep = 1
esp = 0.001
pvar = dim(dist_matrix4)[3]
pall = pvar + pvar*(pvar-1)/2

sel_res = select_rep(N, Y, pho.r)
id.basis = sel_res[[1]]
N.gss = sel_res[[2]]

##### Fitting 
h = rep(1, pvar)


dist_train <- dist_matrix4

Q_rep <- array(0,dim = c(N.gss, N.gss, pvar)) 
for (k in 1:pvar) {
  Q_rep[,,k] <- MetricKern(dist_train[id.basis, id.basis, k]/h[k])
}

dist_c <- array(0, dim = c(N, N.gss, pvar))
for (i in 1:pvar) {
  dist_c[, , i] <- rbind(dist_train[id.basis, id.basis, i],
                         dist_train[-id.basis, id.basis, i])
}

ker_patch <- array(0, dim = c(N-N.gss, N.gss, pvar))
for (k in 1:pvar) {
  ker_patch[,,k] <- MetricKtest(dist_c[,,k]/h[k], N.gss, N-N.gss)
}

Ker4 <- array(0, dim = c(N, N, pvar))
Ker4[id.basis, id.basis, ] <- Q_rep
Ker4[-id.basis, id.basis, ] <- ker_patch

## Tensor product
Q_m4 <- Tensor_f(Ker4)

s1 = s
fit_m4 <- gss:::mspreg1(s = s, r = Q_m4[,id.basis,], id.basis=id.basis, method = 'm',
                        y = Y.all,wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE)

theta_m4 <- fit_m4$theta
c_m4 <- matrix(0, N, 1)
c_m4[id.basis,1] <- fit_m4$c
d_m4 <- fit_m4$d
d_m4

##  fitting
h_K_m <- matrix(0, N, N)
for (i in 1:pall) {
  h_K_m <- 10^theta_m4[i] * Q_m4[, , i] + h_K_m
}

d_m4 <- as.matrix(d_m4) 
round(d_m4, 4)
hat_f_m4 <- h_K_m %*% c_m4
hat_f_m_p <- h_K_m %*% c_m4 + s %*% d_m4  # hat_f_m_p <- h_K_m %*% c_m4 + s %*% d_m4

## MSE
round(mean(((hat_f_m_p - Y.all))^2), 4)


#### Test whether the scalar variable is significant ####

nlambda <- 10^fit_m4$nlambda
M <- matrix(0, N, N)
for (i in 1:pall) {
  M = 10^theta_m4[i] * Q_m4[,,i] + M
}

M <- M + nlambda * diag(N)
sigma_hat <- var(Y.all - hat_f_m4)
cov_d <- solve(t(s1) %*% solve(M) %*% s1) %*% t(s1) %*% solve(M) %*% solve(M) %*%
  s1 %*% solve(t(s1) %*% solve(M) %*% s1) * sigma_hat[1,1]
SDest <- round((diag(cov_d))^(0.5),4)
SDest
# round(2*pnorm(-abs(d_m4/SDest)), 4)


## cosine diagnosis
Cos_diag <- matrix(0, 1, pall)
for (i in 1:pall) {
  Cos_diag[1, i] <- sum((Q_m4[,,i] %*% c_m4 * 10^theta_m4[i]) * hat_f_m4) / sum(hat_f_m4^2)
}
Cos_diag <- round(Cos_diag, 4)
Cos_diag

## single fit ##
fit_fs <- matrix(0, N, 3)
sigSeq <- c(7, 9)
siglen <- length(sigSeq)
for (i in 1:siglen) {
  fit_fs[,i] = (10^theta_m4[sigSeq[i]] * Q_m4[, , sigSeq[i]]) %*% c_m4
}

for (i in c(1:pall)) {
  fit_fs[,siglen+1] = (10^theta_m4[i] * Q_m4[, , i]) %*% c_m4 + fit_fs[,siglen+1]
}


###### MRD ########


n.fit <- length(Y.all)

## fit MRD

fitMRD.1 <- round(mean(((hat_f_m_p - Y.all)/Y.all)^2), 4)
fitMRD.1








######################## Simulation hippocampus data ######################

f.real <- matrix(0, m, 1)
M <- matrix(0, m, m)
theta.real <- c(0.5, 0.5, 0.8, 0.8, 0, 0, 4, 0, 4, 0)
d.real <- d_m4
set.seed(1)
c.real <- c_m4

for (i in c(1:pall)) {
  M <- theta.real[i] * Q_m4[,,i] + M
}

f.real.rm <- c(M %*% c.real)
f.real <- c(M %*% c.real + s %*% d.real)

## real cosine diagnosis

Cos_diag_sim <- matrix(0, 1, pall)
for (i in 1:pall) {
  Cos_diag_sim[1, i] <- sum((Q_m4[,,i] %*% c.real * theta.real[i]) * f.real.rm) / sum(f.real.rm^2)
}
Cos_diag_sim <- round(Cos_diag_sim, 4)
Cos_diag_sim



######################### cross validation ###############################
gcv <- function(nlambda, theta, y, n, Q_fit, Z) {
  Sigma <- matrix(0, n, n)
  for (i in 1:dim(Q_fit)[3]) {
    Sigma <- Q_fit[,,i] * 10^theta[i] + Sigma
  } 
  M <- Sigma + 0.5 * 10^nlambda * diag(n)
  H <- 10^nlambda * solve(M) %*% (diag(n) - Z %*% solve(t(Z) %*% solve(M) %*% Z) %*% t(Z) %*% solve(M))
  gcv_res <- ((t(y) %*% H %*% H %*% y) / n) / (sum(diag(H)) / n)^2
  return(gcv_res)
}

fitFun_cv <- function(sd){
  
  set.error <- rnorm(m, 0, sd)
  Y.real <- c(f.real) + set.error
  N = length(Y.real)
  
  ######################## Selection of representative samples (split between training and testing) #######################
  
  res.j = select_rep(N, Y.real, pho.r)
  id.j = res.j[[1]]
  N.j = res.j[[2]]
  
  ######
  
  dist_train <- dist_matrix4
  
  Q_1 <- array(0, dim = c(N.j, N.j, pvar))
  for (k in 1:pvar) {
    Q_1[,,k] <- MetricKern(dist_train[id.j, id.j, k]/h[k])  # Bandwidth option 1
  }
  
  d_1 <- array(0, dim = c(N, N.j, pvar))
  for (i in 1:pvar) {
    d_1[, , i] <- rbind(dist_train[id.j, id.j, i],
                        dist_train[-id.j, id.j, i])
  }
  
  Q_patch <- array(0, dim = c(N - N.j, N.j, pvar))
  for (k in 1:pvar) {
    Q_patch[,,k] <- MetricKtest(d_1[,,k]/h[k], N.j, N - N.j)
  }
  
  Kermat_train <- array(0, dim = c(N, N, pall))
  Kermat_train[id.j, id.j, 1:pvar] <- Q_1
  Kermat_train[-id.j, id.j, 1:pvar] <- Q_patch
  
  ##### Tensor product #####
  
  Kermat_train <- Tensor_f(Kermat_train[,,1:pvar])
  
  fit_Tensor <- gss:::mspreg1(s = s, r = Kermat_train[, id.j, ], id.basis = id.j, method = 'm',
                              y = Y.real, wt = 1, alpha = 1.4, random = NULL, varht = 1, skip.iter=FALSE)
  
  theta_fit <- fit_Tensor$theta
  c_fit <- matrix(0, N, 1)
  c_fit[id.j,1] <- fit_Tensor$c
  d_fit <- fit_Tensor$d
  
  ###  Fitting
  MSE_01 = c()
  h_f <- matrix(0, N, 1)
  for (i in 1:siglen) {
    h_f <- (10^theta_fit[sigSeq[i]] * Kermat_train[, , sigSeq[i]]) %*% c_fit
    MSE_01[i] <- mean((h_f - fit_fs[,i])^2)
  }
  
  h_f <- matrix(0, N, 1)
  for (i in c(1:pall)) {
    h_f <- (10^theta_fit[i] * Kermat_train[, , i]) %*% c_fit + h_f
  }
  
  MSE_01[siglen+1] <- mean((h_f - fit_fs[,siglen+1])^2)
  
  ## cosine diagnosis
  h_K <- matrix(0, N, N)
  for (i in 1:pall) {
    h_K <- 10^theta_fit[i] * Kermat_train[, , i] + h_K
  }
  
  hat_f <- h_K %*% c_fit
  
  Cos_diag_1 <- matrix(0, 1, pall)
  for (i in 1:pall) {
    Cos_diag_1[1, i] <- sum((Kermat_train[,,i] %*% c_fit * 10^theta_fit[i]) * hat_f) / sum(hat_f^2)
  }
  Cos_diag_1 <- round(Cos_diag_1, 4)
  
  ######### Gaussian #############
  
  
  MSE_02_all <- list()
  Cos_diag_2_all <- list()
  gcv_list <- list()
  
  h_list <- c(0.01, 0.04, 0.08, 0.25, 0.5)
  for (h in 1:length(h_list)) {
    Q_fit <- array(0, dim = c(m, m, pall))
    for (k in 1:length(c_reg)) {
      for (i in 1:m) {
        for (j in 1:m) {
          Q_fit[i,j, k] = kernalGaussian(Xfun[i, 15000 + which(label_right[,2] == c_reg[k])], 
                                         Xfun[j, 15000 + which(label_right[,2] == c_reg[k])],para = h_list[h])
        }
      }
    }
    
    Q_fit <- Tensor_f(Q_fit[, , 1:pvar])
    
    fit_gauss.j <- gss:::mspreg1(s = s, r = Q_fit, id.basis = 1:N, method = 'm',
                                 y = Y.real, wt = 1, alpha = 1.4, random = NULL, varht = 1, skip.iter = FALSE)
    
    theta_fit <- fit_gauss.j$theta
    c_fit <- fit_gauss.j$c
    d_fit <- fit_gauss.j$d
    
    ###  Fitting
    MSE_02 = c()
    h_f <- matrix(0, N, 1)
    for (i in 1:siglen) {
      h_f <- (10^theta_fit[sigSeq[i]] * Q_fit[, , sigSeq[i]]) %*% c_fit
      MSE_02[i] <- mean((h_f - fit_fs[,i])^2)
    }
    
    h_f <- matrix(0, N, 1)
    for (i in c(1:pall)) {
      h_f <- (10^theta_fit[i] * Q_fit[, , i]) %*% c_fit + h_f
    }
    
    MSE_02[siglen+1] <- mean((h_f - fit_fs[,siglen+1])^2)
    
    ## cosine diagnosis
    h_K <- matrix(0, N, N)
    for (i in 1:pall) {
      h_K <- 10^theta_fit[i] * Q_fit[, , i] + h_K
    }
    
    hat_f <- h_K %*% c_fit
    
    Cos_diag_2 <- matrix(0, 1, pall)
    for (i in 1:pall) {
      Cos_diag_2[1, i] <- sum((Q_fit[,,i] %*% c_fit * 10^theta_fit[i]) * hat_f) / sum(hat_f^2)
    }
    Cos_diag_2 <- round(Cos_diag_2, 4)
    MSE_02_all[[h]] <- MSE_02
    Cos_diag_2_all[[h]] <- Cos_diag_2
    gcv_list[[h]] <- gcv(fit_gauss.j$nlambda, theta_fit, Y.real, m, Q_fit, s)
  }
  min_gcv <- which.min(gcv_list)
  MSE_02 <- MSE_02_all[[min_gcv]]
  Cos_diag_2 <- Cos_diag_2_all[[min_gcv]]
  
  ######### Laplacian kernel #############
  
  MSE_03_all <- list()
  Cos_diag_3_all <- list()
  gcv_list <- list()
  
  h_list <- c(0.01, 0.04, 0.08, 0.25, 0.5)
  
  for (h in 1:length(h_list)) {
    Q_fit <- array(0, dim = c(m, m, pall))
    
    for (k in 1:length(c_reg)) {
      for (i in 1:m) {
        for (j in 1:m) {
          Q_fit[i,j, k] = kernalLap(Xfun[i, 15000 + which(label_right[,2] == c_reg[k])], 
                                    Xfun[j, 15000 + which(label_right[,2] == c_reg[k])],para = h_list[h])
        }
      }
    }
    
    Q_fit <- Tensor_f(Q_fit[, , 1:pvar])
    
    fit_lap <- gss:::mspreg1(s = s, r = Q_fit, id.basis = 1:N, method = 'm',
                             y = Y.real, wt = 1, alpha = 1.4, random = NULL, varht = 1, skip.iter = FALSE)
    
    theta_fit <- fit_lap$theta
    c_fit <- fit_lap$c
    d_fit <- fit_lap$d
    
    ###  Fitting
    MSE_03 = c()
    h_f <- matrix(0, N, 1)
    for (i in 1:siglen) {
      h_f <- (10^theta_fit[sigSeq[i]] * Q_fit[, , sigSeq[i]]) %*% c_fit
      MSE_03[i] <- mean((h_f - fit_fs[,i])^2)
    }
    
    h_f <- matrix(0, N, 1)
    for (i in c(1:pall)) {
      h_f <- (10^theta_fit[i] * Q_fit[, , i]) %*% c_fit + h_f
    }
    
    MSE_03[siglen+1] <- mean((h_f - fit_fs[,siglen+1])^2)
    
    ## cosine diagnosis
    
    h_K <- matrix(0, N, N)
    for (i in 1:pall) {
      h_K <- 10^theta_fit[i] * Q_fit[, , i] + h_K
    }
    
    hat_f <- h_K %*% c_fit
    
    Cos_diag_3 <- matrix(0, 1, pall)
    for (i in 1:pall) {
      Cos_diag_3[1, i] <- sum((Q_fit[,,i] %*% c_fit * 10^theta_fit[i]) * hat_f) / sum(hat_f^2)
    }
    Cos_diag_3 <- round(Cos_diag_3, 4)
    MSE_03_all[[h]] <- MSE_03
    Cos_diag_3_all[[h]] <- Cos_diag_3
    gcv_list[[h]] <- gcv(fit_lap$nlambda, theta_fit, Y.real, m, Q_fit, s)
  }
  min_gcv <- which.min(gcv_list)
  MSE_03 <- MSE_03_all[[min_gcv]]
  Cos_diag_3 <- Cos_diag_3_all[[min_gcv]]
  
  return(round(c(MSE_01, MSE_02, MSE_03, Cos_diag_1, Cos_diag_2, Cos_diag_3), 4))
  
}


set.seed(1)
registerDoMC(30)

result.NE2 = foreach(rep = 1:200, .combine = cbind) %dopar%{
  a <- fitFun_cv(0.05)
  b <- fitFun_cv(0.10)
  c <- fitFun_cv(0.15)
  c(a,b,c)
}


mean_left4 = rowMeans(result.NE2)
sd_left4 = apply(result.NE2, 1, sd)
mat = matrix(c(c(1:9), c(40:48), c(79:87)), 9, 3, byrow = TRUE)

mseRes <- data.frame(Sigma = matrix(c(0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.15, 0.15, 0.15), 9, 1),
                     K = rep(c("$mathcal{K}$", "$mathcal{K}_g$", "$mathcal{K}_l$"), 3),
                     
                     pi_7 = round(mean_left4[mat[,1]], 4)*100,
                     pi_7.sd = round(sd_left4[mat[,1]], 4)*100,
                     
                     pi_9 = round(mean_left4[mat[,2]], 4)*100,
                     pi_9.sd = round(sd_left4[mat[,2]], 4)*100,
                     
                     pi_all = round(mean_left4[mat[,3]], 4)*100,
                     pi_all.sd = round(sd_left4[mat[,3]], 4)*100,
                     
                     stringsAsFactors=FALSE)

library(xtable)
# The results in Table 3 were obtained
xtable(mseRes,digits=2)






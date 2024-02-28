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

#Y.train = Y.all
#s = Xscalar.all   


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


c_reg = c(1, 2, 3) 
c_l = length(c_reg)

dist_matrix4 <- array(0,dim = c(m, m, c_l))
for (i in 1:c_l) {
  dist_matrix4[,,i] <- as.matrix(dist(Xfun[,15000 + which(label_right[,2] == c_reg[i])], 
                                      diag = T, upper = T, method = 'euclidean')) 
}

########### 参数选择 ##########
# set.seed(3)
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

##### fit (data to a model)
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

#s1 = matrix(s[,1], length(s[,1]), 1)
s1 = s
fit_m4 <- gss:::mspreg1(s = s1, r = Q_m4[,id.basis,], id.basis=id.basis, method = 'm',
                        y = Y.all, wt=1, alpha = 1.4, random = NULL, varht=1, skip.iter=FALSE)

theta_m4 <- fit_m4$theta
c_m4 <- matrix(0, N, 1)
c_m4[id.basis,1] <- fit_m4$c
d_m4 <- fit_m4$d
d_m4

###  Fitting
h_K_m <- matrix(0, N, N)
for (i in 1:pall) {
  h_K_m <- 10^theta_m4[i] * Q_m4[, , i] + h_K_m
}

d_m4 <- as.matrix(d_m4) 
round(d_m4, 4)
hat_f_m4 <- h_K_m %*% c_m4
hat_f_m_p <- h_K_m %*% c_m4 + s1 %*% d_m4

## MSE
round(mean(((hat_f_m_p - Y.all))^2), 4)


#### Check if the scalar variable is significant ####

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


## cosine diagnosis
Cos_diag <- matrix(0, 1, pall)
for (i in 1:pall) {
  Cos_diag[1, i] <- sum((Q_m4[,,i] %*% c_m4 * 10^theta_m4[i]) * hat_f_m4) / sum(hat_f_m4^2)
}
Cos_diag <- round(Cos_diag, 4)
Cos_diag



n.fit <- length(Y.all)


## fit MRD

fitMRD.1 <- round(mean(((hat_f_m_p - Y.all)/Y.all)^2), 4)










library(gss)
library("doMC")
library("foreach")
source("Manifold_simulation.R")
source("Klein_simulation.R")

registerDoMC(20)
###########  #############
# Simulation needs to be parallel computation, 
# this paper uses 20 cores for parallel computation, 
# the reader can also according to their own needs 
# and hardware conditions to set the number of cores that need to be parallel computation
########################
########## sub-functions 

### Manifold data
simulation_M_rep <- function(n.fit, n.test, h, scale, sigma){
  
  set.seed(2022)
  
  n <- n.fit + n.test
  list_M <- Manifold_f(n)
  
  Manifold_data <- list_M$Manifold_data    ## Generation of 300 * 3 matrices
  
  X <- list_M$X 
  
  y_M <- c()
  for (i in 1:n) {
    y_M[i] <- f_m(Manifold_data[i,1], Manifold_data[i,2], Manifold_data[i,3])
  }
  
  scalar <- rnorm(n, 0,sd = 1)
  y_M <- y_M + scalar
  f_M <- y_M
  y_M <- y_M + rnorm(n, 0,sd = sigma) 
  y_M.test <- y_M[-c(1:n.fit)]
  
  dist.matrix <- array(0,dim = c(n, n, 6))
  
  dist.matrix[,,1] <- as.matrix(dist(X[,1:100], diag = T, upper = T, method = 'manhattan'))  
  dist.matrix[,,2] <- as.matrix(dist(X[,101:250], diag = T, upper = T, method = 'manhattan'))
  dist.matrix[,,3] <- as.matrix(dist(X[,251:450], diag = T, upper = T, method = 'manhattan'))
  
  dist.matrix[,,4] <- dist.matrix[,,1] + dist.matrix[,,2]
  dist.matrix[,,5] <- dist.matrix[,,1] + dist.matrix[,,3]
  dist.matrix[,,6] <- dist.matrix[,,2] + dist.matrix[,,3]
  
  result_rep = foreach(i.rep = 1:200, .combine = cbind) %dopar%
{ 
  ###################### Selection of representative samples ################################
  
  fit_error_m2 <- predict_error_m2 <- 0
  
  for (l in 1:length(scale)) 
    {
    N = n.fit
    nq0 = floor(scale[l] * N)+1
    nk = 20  
    
    Y = y_M[1:n.fit]
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
      loc = loc[sample(1:ninterv[i])]
      id.basis = c(id.basis, loc[1:nnum[i]])
    }
    
    id.basis = sort(unique(id.basis))
    N.gss = length(id.basis)
    
    ######
    dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
    dist.test <- dist.matrix[(n.fit+1):n, , ]
   
    Q_rep <- array(0, dim = c(N.gss, N.gss, 6))
    for (k in 1:6) {
      if(k==4){timeend_1 = Sys.time()}
      Q_rep[,,k] <- MetricKern(dist.train[id.basis,id.basis,k]/h)
    }
    
    d_change <- array(0, dim = c(n.fit, N.gss, 6))
    for (i in 1:6) {
      d_change[, , i] <- rbind(dist.train[id.basis, id.basis, i],
                               dist.train[-id.basis, id.basis, i])
    }
    
    Q_patch <- array(0, dim = c(N-N.gss, N.gss, 6))
    for (k in 1:6) {
      if(k==4){timeend_2 = Sys.time()}
      Q_patch[,,k] <- MetricKtest(d_change[,,k]/h, N.gss, N-N.gss)
    }
    
    Kermat_train <- array(0, dim = c(n.fit, n.fit, 6))
    Kermat_train[id.basis, id.basis, ] <- Q_rep
    Kermat_train[-id.basis, id.basis, ] <- Q_patch
    
    timestart2 <- Sys.time()
    Kermat_train_tensor <- Kermat_train
    
    Kermat_train_tensor[,,4] <- Kermat_train[,,1] * Kermat_train[,,2]
    Kermat_train_tensor[,,5] <- Kermat_train[,,1] * Kermat_train[,,3]
    Kermat_train_tensor[,,6] <- Kermat_train[,,2] * Kermat_train[,,3]
    
    s = matrix(scalar[1:n.fit], n.fit, 1)
    fit_ms_m <- tryCatch(gss:::mspreg1(s = s,r=Kermat_train_tensor[,id.basis,],id.basis=id.basis,method = 'm',
                                       y = y_M[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                         error=function(x) NA)
    if((any(is.na(fit_ms_m))))
    { 
      fit_error_m2[l] = NA
      predict_error_m2[l] =NA
    }else{
      
      theta_m <- fit_ms_m$theta
      c_m <- matrix(0, n.fit, 1)
      c_m[id.basis,1] <- fit_ms_m$c
      d_m <- fit_ms_m$d
      
      ## fitting error
      hat_m <- matrix(0, n.fit, 1)
      M <- matrix(0, n.fit, n.fit)
      
      for (i in 1:6) {
        M <- 10^theta_m[i] * Kermat_train_tensor[,,i] + M
      }
      
      hat_m <- c(M %*% c_m) + scalar[1:n.fit] * d_m
      fit_error_m2[l] <- mean((hat_m-f_M[1:n.fit])^2)
      
      ## prediction error
      
      dx.test <- array(0, dim = c(n.test + N.gss, N.gss, 6))
      for (i in 1:6) {
        dx.test[, , i] <- rbind(dist.train[id.basis, id.basis, i],
                                dist.test[, id.basis, i])
      }
      
      kermat <- array(0, dim = c(n.test, N.gss, 6))
      for (j in 1:6) {
        kermat[, , j] <- MetricKtest(dx.test[, , j]/h, N.gss, n.test)
      }
      
      kermat[,,4] <- kermat[,,1] * kermat[,,2]
      kermat[,,5] <- kermat[,,1] * kermat[,,3]
      kermat[,,6] <- kermat[,,2] * kermat[,,3]
      
      Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
      Kermat_test[,id.basis,] <- kermat
      
      p_y_m <- matrix(0, n.test, 1)
      p_K_m <- matrix(0, n.test, n.fit)
      for (i in 1:6) {
        p_K_m <- 10^theta_m[i] * Kermat_test[,,i] + p_K_m
      }
      p_y_m <- c(p_K_m %*% c_m) + scalar[(n.fit+1):n]*d_m
      predict_error_m2[l] <- mean((p_y_m-y_M.test)^2)
    }
  }
  
 
  c(predict_error_m2, fit_error_m2) 
}

  return(result_rep)
}

############## Klein Bottle

simulation_K_rep <- function(n.fit, n.test, h, scale, sigma){
  set.seed(2022)
  n <- n.fit + n.test
  
  list_K <- Klein_bottle_f(n)
  
  Klein_data <- list_K$Klein_data ## Generation of 300 * 3 matrices
  
  X <- list_K$X 
  
  y_K <- c()
  for (i in 1:n) {
    y_K[i] <- f_k(Klein_data[i,1], Klein_data[i,2], Klein_data[i,3])
  }
  
  scalar <- rnorm(n, 0,sd = 1)
  y_K <- y_K + scalar
  f_K <- y_K
  y_K <- y_K + rnorm(n, 0, sd = sigma) 
  y_K.test <- y_K[-c(1:n.fit)]
  
  
  dist.matrix <- array(0,dim = c(n, n, 6))
  
  dist.matrix[,,1] <- as.matrix(dist(X[,1:100], diag = T, upper = T, method = 'manhattan')) 
  dist.matrix[,,2] <- as.matrix(dist(X[,101:250], diag = T, upper = T, method = 'manhattan'))
  dist.matrix[,,3] <- as.matrix(dist(X[,251:450], diag = T, upper = T, method = 'manhattan'))
  
  dist.matrix[,,4] <- dist.matrix[,,1] + dist.matrix[,,2]
  dist.matrix[,,5] <- dist.matrix[,,1] + dist.matrix[,,3]
  dist.matrix[,,6] <- dist.matrix[,,2] + dist.matrix[,,3]
  
  
  result_rep = foreach(i.rep = 1:200, .combine = cbind) %dopar%
    { 
  ###################### Selection of representative samples ################################
  
  fit_error_m2 <- predict_error_m2 <- 0
  
  for (l in 1:length(scale)) {
    N = n.fit
    nq0 = floor(scale[l] * N)+1
    nk = 20  
    
    Y = y_K[1:n.fit]
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
      loc = loc[sample(1:ninterv[i])]
      id.basis = c(id.basis, loc[1:nnum[i]])
    }
    
    id.basis = sort(unique(id.basis))
    N.gss = length(id.basis)
    
    ######
    dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
    dist.test <- dist.matrix[(n.fit+1):n, , ]
    timeend_1<- timeend_2 <- timeend_3 <- timeend_4 <- timeend_5 <-0
    
    timestart_1 <- Sys.time()
    Q_rep <- array(0, dim = c(N.gss, N.gss, 6))
    for (k in 1:6) {
      if(k==4){timeend_1 = Sys.time()}
      Q_rep[,,k] <- MetricKern(dist.train[id.basis,id.basis,k]/h)
    }
    timeend_4 <- Sys.time()
    
    d_change <- array(0, dim = c(n.fit, N.gss, 6))
    for (i in 1:6) {
      d_change[, , i] <- rbind(dist.train[id.basis, id.basis, i],
                               dist.train[-id.basis, id.basis, i])
    }
    
    timestart_2 <- Sys.time()
    Q_patch <- array(0, dim = c(N-N.gss, N.gss, 6))
    for (k in 1:6) {
      if(k==4){timeend_2 = Sys.time()}
      Q_patch[,,k] <- MetricKtest(d_change[,,k]/h, N.gss, N-N.gss)
    }
    
    timeend1 <- Sys.time()
    
    Kermat_train <- array(0, dim = c(n.fit, n.fit, 6))
    Kermat_train[id.basis, id.basis, ] <- Q_rep
    Kermat_train[-id.basis, id.basis, ] <- Q_patch
    Kermat_train_tensor <- Kermat_train
    
    Kermat_train_tensor[,,4] <- Kermat_train[,,1] * Kermat_train[,,2]
    Kermat_train_tensor[,,5] <- Kermat_train[,,1] * Kermat_train[,,3]
    Kermat_train_tensor[,,6] <- Kermat_train[,,2] * Kermat_train[,,3]
    
    s = matrix(scalar[1:n.fit], n.fit, 1)
    fit_ms_k <- tryCatch(gss:::mspreg1(s = s,r=Kermat_train_tensor[,id.basis,],id.basis=id.basis,method = 'm',
                                       y = y_K[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                         error = function(x) NA)
    if((any(is.na(fit_ms_k))))
    { 
      fit_error_m2[l] = NA
      predict_error_m2[l] =NA
    }else{
      theta_k <- fit_ms_k$theta
      c_k <- matrix(0, n.fit, 1)
      c_k[id.basis,1] <- fit_ms_k$c
      d_k <- fit_ms_k$d
      
      ## fitting error
      hat_k <- matrix(0, n.fit, 1)
      K <- matrix(0, n.fit, n.fit)
      
      for (i in 1:6) {
        K <- 10^theta_k[i] * Kermat_train_tensor[,,i] + K
      }
      
      hat_k <- c(K %*% c_k) + scalar[1:n.fit] * d_k
      fit_error_m2[l] <- mean((hat_k-f_K[1:n.fit])^2)
      
      ## prediction error
      dx.test <- array(0, dim = c(n.test + N.gss, N.gss, 6))
      for (i in 1:6) {
        dx.test[, , i] <- rbind(dist.train[id.basis, id.basis, i],
                                dist.test[, id.basis, i])
      }
      
      kermat <- array(0, dim = c(n.test, N.gss, 6))
      for (j in 1:6) {
        kermat[, , j] <- MetricKtest(dx.test[, , j]/h, N.gss, n.test)
      }
      
      kermat[,,4] <- kermat[,,1] * kermat[,,2]
      kermat[,,5] <- kermat[,,1] * kermat[,,3]
      kermat[,,6] <- kermat[,,2] * kermat[,,3]
      
      Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
      Kermat_test[,id.basis,] <- kermat
      
      p_y_k <- matrix(0, n.test, 1)
      p_K_k <- matrix(0, n.test, n.fit)
      for (i in 1:6) {
        p_K_k <- 10^theta_k[i] * Kermat_test[,,i] + p_K_k
      }
      p_y_k <- c(p_K_k %*% c_k) + scalar[(n.fit+1):n]*d_k
      predict_error_m2[l] <- mean((p_y_k-y_K.test)^2)
    }
  }
  
  
  
  c(predict_error_m2, fit_error_m2) 
}
  
  return(result_rep)
  
}

####################################
set.seed(2022)


sim_m1_200 = simulation_M_rep(n.fit = 200, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.05)

sim_m2_200 = simulation_M_rep(n.fit = 200, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.1)

sim_m3_200 = simulation_M_rep(n.fit = 200, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.15)


sim_m1 = simulation_M_rep(n.fit = 300, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.05)

sim_m2 = simulation_M_rep(n.fit = 300, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.1)

sim_m3 = simulation_M_rep(n.fit = 300, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.15)



sim_m1_400 = simulation_M_rep(n.fit = 400, n.test = 100, h=0.1, 
                              scale=c(1/4), sigma=0.05)

sim_m2_400 = simulation_M_rep(n.fit = 400, n.test = 100, h=0.1, 
                              scale=c(1/4), sigma=0.1)

sim_m3_400 = simulation_M_rep(n.fit = 400, n.test = 100, h=0.1, 
                              scale=c(1/4), sigma=0.15)


data = rbind(sim_m1_200, sim_m2_200, sim_m3_200, sim_m1, sim_m2, sim_m3, sim_m1_400, sim_m2_400, sim_m3_400)
save(data,file="simu_rep_Manifold.Rdata")

##########################################################

set.seed(2022)

sim_k1_200 = simulation_K_rep(n.fit = 200, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.05)

sim_k2_200 = simulation_K_rep(n.fit = 200, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.1)

sim_k3_200 = simulation_K_rep(n.fit = 200, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.15)

sim_k1 = simulation_K_rep(n.fit = 300, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.05)

sim_k2 = simulation_K_rep(n.fit = 300, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.1)

sim_k3 = simulation_K_rep(n.fit = 300, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.15)


sim_k1_400 = simulation_K_rep(n.fit = 400, n.test = 100, h=1, 
                              scale=c(1/4), sigma=0.05)

sim_k2_400 = simulation_K_rep(n.fit = 400, n.test = 100, h=1, 
                              scale=c(1/4), sigma=0.1)

sim_k3_400 = simulation_K_rep(n.fit = 400, n.test = 100, h=1, 
                              scale=c(1/4), sigma=0.15)


data = rbind(sim_k1_200, sim_k2_200, sim_k3_200, sim_k1, sim_k2, sim_k3, sim_k1_400, sim_k2_400, sim_k3_400)
save(data,file="simu_rep_Klein.Rdata")



##################################################################
load("simu_rep_Manifold.Rdata")
res = cbind(apply(data, 1, mean), apply(data, 1, sd))

load("simu_rep_Klein.Rdata")
res2 = cbind(apply(data, 1, mean), apply(data, 1, sd))

res = round(res, 4)
res2 = round(res2, 4)


resu = NULL


i=1
j=c(2,4,6,1,3,5)
tmp1 = res[j ,i]
tmp2 = res[j, i+1]

resu = rbind(resu, paste("200 &SO(3)& MSE &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")\\",sep=""))

resu = rbind(resu, paste("&& PSE &", tmp1[4], "(", tmp2[4],")&",
                         tmp1[5], "(", tmp2[5],")&", tmp1[6], "(", tmp2[6],")\\",sep=""))

tmp1 = res2[j ,i]
tmp2 = res2[j, i+1]

resu = rbind(resu, paste("&Klein & MSE &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")\\",sep=""))

resu = rbind(resu, paste("&& PSE &", tmp1[4], "(", tmp2[4],")&",
                         tmp1[5], "(", tmp2[5],")&", tmp1[6], "(", tmp2[6],")\\",sep=""))


j=c(8,10,12,7,9,11)
tmp1 = res[j ,i]
tmp2 = res[j, i+1]

resu = rbind(resu, paste("300&SO(3)& MSE &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")\\",sep=""))

resu = rbind(resu, paste("&& PSE &", tmp1[4], "(", tmp2[4],")&",
                         tmp1[5], "(", tmp2[5],")&", tmp1[6], "(", tmp2[6],")\\",sep=""))

tmp1 = res2[j ,i]
tmp2 = res2[j, i+1]

resu = rbind(resu, paste("&Klein & MSE &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")\\",sep=""))

resu = rbind(resu, paste("&& PSE &", tmp1[4], "(", tmp2[4],")&",
                         tmp1[5], "(", tmp2[5],")&", tmp1[6], "(", tmp2[6],")\\",sep=""))




j=c(14,16,18,13,15,17)
tmp1 = res[j ,i]
tmp2 = res[j, i+1]

resu = rbind(resu, paste("400&SO(3)& MSE &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")\\",sep=""))

resu = rbind(resu, paste("&& PSE &", tmp1[4], "(", tmp2[4],")&",
                         tmp1[5], "(", tmp2[5],")&", tmp1[6], "(", tmp2[6],")\\",sep=""))

tmp1 = res2[j ,i]
tmp2 = res2[j, i+1]

resu = rbind(resu, paste("&Klein & MSE &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")\\",sep=""))

resu = rbind(resu, paste("&& PSE &", tmp1[4], "(", tmp2[4],")&",
                         tmp1[5], "(", tmp2[5],")&", tmp1[6], "(", tmp2[6],")\\",sep=""))



resu[c(1,3,5,7,9,11) ,]


write.table(resu[c(1,3,5,7,9,11) ,], file="simu-rep-comp.txt", col.names=F, row.names=F)



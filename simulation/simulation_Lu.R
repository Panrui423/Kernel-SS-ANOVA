library(gss)
library("doMC")
library("foreach")
source("Mat_kernel.R")
source("CommonMethods.R")
source("Manifold_simulation.R")
source("Klein_simulation.R")

registerDoMC(20)
###########  #############
# warning! It takes about 3 days to run.
############################


simulationLu_M <- function(n.fit, n.test, sigma){
  
  n <- n.fit + n.test
  
  list_M <- Manifold_f(n)
  
  Manifold_data <- list_M$Manifold_data    ## Generation of 300 * 3 matrices
  
  X <- list_M$X 
  
  y_M <- c()
  for (i in 1:n) {
    y_M[i] <- f_m(Manifold_data[i,1], Manifold_data[i,2], Manifold_data[i,3])
  }
  y_M <- y_M + rnorm(n, 0,sd = sigma) 
  scalar <- rnorm(n, 0,sd = 1)
  y_M <- y_M + scalar
  y_M.train <- y_M[1:n.fit]
  y_M.test <- y_M[-c(1:n.fit)]
  
  dist.matrix <- array(0,dim = c(n, n, 3))
  
  dist.matrix[,,1] <- as.matrix(dist(X[,1:100], diag = T, upper = T, method = 'manhattan'))  
  dist.matrix[,,2] <- as.matrix(dist(X[,101:250], diag = T, upper = T, method = 'manhattan'))
  dist.matrix[,,3] <- as.matrix(dist(X[,251:450], diag = T, upper = T, method = 'manhattan'))
  
  ###############
  
  dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
  dist.test <- dist.matrix[(n.fit+1):n, 1:n.fit, ]
  train.num <- n.fit
  test.num <- n.test
  pvar <- 3
  pall <- 6
  
  start = Sys.time()
  #print(start)
  Kernel_train = array(0, dim = c(train.num, train.num, pvar))
  Kernel_test = array(0, dim = c(test.num, train.num, pvar))
  for(i in 1:pvar)
  {
    dist_train = dist.train[ , , i]
    Kernel_Lu  = gen_kernel_matrix(D=dist_train, sample=TRUE, n_size=45, lambda = 1)
    
    test_mat = foreach(k = 1:test.num, .combine = rbind) %dopar%
      {
        newmat  = dist.test[k, , i]
        mat = Newkernel_solve(d=newmat, Kernel_matrix=Kernel_Lu, n=train.num)
        mat[train.num+1, 1:train.num]
      }
    Kernel_train[,,i] = Kernel_Lu
    Kernel_test[,,i] = test_mat
  }
  
  Kernel_train = Tensor_f(Kernel_train)
  Kernel_test = Tensor_f(Kernel_test)
  
  s_train = matrix(scalar[1:n.fit], n.fit, 1)
  fit_lu <- gss:::mspreg1(s = s_train, r = Kernel_train, id.basis = 1:train.num, method = 'm',
                          y = y_M.train, wt = 1, alpha = 1.4, random = NULL, varht = 1, skip.iter = FALSE)
  
  theta_fit <- fit_lu$theta
  c_fit <- fit_lu$c
  d_fit <- fit_lu$d
  
  ## fitting error
  hat_m <- matrix(0, train.num, 1)
  M.com <- computeSigma(theta_fit, Kernel_train, pall, train.num)
  
  hat_m <- c(M.com %*% c_fit) + s_train %*% d_fit
  fit_error_4 <- mean((hat_m - y_M.train)^2)
  
  ## prediction error
  s_test <- scalar[-(1:n.fit)]
  p_y_m <- matrix(0, test.num, 1)
  p_K_m <- matrix(0, test.num, train.num)
  for (i in 1:pall) {
    p_K_m <- 10^theta_fit[i] * Kernel_test[,,i] + p_K_m
  }
  p_y_m <- c(p_K_m %*% as.matrix(c_fit)) + s_test * d_fit
  predict_error_m4 <- c(mean((p_y_m - y_M.test)^2),  mean(abs((p_y_m - y_M.test)/y_M.test)))
  
  end = Sys.time()
  
  res = c(predict_error_m4, fit_error_4, end-start)
  res
  
}

###################
simulationLu_K <- function(n.fit, n.test, sigma){
  
  n <- n.fit + n.test
  
  list_K <- Klein_bottle_f(n)
  
  Klein_data <- list_K$Klein_data ## Generation of 300 * 3 matrices
  
  X <- list_K$X 
  
  y_K <- c()
  for (i in 1:n) {
    y_K[i] <- f_k(Klein_data[i,1], Klein_data[i,2], Klein_data[i,3])
  }
  y_K <- y_K + rnorm(n, 0, sd = sigma) 
  scalar <- rnorm(n, 0,sd = 1)
  y_K <- y_K + scalar
  y_K.train <- y_K[1:n.fit]
  y_K.test <- y_K[-c(1:n.fit)]
  
  
  dist.matrix <- array(0,dim = c(n, n, 6))
  
  dist.matrix[,,1] <- as.matrix(dist(X[,1:100], diag = T, upper = T, method = 'manhattan')) 
  dist.matrix[,,2] <- as.matrix(dist(X[,101:250], diag = T, upper = T, method = 'manhattan'))
  dist.matrix[,,3] <- as.matrix(dist(X[,251:450], diag = T, upper = T, method = 'manhattan'))
  
  ###############
  
  
  dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
  dist.test <- dist.matrix[(n.fit+1):n, 1:n.fit, ]
  train.num <- n.fit
  test.num <- n.test
  pvar <- 3
  pall <- 6
  
  start = Sys.time()
  
  Kernel_train = array(0, dim = c(train.num, train.num, pvar))
  Kernel_test = array(0, dim = c(test.num, train.num, pvar))
  for(i in 1:pvar)
  {
    dist_train = dist.train[ , , i]
    Kernel_Lu  = gen_kernel_matrix(D=dist_train, sample=TRUE, n_size=45, lambda = 1)
    
    test_mat = foreach(k = 1:test.num, .combine = rbind) %dopar%
      {
        newmat  = dist.test[k, , i]
        mat = Newkernel_solve(d=newmat, Kernel_matrix=Kernel_Lu, n=train.num)
        mat[train.num+1, 1:train.num]
      }
    Kernel_train[,,i] = Kernel_Lu
    Kernel_test[,,i] = test_mat
  }
  
  Kernel_train = Tensor_f(Kernel_train)
  Kernel_test = Tensor_f(Kernel_test)
  
  s_train = matrix(scalar[1:n.fit], n.fit, 1)
  fit_lu <- gss:::mspreg1(s = s_train, r = Kernel_train, id.basis = 1:train.num, method = 'm',
                          y = y_K.train, wt = 1, alpha = 1.4, random = NULL, varht = 1, skip.iter = FALSE)
  
  theta_fit <- fit_lu$theta
  c_fit <- fit_lu$c
  d_fit <- fit_lu$d
  
  ## fitting error
  hat_m <- matrix(0, train.num, 1)
  M.com <- computeSigma(theta_fit, Kernel_train, pall, train.num)
  
  hat_m <- c(M.com %*% c_fit) + s_train %*% d_fit
  fit_error_4 <- mean((hat_m - y_K.train)^2)
  
  ## prediction error
  s_test = scalar[-(1:n.fit)]
  p_y_m <- matrix(0, test.num, 1)
  p_K_m <- matrix(0, test.num, train.num)
  for (i in 1:pall) {
    p_K_m <- 10^theta_fit[i] * Kernel_test[,,i] + p_K_m
  }
  p_y_m <- c(p_K_m %*% as.matrix(c_fit)) + s_test * d_fit
  predict_error_m4 <- c(mean((p_y_m - y_K.test)^2),  mean(abs((p_y_m - y_K.test)/y_K.test)))
  
  end = Sys.time()
  
  res = c(predict_error_m4, fit_error_4, end - start)
  res
  
}


###### Manifold ##############
set.seed(2022)


result_Lu_M = NULL
for(i.rep in 1:200)
{
  resu = simulationLu_M(n.fit = 100, n.test = 100, sigma=0.05)
  result_Lu_M = cbind(result_Lu_M, resu)
  write.table(result_Lu_M, file="result_Lu_M100_05.txt", col.names=F, row.names=F)
  print(c(i.rep, resu[1:4]))
  
}



result = result_Lu_M
round(apply(result, 1, function(x){mean(x)}), 4)

############

set.seed(2022)


result_Lu_M = NULL
for(i.rep in 1:200)
{
  resu = simulationLu_M(n.fit = 100, n.test = 100, sigma=0.10)
  result_Lu_M = cbind(result_Lu_M, resu)
  write.table(result_Lu_M, file="result_Lu_M100.txt", col.names=F, row.names=F)
  print(c(i.rep, resu[1:4]))
  
}



result = result_Lu_M
round(apply(result, 1, function(x){mean(x)}), 4)


#########
set.seed(2022)


result_Lu_M = NULL
for(i.rep in 1:200)
{
  resu = simulationLu_M(n.fit = 100, n.test = 100, sigma=0.15)
  result_Lu_M = cbind(result_Lu_M, resu)
  write.table(result_Lu_M, file="result_Lu_M100_15.txt", col.names=F, row.names=F)
  print(c(i.rep, resu[1:4]))
  
}



result = result_Lu_M
round(apply(result, 1, function(x){mean(x)}), 4)


###### Klein bottle ##############
set.seed(2022)

result_Lu_K = NULL
for(i.rep in 1:200)
  #result_Lu_K = foreach(i.rep = 1:200, .combine = cbind) %dopar%
{
  resu = simulationLu_K(n.fit = 100, n.test = 100, sigma=0.05)
  # resu
  result_Lu_K = cbind(result_Lu_K, resu)
  write.table(result_Lu_K, file="result_Lu_K100_05.txt", col.names=F, row.names=F)
  print(c(i.rep, resu[1:4]))
  
}



result = result_Lu_K
round(apply(result, 1, function(x){mean(x)}), 4)

##################

set.seed(2022)

result_Lu_K = NULL
for(i.rep in 1:200)
#result_Lu_K = foreach(i.rep = 1:200, .combine = cbind) %dopar%
{
  resu = simulationLu_K(n.fit = 100, n.test = 100, sigma=0.10)
 # resu
  result_Lu_K = cbind(result_Lu_K, resu)
  write.table(result_Lu_K, file="result_Lu_K100.txt", col.names=F, row.names=F)
  print(c(i.rep, resu[1:4]))
  
}



result = result_Lu_K
round(apply(result, 1, function(x){mean(x)}), 4)

###########
set.seed(2022)

result_Lu_K = NULL
for(i.rep in 1:200)
  #result_Lu_K = foreach(i.rep = 1:200, .combine = cbind) %dopar%
{
  resu = simulationLu_K(n.fit = 100, n.test = 100, sigma=0.15)
  # resu
  result_Lu_K = cbind(result_Lu_K, resu)
  write.table(result_Lu_K, file="result_Lu_K100_15.txt", col.names=F, row.names=F)
  print(c(i.rep, resu[1:4]))
  
}



result = result_Lu_K
round(apply(result, 1, function(x){mean(x)}), 4)
##############################################################
result_Lu_M = read.table("result_Lu_M100_05.txt")
result_Lu_K = read.table("result_Lu_K100_05.txt")

result1 = result_Lu_M
result2 = result_Lu_K
est05 = cbind(round(apply(result1, 1, function(x){mean(x)}), 4), round(apply(result1, 1, function(x){sd(x)}), 4),
              round(apply(result2, 1, function(x){mean(x)}), 4), round(apply(result2, 1, function(x){sd(x)}), 4) )



result_Lu_M = read.table("result_Lu_M100.txt")
result_Lu_K = read.table("result_Lu_K100.txt")

result1 = result_Lu_M
result2 = result_Lu_K
est10 = cbind(round(apply(result1, 1, function(x){mean(x)}), 4), round(apply(result1, 1, function(x){sd(x)}), 4),
              round(apply(result2, 1, function(x){mean(x)}), 4), round(apply(result2, 1, function(x){sd(x)}), 4) )


result_Lu_M = read.table("result_Lu_M100_15.txt")
result_Lu_K = read.table("result_Lu_K100_15.txt")

result1 = result_Lu_M
result2 = result_Lu_K
est15 = cbind(round(apply(result1, 1, function(x){mean(x)}), 4), round(apply(result1, 1, function(x){sd(x)}), 4),
              round(apply(result2, 1, function(x){mean(x)}), 4), round(apply(result2, 1, function(x){sd(x)}), 4) )

################################  table  #####################
load("simu_runtime100.Rdata")

result = result100
est100 = matrix(round(apply(result, 1, function(x){mean(x,na.rm = TRUE)}), 4),12,)
sd100 = matrix(round(apply(result, 1, function(x){sd(x,na.rm = TRUE)}), 4),12,)

est100 = est100[c(1, 3, 5, 7), ]
dimnames(est100)[[2]] <- c("SO(3)", "SO(3)","SO(3)", "Klein", "Klein", "Klein")
dimnames(est100)[[1]] <- c("Our(1/4)", "Our(all)", "Gaussian", "Laplacian")


sd100 = sd100[c(1, 3, 5, 7), ]
dimnames(sd100)[[2]] <- c("SO(3)", "SO(3)","SO(3)", "Klein", "Klein", "Klein")
dimnames(sd100)[[1]] <- c("Our(1/4)", "Our(all)", "Gaussian", "Laplacian")

est100





resu = NULL

tmp1 = est100[ , 1]
tmp2 = sd100[ , 1]

resu = rbind(resu, paste(0.05, "& SO(3) &", 100, "&", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", est05[1,1],"(",est05[1,2],")\\",sep=""))

tmp1 = est100[ , 4]
tmp2 = sd100[ , 4]

resu = rbind(resu, paste(  "   & Klein &",  "   &", tmp1[1], "(", tmp2[1],")&",
                           tmp1[2], "(", tmp2[2],")&",est05[1,3],"(",est05[1,4],")\\",sep=""))


tmp1 = est100[ , 2]
tmp2 = sd100[ , 2]

resu = rbind(resu, paste(0.10, "& SO(3) &", 100, "&", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", est10[1,1],"(",est10[1,2],")\\",sep=""))

tmp1 = est100[ , 5]
tmp2 = sd100[ , 5]

resu = rbind(resu, paste(  "   & Klein &",  "   &", tmp1[1], "(", tmp2[1],")&",
                           tmp1[2], "(", tmp2[2],")&", est10[1,3],"(",est10[1,4],")\\",sep=""))


tmp1 = est100[ , 3]
tmp2 = sd100[ , 3]

resu = rbind(resu, paste(0.15, "& SO(3) &", 100, "&", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", est15[1,1],"(",est15[1,2],")\\",sep=""))

tmp1 = est100[ , 6]
tmp2 = sd100[ , 6]

resu = rbind(resu, paste(  "   & Klein &",  "   &", tmp1[1], "(", tmp2[1],")&",
                           tmp1[2], "(", tmp2[2],")&", est15[1,3],"(",est15[1,4],")\\",sep=""))


resu

write.table(resu, file = "Lucomp.txt", col.names = F, row.names = F)


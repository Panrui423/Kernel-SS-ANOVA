library(gss)

Rcpp::sourceCpp('Kernelmatrix.cpp')

set.seed(2021)
#########  Klein bottle ###########

c0_k = 0.0328

a_l <- function(l, t){
  if(l%%2 != 0){
    return(cos(l*pi*t/10) / 5^(1/2))
  } 
  else{
    return(sin((l-1)*pi*t/10) / 5^(1/2))
  }
}

count_k <- function(){
  u_i1_k <- runif(3, 0, 2*pi)
  u_i2_k <- runif(3, 0, 2*pi)
  
  X <- c()
  X_il_w <- c()
  
  for (i in 1:3) {
    nu_i1_k <- (2*cos(u_i2_k[i])+1) * cos(u_i1_k[i])
    nu_i2_k <- (2*cos(u_i2_k[i])+1) * sin(u_i1_k[i])
    nu_i3_k <- (2*sin(u_i2_k[i])) * cos(u_i1_k[i]/2)
    nu_i4_k <- (2*sin(u_i2_k[i])) * sin(u_i1_k[i]/2)
    
    X_1i <- function(t){
      a <- (nu_i1_k * a_l(1,t) + nu_i2_k * a_l(2,t) + nu_i3_k * a_l(3,t) + nu_i4_k * a_l(4,t)) / c0_k
      a  
    } 
    
    t <- seq(0,1,length = (50 + i*50))
    X <- append(X, X_1i(t))
    
    inte_func_k <- function(t){ 
      X_1i(t)^2 * t/500
    }
    
    X_il_w[i] <- integrate(inte_func_k, 0, 1)$value
  }
  
  X <- matrix(X, 1, 450)
  X_il_w <- matrix(X_il_w, 1, 3)
  
  return(list(X = X, X_il_w = X_il_w))
}


Klein_bottle_f <- function(n){
  Klein_data <- matrix(0, n, 3) 
  X <- matrix(0, n, 100+150+200)
  for (i in 1:n) {
    a <- count_k()
    Klein_data[i,] <- a$X_il_w
    X[i,] <- a$X
  }
  return(list(Klein_data = Klein_data, X = X))
}

###### f_1 f_2 f_12 f_13

f_1 <- function(x){
  2 * sin(5*pi*x) * cos(2*pi*x^2) + 2*gamma(1 + x/2)
}

f_2 <- function(x){
  2 * sin(5*pi*x) * cos(2*pi*x^2) + 2*gamma(1 + x/2)
}

f_12 <- function(x1, x2){
  2 * sin(5*pi*x1) * cos(2*pi*x2^2) + 2*gamma(1 + x1/4 + x2/4)
}

f_13 <- function(x1, x3){
  2 * sin(5*pi*x1) * cos(2*pi*x3^2) + 2*gamma(1 + x1/4 + x3/4)
}

f_k <- function(x1, x2, x3){
  (6*f_1(x1) + f_2(x2) + f_12(x1, x2) + f_13(x1, x3))/10
} 

kernalGaussian <- function(s,t,para){
  h <- length(s)*para
  kernalX <-exp(-sum((s-t)^2)/(2*h^2))
  return(kernalX)
}

kernalLap_r <- function(X, para)
{
  n <- dim(X)[1]
  p <- dim(X)[2]
  Lap.matrix = matrix(0, n, n)
  h <- p*para
  for (i in 1:n) {
    for (j in 1:i) {
      Lap.matrix[i,j] = exp(-sum((X[i,] - X[j,])^2)^0.5/(2*h)) 
      Lap.matrix[j,i] = Lap.matrix[i,j]
    }
  }
  return(Lap.matrix)
}

simulation_K <- function(n.fit, n.test, h, scale, sigma){
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
  y_K.test <- y_K[-c(1:n.fit)]
  
  
  dist.matrix <- array(0,dim = c(n, n, 6))
  
  dist.matrix[,,1] <- as.matrix(dist(X[,1:100], diag = T, upper = T, method = 'manhattan')) 
  dist.matrix[,,2] <- as.matrix(dist(X[,101:250], diag = T, upper = T, method = 'manhattan'))
  dist.matrix[,,3] <- as.matrix(dist(X[,251:450], diag = T, upper = T, method = 'manhattan'))
  
  dist.matrix[,,4] <- dist.matrix[,,1] + dist.matrix[,,2]
  dist.matrix[,,5] <- dist.matrix[,,1] + dist.matrix[,,3]
  dist.matrix[,,6] <- dist.matrix[,,2] + dist.matrix[,,3]
  
  
  
  ###################### Selection of representative samples ################################
  
  fit_error_m1 <- predict_error_m1 <- fit_error_m2 <- predict_error_m2 <- 
    runningtime1 <- runningtime2 <- c()
  
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
    
    ######  Method1
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
    runningtime1[l] <- (timeend1-timestart_2) + (timeend_4-timestart_1)
    
    Kermat_train <- array(0, dim = c(n.fit, n.fit, 6))
    Kermat_train[id.basis, id.basis, ] <- Q_rep
    Kermat_train[-id.basis, id.basis, ] <- Q_patch
    
    s = matrix(scalar[1:n.fit], n.fit, 1)
    fit_ms_k <- tryCatch(gss:::mspreg1(s = s,r=Kermat_train[,id.basis,],id.basis=id.basis,method = 'm',
                                       y = y_K[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                         error = function(x) NA)
    
    if((any(is.na(fit_ms_k))))
    { 
      fit_error_m1[l] = NA
      predict_error_m1[l] =NA
    }else{
      theta_k <- fit_ms_k$theta
      c_k <- matrix(0, n.fit, 1)
      c_k[id.basis,1] <- fit_ms_k$c
      d_k <- fit_ms_k$d
      
      ## fitting error
      hat_k <- matrix(0, n.fit, 1)
      K <- matrix(0, n.fit, n.fit)
      
      for (i in 1:6) {
        K <- 10^theta_k[i] * Kermat_train[,,i] + K
      }
      
      hat_k <- c(K %*% c_k) + scalar[1:n.fit] * d_k
      fit_error_m1[l] <- mean((hat_k-y_K[1:n.fit])^2)
      
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
      
      Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
      Kermat_test[,id.basis,] <- kermat
      
      p_y_k <- matrix(0, n.test, 1)
      p_K_k <- matrix(0, n.test, n.fit)
      for (i in 1:6) {
        p_K_k <- 10^theta_k[i] * Kermat_test[,,i] + p_K_k
      }
      p_y_k <- c(p_K_k %*% c_k) + scalar[(n.fit+1):n]*d_k
      predict_error_m1[l] <- mean((p_y_k-y_K.test)^2)
    }
    
    ######  Method2 #######
    
    timestart2 <- Sys.time()
    Kermat_train_tensor <- Kermat_train
    
    Kermat_train_tensor[,,4] <- Kermat_train[,,1] * Kermat_train[,,2]
    Kermat_train_tensor[,,5] <- Kermat_train[,,1] * Kermat_train[,,3]
    Kermat_train_tensor[,,6] <- Kermat_train[,,2] * Kermat_train[,,3]
    
    timeend2 <- Sys.time()
    runningtime2[l] <-(timeend_1-timestart_1) + 
      (timeend_2-timestart_2) + 
      (timeend2 - timestart2)
    
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
      fit_error_m2[l] <- mean((hat_k-y_K[1:n.fit])^2)
      
      ## prediction error
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
  
  ######################### No representative sample selected #####################
  
  #### Method1 ######
  dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
  dist.test <- dist.matrix[(n.fit+1):n, , ]
  timeend_1<- timeend_2 <- timeend_3 <- 0
  
  timestart_1 <- Sys.time()
  Q_K <- array(0, dim = c(n.fit, n.fit, 6))
  for (k in 1:6) {
    if(k==4){timeend_1 = Sys.time()}
    Q_K[,,k] <- MetricKern(dist.matrix[1:n.fit, 1:n.fit, k]/h)
  }
  
  time2 <- Sys.time()
  runningtime3 <- time2-timestart_1
  
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_k <- tryCatch(gss:::mspreg1(s = s, r=Q_K[,1:n.fit,],id.basis=1:n.fit,method = 'm',
                            y = y_K[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                       error = function(x) NA)
  if((any(is.na(fit_ms_k))))
  { 
    fit_error_m1_no = NA
    predict_error_m1_no =NA
  }else{
  theta_k <- fit_ms_k$theta
  c_k <- fit_ms_k$c
  d_k <- fit_ms_k$d
  
  ## fitting error
  hat_k <- matrix(0, n.fit, 1)
  K <- matrix(0, n.fit, n.fit)
  
  for (i in 1:6) {
    K <- 10^theta_k[i] * Q_K[,,i] + K
  }
  
  hat_k <- c(K %*% matrix(c_k,n.fit,1)) + scalar[1:n.fit]*d_k
  fit_error_m1_no <- mean((hat_k-y_K[1:n.fit])^2)
  
  ## prediction error
  dx.test <- array(0, dim = c(n.test + n.fit, n.fit, 6))
  for (i in 1:6) {
    dx.test[, , i] <- rbind(dist.train[1:n.fit, 1:n.fit, i],
                            dist.test[, 1:n.fit, i])
  }
  
  kermat <- array(0, dim = c(n.test, n.fit, 6))
  for (j in 1:6) {
    kermat[, , j] <- MetricKtest(dx.test[, , j]/h, n.fit, n.test)
  }
  
  Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
  Kermat_test[,1:n.fit,] <- kermat
  
  p_y_k <- matrix(0, n.test, 1)
  p_K_k <- matrix(0, n.test, n.fit)
  for (i in 1:6) {
    p_K_k <- 10^theta_k[i] * Kermat_test[,,i] + p_K_k
  }
  
  p_y_k <- c(p_K_k %*% matrix(c_k, n.fit, 1)) + scalar[(n.fit+1):n] * d_k
  
  predict_error_m1_no <- mean((p_y_k-y_K.test)^2)
  }
  
  #### Method2 ######
  
  timestart4 <- Sys.time()
  Q_K[,,4] <- Q_K[,,1] * Q_K[,,2]
  Q_K[,,5] <- Q_K[,,1] * Q_K[,,3]
  Q_K[,,6] <- Q_K[,,2] * Q_K[,,3]
  
  timeend4 <- Sys.time()
  
  runningtime4 <- (timeend_1 - timestart_1) + (timeend4 - timestart4)
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_k <- tryCatch(gss:::mspreg1(s = s, r=Q_K[,1:n.fit,],id.basis=1:n.fit,method = 'm',
                            y = y_K[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                            error=function(x) NA)
  if((any(is.na(fit_ms_k))))
  { 
    fit_error_m2_no = NA
    predict_error_m2_no =NA
  }else{
  theta_k <- fit_ms_k$theta
  c_k <- fit_ms_k$c
  d_k <- fit_ms_k$d
  
  ## fitting error
  hat_k <- matrix(0, n.fit, 1)
  K <- matrix(0, n.fit, n.fit)
  
  for (i in 1:6) {
    K <- 10^theta_k[i] * Q_K[,,i] + K
  }
  
  hat_k <- c(K %*% matrix(c_k,n.fit,1)) + scalar[1:n.fit]*d_k
  fit_error_m2_no <- mean((hat_k-y_K[1:n.fit])^2)
  
  ##  prediction error
  kermat[,,4] <- kermat[,,1] * kermat[,,2]
  kermat[,,5] <- kermat[,,1] * kermat[,,3]
  kermat[,,6] <- kermat[,,2] * kermat[,,3]
  
  Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
  Kermat_test[,1:n.fit,] <- kermat
  
  p_y_k <- matrix(0, n.test, 1)
  p_K_k <- matrix(0, n.test, n.fit)
  for (i in 1:6) {
    p_K_k <- 10^theta_k[i] * Kermat_test[,,i] + p_K_k
  }
  
  p_y_k <- c(p_K_k %*% matrix(c_k, n.fit, 1)) + scalar[(n.fit+1):n] * d_k
  
  predict_error_m2_no <- mean((p_y_k-y_K.test)^2)
  }
  #### Method3 ######
  Q <- array(0,dim = c(n,n,6))
  for (i in 1:n) {
    for (j in 1:n) {
      Q[i,j,1] = kernalGaussian(X[i,1:100], X[j,1:100], para = 5)
    }
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      Q[i,j,2] = kernalGaussian(X[i,101:250], X[j,101:250], para = 5)
    }
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      Q[i,j,3] = kernalGaussian(X[i,251:450], X[j,251:450], para = 5)
    }
  }
  
  Q[,,4] <- Q[,,1] * Q[,,2]
  Q[,,5] <- Q[,,1] * Q[,,3]
  Q[,,6] <- Q[,,2] * Q[,,3]
  
  Q_K <- Q[1:n.fit, 1:n.fit, ]
  
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_k <- tryCatch(gss:::mspreg1(s = s, r=Q_K, id.basis=1:n.fit, method = 'm',
                            y = y_K[1:n.fit], wt=1, alpha = 1.4, random = NULL,varht=1,skip.iter=FALSE),
                            error=function(x) NA)
  if((any(is.na(fit_ms_k))))
  { 
    fit_error_m3 = NA
    predict_error_m3 =NA
  }else{
    theta_k <- fit_ms_k$theta
    c_k <- fit_ms_k$c
    d_k <- fit_ms_k$d
    
    ## fitting error
    hat_k <- matrix(0, n.fit, 1)
    K <- matrix(0, n.fit, n.fit)
    
    for (i in 1:6) {
      K <- 10^theta_k[i] * Q_K[,,i] + K
    }
    
    hat_k <- c(K %*% matrix(c_k,n.fit,1)) + scalar[1:n.fit]*d_k
    fit_error_m3 <- mean((hat_k-y_K[1:n.fit])^2)
    
    ##  prediction error
    
    Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
    Kermat_test <- Q[(n.fit+1):n, 1:n.fit,]
    
    p_y_k <- matrix(0, n.test, 1)
    p_K_k <- matrix(0, n.test, n.fit)
    for (i in 1:6) {
      p_K_k <- 10^theta_k[i] * Kermat_test[,,i] + p_K_k
    }
    
    p_y_k <- c(p_K_k %*% matrix(c_k, n.fit, 1)) + scalar[(n.fit+1):n] * d_k
    
    predict_error_m3 <- mean((p_y_k - y_K.test)^2)
    
  }
  
  #### Method4 ######
  Q <- array(0,dim = c(n,n,6))
  
  b = 2e-02
  
  Q[,,1] = kernalLap_r(X[,1:100], b)
  Q[,,2] = kernalLap_r(X[,101:250], b)
  Q[,,3] = kernalLap_r(X[,251:450],b)

  Q[,,4] <- Q[,,1] * Q[,,2]
  Q[,,5] <- Q[,,1] * Q[,,3]
  Q[,,6] <- Q[,,2] * Q[,,3]
  
  Q_K <- Q[1:n.fit, 1:n.fit, ]
  
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_k <- tryCatch(gss:::mspreg1(s = s, r=Q_K,id.basis=1:n.fit,method = 'm',
                            y = y_K[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                            error=function(x) NA)
  if((any(is.na(fit_ms_k))))
  { 
    fit_error_m4 = NA
    predict_error_m4 = NA
  }else{
    theta_k <- fit_ms_k$theta
    c_k <- fit_ms_k$c
    d_k <- fit_ms_k$d
    
    ## fitting error
    hat_k <- matrix(0, n.fit, 1)
    K <- matrix(0, n.fit, n.fit)
    
    for (i in 1:6) {
      K <- 10^theta_k[i] * Q_K[,,i] + K
    }
    
    hat_k <- c(K %*% matrix(c_k,n.fit,1)) + scalar[1:n.fit]*d_k
    fit_error_m4 <- mean((hat_k - y_K[1:n.fit])^2)
    
    ## prediction error
    
    Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
    Kermat_test <- Q[(n.fit+1):n, 1:n.fit,]
    
    p_y_k <- matrix(0, n.test, 1)
    p_K_k <- matrix(0, n.test, n.fit)
    for (i in 1:6) {
      p_K_k <- 10^theta_k[i] * Kermat_test[,,i] + p_K_k
    }
    
    p_y_k <- c(p_K_k %*% matrix(c_k, n.fit, 1)) + scalar[(n.fit+1):n] * d_k
    
    predict_error_m4 <- mean((p_y_k - y_K.test)^2)
    
  }
  
  mat <- matrix(c(runningtime1/60, runningtime2/60, 
                  predict_error_m1, fit_error_m1, predict_error_m2, fit_error_m2, 
                  as.double(runningtime3, units = 'mins'), as.double(runningtime4, units = 'mins') , 
                  predict_error_m1_no, fit_error_m1_no, predict_error_m2_no, fit_error_m2_no,
                  fit_error_m3, predict_error_m3,
                  fit_error_m4, predict_error_m4), (10+length(scale)*6), 1)
  
  return(mat)

}


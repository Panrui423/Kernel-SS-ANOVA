library(gss)

Rcpp::sourceCpp('Kernelmatrix.cpp')

set.seed(2021)
################  manifold data ############
c0 = 0.0865

a_l <- function(l, t){
  if(l%%2 != 0){
    return(cos(l*pi*t/10) / 5^(1/2))
  } 
  else{
    return(sin((l-1)*pi*t/10) / 5^(1/2))
  }
}

R <- function(r, theta){
  c <- c(cos(theta), r[3]*sin(theta), -r[2]*sin(theta), 
         -r[3]*sin(theta), cos(theta), r[1]*sin(theta), 
         r[2]*sin(theta), -r[1]*sin(theta), cos(theta))
  (1-cos(theta)) * (r %*% t(r)) + matrix(c, 3, 3)
}

e_2 <- c(0, 1, 0)
e_3 <- c(0, 0, 1)

count_m <- function(){
  u_i1 <- runif(3, 0, 2*pi)
  u_i2 <- runif(3, 0, 2*pi)
  u_i3 <- runif(3, 0, 2*pi)
  
  X <- c()
  X_il_w <- c()
  for (i in 1:3) {
    v <- c(R(e_3,u_i1[i]) * R(e_2,u_i2[i]) * R(e_3,u_i3[i]))
    
    X_1i <-function(t){
      (v[1] * a_l(1,t) + v[2] * a_l(2,t) + v[3] * a_l(3,t) + 
         v[4] * a_l(4,t) + v[5] * a_l(5,t) + v[6] * a_l(6,t) + 
         v[7] * a_l(7,t) + v[8] * a_l(8,t) + v[9] * a_l(9,t))/c0 
    }
    t <- seq(0,1,length = (50 + i*50))
    X <- append(X, X_1i(t))
    
    inte_func_m <-function(t){
      X_1i(t)^2 * t/50
    }
    
    X_il_w[i] <- integrate(inte_func_m, 0, 1)$value
  }
  
  X <- matrix(X, 1, 450)
  X_il_w <- matrix(X_il_w, 1, 3)
  return(list(X = X, X_il_w = X_il_w))
}


Manifold_f <- function(n){
  Manifold_data <- matrix(0, n, 3)
  X <- matrix(0, n, 100+150+200)
  for (i in 1:n) {
    a <- count_m()
    Manifold_data[i,] <- a$X_il_w
    X[i,] <- a$X
  }
  return(list(Manifold_data = Manifold_data, X = X))
}

##### f_1 f_2 f_12 f_13

f_1 <- function(x){
  2 * sin(pi*5*x) * cos(pi*2*x^2) + 2*gamma(1 + x/2)
}

f_2 <- function(x){
  2 * sin(pi*5*x) * cos(pi*2*x^2) + 2*gamma(1 + x/2)
}

f_12 <- function(x1, x2){
  2 * sin(pi*5*x1) * cos(pi*2*x2^2) + 2*gamma(1 + x1/4 + x2/4)
}

f_13 <- function(x1, x3){
  2 * sin(pi*5*x1) * cos(pi*2*x3^2) + 2*gamma(1 + x1/4 + x3/4)
}


f_m <- function(x1, x2, x3){
  ( f_1(x1) * f_2(x2) + f_12(x1, x2) + f_13(x1, x3) )/5
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

### Manifold data
simulation_M <- function(n.fit, n.test, h, scale, sigma){
  
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
  y_M.test <- y_M[-c(1:n.fit)]
  
  dist.matrix <- array(0,dim = c(n, n, 6))
  
  dist.matrix[,,1] <- as.matrix(dist(X[,1:100], diag = T, upper = T, method = 'manhattan'))  
  dist.matrix[,,2] <- as.matrix(dist(X[,101:250], diag = T, upper = T, method = 'manhattan'))
  dist.matrix[,,3] <- as.matrix(dist(X[,251:450], diag = T, upper = T, method = 'manhattan'))
  
  dist.matrix[,,4] <- dist.matrix[,,1] + dist.matrix[,,2]
  dist.matrix[,,5] <- dist.matrix[,,1] + dist.matrix[,,3]
  dist.matrix[,,6] <- dist.matrix[,,2] + dist.matrix[,,3]
  
  
  ###################### Selection of representative samples ################################
  
 fit_error_m2 <- predict_error_m2 <- runningtime1 <- runningtime2 <- runningtime3 <- runningtime4  <- c()
  
  for (l in 1:length(scale)) {
    timestart_M1 <- Sys.time()
    
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
    
    dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
    dist.test <- dist.matrix[(n.fit+1):n, , ]
    
    Q_rep <- array(0, dim = c(N.gss, N.gss, 6))
    for (k in 1:6) {
      Q_rep[,,k] <- MetricKern(dist.train[id.basis,id.basis,k]/h)
    }
    
    d_change <- array(0, dim = c(n.fit, N.gss, 6))
    for (i in 1:6) {
      d_change[, , i] <- rbind(dist.train[id.basis, id.basis, i],
                               dist.train[-id.basis, id.basis, i])
    }
    
    Q_patch <- array(0, dim = c(N-N.gss, N.gss, 6))
    for (k in 1:6) {
      Q_patch[,,k] <- MetricKtest(d_change[,,k]/h, N.gss, N-N.gss)
    }
    
    Kermat_train <- array(0, dim = c(n.fit, n.fit, 6))
    Kermat_train[id.basis, id.basis, ] <- Q_rep
    Kermat_train[-id.basis, id.basis, ] <- Q_patch
    
    ######  Method2 #######
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
      fit_error_m2[l] <- mean((hat_m-y_M[1:n.fit])^2)
      
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
    timeend_M1 <- Sys.time()
    
    runningtime1[l] = timeend_M1 - timestart_M1
  }
  
  ####### No representative sample selected ########
  timestart_M2 <- Sys.time()
  predict_error_m2_no <- 0
  fit_error_m2_no <- 0
  if(TRUE)
  {
  dist.train <- dist.matrix[1:n.fit, 1:n.fit, ]
  dist.test <- dist.matrix[(n.fit+1):n, , ]
  
  Q_M <- array(0, dim = c(n.fit, n.fit, 6))
  for (k in 1:6) {
    Q_M[,,k] <- MetricKern(dist.matrix[1:n.fit, 1:n.fit, k]/h)
  }
  
  #### Method2 ######
  
  Q_M[,,4] <- Q_M[,,1] * Q_M[,,2]
  Q_M[,,5] <- Q_M[,,1] * Q_M[,,3]
  Q_M[,,6] <- Q_M[,,2] * Q_M[,,3]
  
  
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_m <- tryCatch(gss:::mspreg1(s = s, r=Q_M,id.basis=1:n.fit,method = 'm',
                                     y = y_M[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                       error=function(x) NA)
  if((any(is.na(fit_ms_m))))
  { 
    fit_error_m2_no = NA
    predict_error_m2_no =NA
  }else{
    theta_m <- fit_ms_m$theta
    c_m <- fit_ms_m$c
    d_m <- fit_ms_m$d
    
    ## fitting error
    hat_m <- matrix(0, n.fit, 1)
    M <- matrix(0, n.fit, n.fit)
    
    for (i in 1:6) {
      M <- 10^theta_m[i] * Q_M[,,i] + M
    }
    
    hat_m <- c(M %*% matrix(c_m,n.fit,1)) + scalar[1:n.fit]*d_m
    fit_error_m2_no <- mean((hat_m-y_M[1:n.fit])^2)
    
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
    
    kermat[,,4] <- kermat[,,1] * kermat[,,2]
    kermat[,,5] <- kermat[,,1] * kermat[,,3]
    kermat[,,6] <- kermat[,,2] * kermat[,,3]
    
    Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
    Kermat_test[,1:n.fit,] <- kermat
    
    p_y_m <- matrix(0, n.test, 1)
    p_K_m <- matrix(0, n.test, n.fit)
    for (i in 1:6) {
      p_K_m <- 10^theta_m[i] * Kermat_test[,,i] + p_K_m
    }
    
    p_y_m <- c(p_K_m %*% matrix(c_m, n.fit, 1)) + scalar[(n.fit+1):n] * d_m
    
    predict_error_m2_no <- mean((p_y_m-y_M.test)^2)
  }
  }
  timeend_M2 <- Sys.time()
  runningtime2[1] = timeend_M2 - timestart_M2
  
  
  #### Method3 ######
  timestart_M3 <- Sys.time()
  Q <- array(0,dim = c(n,n,6))
  for (i in 1:n) {
    for (j in 1:n) {
      Q[i,j,1] = kernalGaussian(X[i,1:100], X[j,1:100], para = 0.03)  
    }
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      Q[i,j,2] = kernalGaussian(X[i,101:250], X[j,101:250], para = 0.03)  
    }
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      Q[i,j,3] = kernalGaussian(X[i,251:450], X[j,251:450], para = 0.03)
    }
  }
  
  Q[,,4] <- Q[,,1] * Q[,,2]
  Q[,,5] <- Q[,,1] * Q[,,3]
  Q[,,6] <- Q[,,2] * Q[,,3]
  
  Q_M <- Q[1:n.fit, 1:n.fit, ]
  
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_m <- tryCatch(gss:::mspreg1(s = s, r=Q_M, id.basis=1:n.fit, method = 'm',
                                     y = y_M[1:n.fit], wt=1, alpha = 1.4, random = NULL,varht=1,skip.iter=FALSE), 
                       error=function(x) NA)
  if((any(is.na(fit_ms_m))))
  { 
    fit_error_m3 = NA
    predict_error_m3 = NA
  }else{
    theta_m <- fit_ms_m$theta
    c_m <- fit_ms_m$c
    d_m <- fit_ms_m$d
    
    ## fitting error
    hat_m <- matrix(0, n.fit, 1)
    M <- matrix(0, n.fit, n.fit)
    
    for (i in 1:6) {
      M <- 10^theta_m[i] * Q_M[,,i] + M
    }
    
    hat_m <- c(M %*% matrix(c_m,n.fit,1)) + scalar[1:n.fit]*d_m
    fit_error_m3 <- mean((hat_m-y_M[1:n.fit])^2)
    
    ## prediction error
    
    Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
    Kermat_test <- Q[(n.fit+1):n, 1:n.fit,]
    
    p_y_m <- matrix(0, n.test, 1)
    p_K_m <- matrix(0, n.test, n.fit)
    for (i in 1:6) {
      p_K_m <- 10^theta_m[i] * Kermat_test[,,i] + p_K_m
    }
    
    p_y_m <- c(p_K_m %*% matrix(c_m, n.fit, 1)) + scalar[(n.fit+1):n] * d_m
    
    predict_error_m3 <- mean((p_y_m-y_M.test)^2)
    
  }
  
  timeend_M3 <- Sys.time()
  runningtime3[1] = timeend_M3 - timestart_M3
  
  
  #### Method4 ######
  timestart_M4 <- Sys.time()
  
  Q <- array(0,dim = c(n,n,6))
  
  b = 0.02
  
  Q[,,1] = kernalLap_r(X[,1:100], b)
  Q[,,2] = kernalLap_r(X[,101:250], b)
  Q[,,3] = kernalLap_r(X[,251:450], b)
  
  Q[,,4] <- Q[,,1] * Q[,,2]
  Q[,,5] <- Q[,,1] * Q[,,3]
  Q[,,6] <- Q[,,2] * Q[,,3]
  
  Q_M <- Q[1:n.fit, 1:n.fit, ]
  
  s = matrix(scalar[1:n.fit], n.fit, 1)
  fit_ms_m <- tryCatch(gss:::mspreg1(s = s, r=Q_M,id.basis=1:n.fit,method = 'm',
                                     y = y_M[1:n.fit],wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE),
                       error=function(x) NA)
  if((any(is.na(fit_ms_m))))
  { 
    fit_error_m4 = NA
    predict_error_m4 = NA
  }else{
    theta_m <- fit_ms_m$theta
    c_m <- fit_ms_m$c
    d_m <- fit_ms_m$d
    
    ## fitting error
    hat_m <- matrix(0, n.fit, 1)
    M <- matrix(0, n.fit, n.fit)
    
    for (i in 1:6) {
      M <- 10^theta_m[i] * Q_M[,,i] + M
    }
    
    hat_m <- c(M %*% matrix(c_m,n.fit,1)) + scalar[1:n.fit]*d_m
    fit_error_m4 <- mean((hat_m-y_M[1:n.fit])^2)
    
    ## prediction error
    
    Kermat_test <- array(0, dim = c(n.test, n.fit, 6))
    Kermat_test <- Q[(n.fit+1):n, 1:n.fit,]
    
    p_y_m <- matrix(0, n.test, 1)
    p_K_m <- matrix(0, n.test, n.fit)
    for (i in 1:6) {
      p_K_m <- 10^theta_m[i] * Kermat_test[,,i] + p_K_m
    }
    
    p_y_m <- c(p_K_m %*% matrix(c_m, n.fit, 1)) + scalar[(n.fit+1):n] * d_m
    
    predict_error_m4 <- mean((p_y_m-y_M.test)^2)
    
  }
  
  timeend_M4 <- Sys.time()
  runningtime4[1] = timeend_M4 - timestart_M4
  
  ## 
  mat <- matrix(c(predict_error_m2, fit_error_m2, predict_error_m2_no, fit_error_m2_no,  predict_error_m3, fit_error_m3,
                  predict_error_m4, fit_error_m4, as.double(runningtime1)/60, as.double(runningtime2, units = 'mins'),
                  as.double(runningtime3)/60, as.double(runningtime4)/60)) # In mins
  
  return(mat)
}








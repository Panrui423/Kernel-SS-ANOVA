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
source("CommonMethods.R")
source("LoadData.R")

######################################

# [parasubiculum_203] <- 1
# [presubiculum_204] <- 2
# [subiculum_205] <- 3
# [CA1_206] <- 4
# [CA3_208] <- 5
# [CA4_209] <- 6
# [GC_DG_210] <- 7
# [molecular_layer_HP_214] <- 8
# [HATA_211] <- 9
# [fimbria_212] <- 10
# [HP_tail_226] <- 11
# [hippocampal_fissure_215] <- 12


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

Ker4 <- Ker_est_rep(N, dist_train, N.gss, id.basis, pvar, h)


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
h_K_m <- computeSigma(theta_m4, Q_m4, pall, N)

d_m4 <- as.matrix(d_m4) 
round(d_m4, 4)
hat_f_m4 <- h_K_m %*% c_m4
hat_f_m_p <- h_K_m %*% c_m4 + s %*% d_m4  # hat_f_m_p <- h_K_m %*% c_m4 + s %*% d_m4

## MSE
round(mean(((hat_f_m_p - Y.all))^2), 4)


#### Test whether the scalar variable is significant ####

nlambda <- 10^fit_m4$nlambda
M <- h_K_m + nlambda * diag(N)
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



Lap.gram <- array(0, dim = c(m, m, pall))

kernalLap <- function(s,t,para)
{
  h <- length(s)*para
  kernalX <-exp(-sum((s-t)^2)^0.5/(2*h))
  return(kernalX)
}

for (k in 1:c_l) {
  for (i in 1:m) {
    for (j in 1:m) {
      Lap.gram[i,j, k] = kernalLap(Xfun[i, 15000 + which(label_right[,2] == c_reg[k])], 
                                   Xfun[j, 15000 + which(label_right[,2] == c_reg[k])],para = 4e-02)
    }
  }
}

Lap.gram <- Tensor_f(Lap.gram[, , 1:pvar])

fit_Lap <- gss:::mspreg1(s = s, r = Lap.gram, id.basis=1:n.fit, method = 'm',
                         y = Y.all, wt=1, alpha = 1.4, random = NULL, varht=1, skip.iter=FALSE)

theta_l <- fit_Lap$theta
c_l <- matrix(fit_Lap$c, n.fit, 1) 
d_l <- fit_Lap$d

## fit MRD

L <- matrix(0, n.fit, n.fit)
for (i in 1:pall) {
  L <- 10^theta_l[i] * Lap.gram[,,i] + L
}

hat_lap <- c(L %*% c_l + s %*% d_l)
fitMRD.3 <- round(mean(((hat_lap - Y.all)/Y.all)^2), 4)


























######################## Simulation hippocampus data ######################

################ calculation condition ####################
# Simulation needs to be parallel computation, 
# this paper uses 10 cores for parallel computation, 
# the reader can also according to their own needs 
# and hardware conditions to set the number of cores that need to be parallel computation
####################################


f.real <- matrix(0, m, 1)
M <- matrix(0, m, m)
theta.real <- c(0.5, 0.5, 0.8, 0.8, 0, 0, 4, 0, 4, 0) # theta.real <- c(0.5, 0.5, 0.8, 0.8, 0, 0, 100, 0, 50, 0)  
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
  beta_hat_1 <- d_fit
  
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
  beta_hat_2_list <- list()
  
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
    beta_hat_2_list[[h]] <- d_fit
    gcv_list[[h]] <- gcv(fit_gauss.j$nlambda, theta_fit, Y.real, m, Q_fit, s)
  }
  min_gcv <- which.min(gcv_list)
  MSE_02 <- MSE_02_all[[min_gcv]]
  beta_hat_2 <- beta_hat_2_list[[min_gcv]]
  Cos_diag_2 <- Cos_diag_2_all[[min_gcv]]
  
  ######### Laplacian kernel #############
  
  MSE_03_all <- list()
  Cos_diag_3_all <- list()
  gcv_list <- list()
  beta_hat_3_list <- list()
  
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
    beta_hat_3_list[[h]] <- d_fit
    gcv_list[[h]] <- gcv(fit_lap$nlambda, theta_fit, Y.real, m, Q_fit, s)
  }
  min_gcv <- which.min(gcv_list)
  MSE_03 <- MSE_03_all[[min_gcv]]
  beta_hat_3 <- beta_hat_3_list[[min_gcv]]
  Cos_diag_3 <- Cos_diag_3_all[[min_gcv]]
  
  return(round(c(MSE_01, MSE_02, MSE_03, Cos_diag_1, Cos_diag_2, Cos_diag_3, beta_hat_1, beta_hat_2, beta_hat_3), 4))
  
}


set.seed(1)
registerDoMC(25)
## warning! long time runing 
result.NE2 = foreach(rep = 1:200, .combine = cbind) %dopar%{
  a <- fitFun_cv(0.05)
  b <- fitFun_cv(0.10)
  c <- fitFun_cv(0.15)
  c(a,b,c)
}


mean_left4 = rowMeans(result.NE2)
sd_left4 = apply(result.NE2, 1, sd)
mat = matrix(c(c(1:9), c(52:60), c(103:111)), 9, 3, byrow = TRUE)

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


########### estimated error of beta #################

# First repeat d.real 9 times to get a 36x1 matrix
d.real_expanded <- matrix(rep(d.real, 9), nrow = 36, ncol = 1)

# then repeat this 36x1 matrix 200 times to get a 36x200 matrix
d_matrix <- matrix(rep(d.real_expanded, 200), nrow = 36, ncol = 200)

beta_res <- (result.NE2[c(40:51, 91:102, 142:153),] - d_matrix)^2

sum_matrix <- matrix(0, nrow = 9, ncol = 200)

for (i in 1:9) {
  rows <- ((i-1) * 4 + 1):(i * 4)  
  sum_matrix[i, ] <- colSums(beta_res[rows, ])
}

EE_matrix <- rowMeans(sum_matrix^(0.5)) # 9*1 


# Plot the chart.
len_rep = 200
seq1 = factor(c(rep("DG-CA23", len_rep), rep("DG-CA1", len_rep), rep("DG-Subi", len_rep), rep("CA23-CA1", len_rep),
                rep("CA23-Subi", len_rep), rep("CA1-Subi", len_rep)), 
              levels = c("DG-CA23","DG-CA1", "DG-Subi", "CA23-CA1", "CA23-Subi", "CA1-Subi"))
data1 = c(result.NE2[5+9,], result.NE2[6+9,], result.NE2[7+9,], result.NE2[8+9,], result.NE2[9+9,], result.NE2[10+9,])
data2 = c(result.NE2[15+9,], result.NE2[16+9,], result.NE2[17+9,], result.NE2[18+9,], result.NE2[19+9,], result.NE2[20+9,])
data3 = c(result.NE2[25+9,], result.NE2[26+9,], result.NE2[27+9,], result.NE2[28+9,], result.NE2[29+9,], result.NE2[30+9,])
data4 = c(result.NE2[12+35+18,], result.NE2[12+36+18,], result.NE2[12+37+18,], result.NE2[12+38+18,], result.NE2[12+39+18,], result.NE2[12+40+18,])
data5 = c(result.NE2[12+45+18,], result.NE2[12+46+18,], result.NE2[12+47+18,], result.NE2[12+48+18,], result.NE2[12+49+18,], result.NE2[12+50+18,])
data6 = c(result.NE2[12+55+18,], result.NE2[12+56+18,], result.NE2[12+57+18,], result.NE2[12+58+18,], result.NE2[12+59+18,], result.NE2[12+60+18,])
data7 = c(result.NE2[24+65+27,], result.NE2[24+66+27,], result.NE2[24+67+27,], result.NE2[24+68+27,], result.NE2[24+69+27,], result.NE2[24+70+27,])
data8 = c(result.NE2[24+75+27,], result.NE2[24+76+27,], result.NE2[24+77+27,], result.NE2[24+78+27,], result.NE2[24+79+27,], result.NE2[24+80+27,])
data9 = c(result.NE2[24+85+27,], result.NE2[24+86+27,], result.NE2[24+87+27,], result.NE2[24+88+27,], result.NE2[24+89+27,], result.NE2[24+90+27,])
real_cos = c(rep(Cos_diag_sim[5], 200), rep(Cos_diag_sim[6], 200), rep(Cos_diag_sim[7], 200),
             rep(Cos_diag_sim[8], 200), rep(Cos_diag_sim[9], 200), rep(Cos_diag_sim[10], 200))
dataf = data.frame(data1, data2, data3,
                   data4, data5, data6,
                   data7, data8, data9, seq1, real_cos)


p1 = ggplot(dataf, aes(seq1, data1, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") +
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) +
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = element_blank(), y = "Cosine diagnosis values",
       title = expression(paste("Metric kernel(", paste(sigma, " = 0.05)"))))

p2 = ggplot(dataf, aes(seq1, data2, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = element_blank(), y = element_blank(), 
       title = expression(paste("Exponential kernel(", paste(sigma, " = 0.05)")))) +
  theme(axis.text.y = element_blank()) 

p3 = ggplot(dataf, aes(seq1, data3, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = element_blank(), y = element_blank(), 
       title = expression(paste("Laplacian kernel(", paste(sigma, " = 0.05)")))) +
  theme(axis.text.y = element_blank()) 

p4 = ggplot(dataf, aes(seq1, data4, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = element_blank(), y = "Cosine diagnosis values", 
       title = expression(paste("Metric kernel(", paste(sigma, " = 0.1)"))))

p5 = ggplot(dataf, aes(seq1, data5, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = element_blank(), y = element_blank(), 
       title = expression(paste("Exponential kernel(", paste(sigma, " = 0.1)")))) +
  theme(axis.text.y = element_blank()) 

p6 = ggplot(dataf, aes(seq1, data6, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = element_blank(), y = element_blank(), 
       title = expression(paste("Laplacian kernel(", paste(sigma, " = 0.1)")))) +
  theme(axis.text.y = element_blank()) 

p7 = ggplot(dataf, aes(seq1, data7, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = "Subregions", y = "Cosine diagnosis values", 
       title = expression(paste("Metric kernel(", paste(sigma, " = 0.15)"))))

p8 = ggplot(dataf, aes(seq1, data8, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = "Subregions", y = element_blank(), 
       title = expression(paste("Exponential kernel(", paste(sigma, " = 0.15)")))) +
  theme(axis.text.y = element_blank()) 

p9 = ggplot(dataf, aes(seq1, data9, fill = seq1)) +
  geom_violin(trim=TRUE) + geom_jitter(shape = 21, size = 0.5) + 
  theme(legend.position="none") + 
  geom_boxplot(width = 0.2)+coord_cartesian(ylim = c(-0.5, 1.5)) + 
  geom_errorbar(aes(ymax = real_cos, ymin = real_cos), colour = "#AA0000") +
  labs(x = "Subregions", y = element_blank(), 
       title = expression(paste("Laplacian kernel(", paste(sigma, " = 0.15)")))) +
  theme(axis.text.y = element_blank()) 


# The result in Figure 2 is obtained
multiplot(p1, p2, p3, cols=3)
multiplot(p4, p5, p6, cols=3)
multiplot(p7, p8, p9, cols=3)



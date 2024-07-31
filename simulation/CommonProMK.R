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


#####################################################
## Tensor product
Q_m4 <- ker_est_all(N, dist_train, id.basis, N.gss, pvar, h)
#####################################################

s1 = s
fit_m4 <- gss:::mspreg1(s = s1, r = Q_m4[,id.basis,], id.basis=id.basis, method = 'm',
                        y = Y.all,wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE)

theta_m4 <- fit_m4$theta
c_m4 <- matrix(0, N, 1)
c_m4[id.basis,1] <- fit_m4$c
d_m4 <- fit_m4$d


###  fitting
h_K_m <- matrix(0, N, N)
for (i in 1:pall) {
  h_K_m <- 10^theta_m4[i] * Q_m4[, , i] + h_K_m
}

d_m4 <- as.matrix(d_m4) 
hat_f_m4 <- h_K_m %*% c_m4
hat_f_m_p <- h_K_m %*% c_m4 + s1 %*% d_m4


## MSE
mse = round(mean(((hat_f_m_p - Y.all))^2), 4)

cat("MSE is: ", mse, "\n")

#### Test whether the scalar variable is significant ####

nlambda <- 10^fit_m4$nlambda
M <- matrix(0, N, N)
for (i in 1:pall) {
  M = 10^theta_m4[i] * Q_m4[,,i] + M
}

M <- M + nlambda * diag(N)
M_inv = solve(M)

#################################
sigma_hat <- var(Y.all - hat_f_m_p)   #### wrong  var(Y.all - hat_f_m4) 
cov_d <- solve(t(s1) %*% M_inv %*% s1) %*% t(s1) %*% M_inv %*% M_inv %*% s1 %*% solve(t(s1) %*% M_inv %*% s1) * sigma_hat[1,1]
SDest <- round((diag(cov_d))^(0.5),4)
# scalar parameter estimates and sds
est = rbind(c(round(d_m4, 4)), SDest)
dimnames(est)[[1]] = c("est", "SDest")
est

cat("Estimates of scale parameters are: ", "\n")
print(est)
cat("\n")

######################## MRD ##############################

n.fit <- length(Y.all)

## fit MRD

fitMRD.1 <- round(mean(((hat_f_m_p - Y.all)/Y.all)^2), 4)
fitMRD.1

cat("MRD is: ", fitMRD.1, "\n")

########################## cosine diagnosis #####################
Cos_diag <- matrix(0, 1, pall)
for (i in 1:pall) {
  Cos_diag[1, i] <- sum((Q_m4[,,i] %*% c_m4 * 10^theta_m4[i]) * hat_f_m4) / sum(hat_f_m4^2)
}
# Cos_diag <- round(Cos_diag, 4)
pnam = c(1:pvar)
for(i in 1:(pvar-1))
  for(j in (i+1):pvar)
    pnam = c(pnam, paste(i,",",j, sep=""))

est = t(matrix(round(Cos_diag, 4)))
dimnames(est)[[1]] = c("Stat")
dimnames(est)[[2]] = pnam
cat("Cosine diag statistics are: ", "\n")
print(est)
cat("\n")
########################

######################## F test ####################

h_K_m <- matrix(0, N, N)
for (i in 1:pall) {
  h_K_m <- 10^theta_m4[i] * Q_m4[, , i] + h_K_m
}

d_m4 <- as.matrix(d_m4) 
hat_f_m_p <- hat_f_m4 + s1 %*% d_m4

##################### F test for the nonparametric f

fit_H0 = lm(Y.all~s1[,-1])
d_H0 <- fit_H0$coef

###  fitting
hat_f_H0 <- s1 %*% matrix(d_H0)

M_with_hipp <- h_K_m + (nlambda) * diag(N)
M_inv_with = solve(M_with_hipp)
Alam_H0 <- s1 %*% solve(t(s1) %*%  s1) %*% t(s1)
Alam_H1 <- diag(N) - (nlambda)*(M_inv_with - M_inv_with %*% s1 %*%  solve(t(s1) %*% M_inv_with %*% s1) %*% t(s1) %*% M_inv_with)

fhat_H0 <- hat_f_H0   #Alam_H0 %*% Y.all
fhat_H1 <- hat_f_m_p   #Alam_H1 %*% Y.all


RSS_H0 <- sum((Y.all - fhat_H0)^2)
RSS_H1 <- sum((Y.all - fhat_H1)^2)


df_H0 <- N - sum(diag(Alam_H0))
df_H1 <- N - sum(diag(Alam_H1))

F_stat <- (RSS_H0 - RSS_H1) * df_H1 / (df_H0 - df_H1) / RSS_H1
# F test statistic
F_stat

times_boot <- 200
ep <- Y.all - fhat_H1
ep <- ep - mean(ep)
bootstrap_samples <- bootstrap_sample(ep, times_boot)

set.seed(2022)

result_boot = foreach(i = 1:times_boot, .combine = cbind) %dopar% 
  {
    
    y_star <- bootstrap_samples[, i] + fhat_H0  
    
    if(TRUE)
    {
      sel_res = select_rep(N, y_star, pho.r)
      id.basis.boot = sel_res[[1]]
      N.gss = sel_res[[2]]
      
      dist_train <- dist_matrix4
      Q_m4_boot <- ker_est_all(N, dist_train, id.basis.boot, N.gss, pvar, h)
      
      
    } else {
      Q_m4_boot <- Q_m4  
      id.basis.boot <- id.basis
    }
    fit_boot_all <- gss:::mspreg1(s = s1, r = Q_m4_boot[,id.basis.boot,], id.basis=id.basis.boot, method = 'm',
                                  y = y_star,wt=1,alpha = 1.4,random = NULL,varht=1,skip.iter=FALSE)
    
    theta_boot_all <- fit_boot_all$theta
    c_boot_all <- matrix(0, N, 1)
    c_boot_all[id.basis.boot,1] <- fit_boot_all$c
    d_boot_all <- fit_boot_all$d
    
    
    ###  fitting
    h_K_boot_all <- matrix(0, N, N)
    for (i in 1:pall) {
      h_K_boot_all <- 10^theta_boot_all[i] * Q_m4_boot[, , i] + h_K_boot_all
    }
    
    hat_f_boot_all <- h_K_boot_all %*% c_boot_all + s1 %*% d_boot_all
    #### without
    
    fit_boot_H0 = lm(y_star~s1[,-1])
    d_boot_H0 <- fit_boot_H0$coef
    hat_f_boot_H0 <- s1 %*% matrix(d_boot_H0)
    ######
    nlam_boot_all <- 10^fit_boot_all$nlambda 
    
    M_with_CA1_boot <- h_K_boot_all + nlam_boot_all * diag(N)
    M_inv_with_b = solve(M_with_CA1_boot)
    
    Alam_H0_boot <- s1 %*% solve(t(s1) %*%  s1) %*% t(s1)
    Alam_H1_boot <- diag(N) - nlam_boot_all*(M_inv_with_b - M_inv_with_b %*% s1 %*% 
                                               solve(t(s1) %*% M_inv_with_b %*% s1) %*% t(s1) %*% M_inv_with_b)
    
    fhat_H0_boot <- hat_f_boot_H0
    fhat_H1_boot <- hat_f_boot_all
    
    
    RSS_H0_boot <- sum((y_star - fhat_H0_boot)^2)
    RSS_H1_boot <- sum((y_star - fhat_H1_boot)^2)
    
    
    df_H0_boot <- N - sum(diag(Alam_H0_boot))
    df_H1_boot <- N - sum(diag(Alam_H1_boot))
    
    F_stat_boot <- (RSS_H0_boot - RSS_H1_boot) * df_H1_boot / (df_H0_boot - df_H1_boot) / RSS_H1_boot
    F_stat_boot
    
  }

p_value <- sum(c(result_boot) - F_stat > 0) / times_boot

cat("F statistic is: ", F_stat, "\n")
cat("P value is: ", p_value, "\n")
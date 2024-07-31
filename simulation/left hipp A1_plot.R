library(gss)
library(stringr)
library(ggplot2)
library(doMC)
library(foreach)
#library(ggpubr)
#library(tidyverse)
library(SpatialPack)
library(plotly)
library(MASS)
library(boot)
library(plot3D)
library(kernlab)

Rcpp::sourceCpp("Kernelmatrix.cpp")
source("CommonMethods.R")
source("LoadData.R")
registerDoMC(20)


####################################

######################################




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
  dist_matrix4[,,i] <- as.matrix(dist(Xfun[,which(label_left[,2] == c_reg[i])], 
                                      diag = T, upper = T, method = 'euclidean')) 
}

########### Parameter selection ##########
# 
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


################ visualization ###############

##### pi_{3,4}  (CA1 - Subi) #####
i = 3
kernmat = as.kernelMatrix(Q_m4[,,i])
kdiag = diag(kernmat)
simmat = kronecker(matrix(kdiag), t(matrix(kdiag)), FUN = "+") - 2*kernmat
fit = cmdscale(simmat, eig=TRUE, k=2)
X_CA1_1 = fit$points[,1]
X_CA1_2 = fit$points[,2]

i = 4
kernmat = as.kernelMatrix(Q_m4[,,i])
kdiag = diag(kernmat)
simmat = kronecker(matrix(kdiag), t(matrix(kdiag)), FUN = "+") - 2*kernmat
fit = cmdscale(simmat, eig=TRUE, k=2)
X_Subi_1 = fit$points[,1]
X_Subi_2 = fit$points[,2]

#### compute f34 ######
i = 10
fi = Q_m4[,,i] %*% c_m4 * 10^theta_m4[i]

#######################
# plot for f34 vs. the first Coordinate

df_with_data <- data.frame(
  X_CA1_1 = X_CA1_1,
  X_CA1_2 = X_CA1_2,
  X_Subi_1 = X_Subi_1,
  X_Subi_2 = X_Subi_2,
  hat_f_m_p = fi
)

###########################
### 3D plot 

# 进行LOESS拟合，使用X_CA1_mean和X_Subi_mean预测hat_f_m_p

# CA1: 1st Coordinate and Subi: 1st Coordinate

loess_fit <- loess(hat_f_m_p ~ X_CA1_1 + X_Subi_1, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_1 = X_CA1_1,
  X_Subi_1 = X_Subi_1,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_1), max(df_with_data$X_CA1_1), length.out = 50)
y <- seq(min(df_with_data$X_Subi_1), max(df_with_data$X_Subi_1), length.out = 50)
grid <- expand.grid(X_CA1_1 = x, X_Subi_1 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_1, y = ~X_Subi_1, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-First Coordinate", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-First Coordinate", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig

# CA1: 1st Coordinate and Subi: 2nd Coordinate

loess_fit <- loess(hat_f_m_p ~ X_CA1_1 + X_Subi_2, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_1 = X_CA1_1,
  X_Subi_2 = X_Subi_2,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_1), max(df_with_data$X_CA1_1), length.out = 50)
y <- seq(min(df_with_data$X_Subi_2), max(df_with_data$X_Subi_2), length.out = 50)
grid <- expand.grid(X_CA1_1 = x, X_Subi_2 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_1, y = ~X_Subi_2, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-First Coordinate", color = "black", showgrid = TRUE),
      yaxis = list(title = "subi-Second Coordinate", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig

# CA1: 2nd Coordinate and Subi: 2nd Coordinate

loess_fit <- loess(hat_f_m_p ~ X_CA1_2 + X_Subi_2, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_2 = X_CA1_2,
  X_Subi_2 = X_Subi_2,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_2), max(df_with_data$X_CA1_2), length.out = 50)
y <- seq(min(df_with_data$X_Subi_2), max(df_with_data$X_Subi_2), length.out = 50)
grid <- expand.grid(X_CA1_2 = x, X_Subi_2 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_2, y = ~X_Subi_2, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-Second Coordinate", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-Second Coordinate", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig


# CA1: 2nd Coordinate and Subi: 1st Coordinate

loess_fit <- loess(hat_f_m_p ~ X_CA1_2 + X_Subi_1, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_2 = X_CA1_2,
  X_Subi_1 = X_Subi_1,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_2), max(df_with_data$X_CA1_2), length.out = 50)
y <- seq(min(df_with_data$X_Subi_1), max(df_with_data$X_Subi_1), length.out = 50)
grid <- expand.grid(X_CA1_2 = x, X_Subi_1 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_2, y = ~X_Subi_1, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-Second Coordinate", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-First Coordinate", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig

#################### Kernel PCA ##
########################################
##### pi_{3,4} CA1-Subi #####
i = 3
kernmat = as.kernelMatrix(Q_m4[,,i])
kpc <- kpca(kernmat, features=3)
pca.my = pcv(kpc)
train.my = t(t(pca.my)%*%t(kernmat))
X_CA1_1 = train.my[,1]
X_CA1_2 = train.my[,2]

i = 4
kernmat = as.kernelMatrix(Q_m4[,,i])
kpc <- kpca(kernmat, features=3)
pca.my = pcv(kpc)
train.my = t(t(pca.my)%*%t(kernmat))
X_Subi_1 = train.my[,1]
X_Subi_2 = train.my[,2]

#### compute f12 ######
i = 10
fi = Q_m4[,,i] %*% c_m4 * 10^theta_m4[i]

#######################
# plot for f34 vs. the first pc

df_with_data <- data.frame(
  X_CA1_1 = X_CA1_1,
  X_CA1_2 = X_CA1_2,
  X_Subi_1 = X_Subi_1,
  X_Subi_2 = X_Subi_2,
  hat_f_m_p = fi
)

###########################
### 3D plot 

# 进行LOESS拟合，使用X_CA1_mean和X_Subi_mean预测hat_f_m_p

# Para: 1st PC and Pre: 1st PC

loess_fit <- loess(hat_f_m_p ~ X_CA1_1 + X_Subi_1, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_1 = X_CA1_1,
  X_Subi_1 = X_Subi_1,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_1), max(df_with_data$X_CA1_1), length.out = 50)
y <- seq(min(df_with_data$X_Subi_1), max(df_with_data$X_Subi_1), length.out = 50)
grid <- expand.grid(X_CA1_1 = x, X_Subi_1 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_1, y = ~X_Subi_1, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-First PC", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-First PC", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig

# CA1: 1st PC and Subi: 2nd PC

loess_fit <- loess(hat_f_m_p ~ X_CA1_1 + X_Subi_2, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_1 = X_CA1_1,
  X_Subi_2 = X_Subi_2,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_1), max(df_with_data$X_CA1_1), length.out = 50)
y <- seq(min(df_with_data$X_Subi_2), max(df_with_data$X_Subi_2), length.out = 50)
grid <- expand.grid(X_CA1_1 = x, X_Subi_2 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_1, y = ~X_Subi_2, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-First PC", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-Second PC", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig

# CA1: 2nd PC and Subi: 2nd PC

loess_fit <- loess(hat_f_m_p ~ X_CA1_2 + X_Subi_2, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_2 = X_CA1_2,
  X_Subi_2 = X_Subi_2,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_2), max(df_with_data$X_CA1_2), length.out = 50)
y <- seq(min(df_with_data$X_Subi_2), max(df_with_data$X_Subi_2), length.out = 50)
grid <- expand.grid(X_CA1_2 = x, X_Subi_2 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_2, y = ~X_Subi_2, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-Second PC", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-Second PC", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig


# CA1: 2nd PC and Subi: 1st PC

loess_fit <- loess(hat_f_m_p ~ X_CA1_2 + X_Subi_1, span = 0.2, data = df_with_data)
fitted_values <- predict(loess_fit)

df_with_fit <- data.frame(
  X_CA1_2 = X_CA1_2,
  X_Subi_1 = X_Subi_1,
  fitted_values = fitted_values
)

# 创建拟合曲面数据
x <- seq(min(df_with_data$X_CA1_2), max(df_with_data$X_CA1_2), length.out = 50)
y <- seq(min(df_with_data$X_Subi_1), max(df_with_data$X_Subi_1), length.out = 50)
grid <- expand.grid(X_CA1_2 = x, X_Subi_1 = y)
z <- t(matrix(predict(loess_fit, newdata = grid),ncol = 50, nrow = 50))


fig <- plot_ly(x = ~x, y = ~y, z = ~ z, type = "surface",color="red", alpha = 0.5, showscale = FALSE) %>%
  add_markers(data = df_with_data, 
              x = ~X_CA1_2, y = ~X_Subi_1, z = ~ hat_f_m_p,
              marker = list(size = 3, color = 'black')) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "CA1-Second PC", color = "black", showgrid = TRUE),
      yaxis = list(title = "Subi-First PC", color = "black", showgrid = TRUE),
      zaxis = list(title = "Fit", color = "black", showgrid = TRUE)
    )
  )

fig
#############################
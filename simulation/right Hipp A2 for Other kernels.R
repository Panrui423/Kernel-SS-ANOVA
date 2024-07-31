library(gss)
library(stringr)
library(ggplot2)
library(doMC)
library(foreach)
library(SpatialPack)
library(plotly)
library(MASS)
library(boot)
library(plot3D)

Rcpp::sourceCpp("Kernelmatrix.cpp")
source("CommonMethods.R")
registerDoMC(20)

############################# Right Analysis 2 ##################################
##############################  Gaussian kernel method #########################

source("LoadData.R")
c_reg = c(1, 2, 3) 
c_l = length(c_reg)

N = length(Y.all)
Y = Y.all
pvar = length(c_reg)
pall = pvar + pvar*(pvar-1)/2
h = 0.1
id.basis = 1:N
#####################################################
## Tensor product

#############  Gaussian Kernel #############
Q_all <- array(0, dim = c(N, N, pvar))
for (k in 1:pvar) {
  for (i in 1:N) {
    for (j in 1:N) {
      Q_all[i,j, k] = kernalGaussian(Xfun[i, which(label_right[,2] == c_reg[k])], 
                                     Xfun[j, which(label_right[,2] == c_reg[k])], para = h)
    }
  }
}

Q_m4 <- Tensor_f(Q_all)
#####################################################
cat("Results with Gaussian kernel are as follow:","\n")
source("CommonProGL.R")


######################################################
########################## Laplace Kernel   method  ###########

source("LoadData.R")
c_reg = c(1, 2, 3)
c_l = length(c_reg)

N = length(Y.all)
Y = Y.all
pvar = length(c_reg)
pall = pvar + pvar*(pvar-1)/2
h = 0.05
id.basis = 1:N
#####################################################
## Tensor product

#############  Laplace Kernel #############
Q_all <- array(0, dim = c(N, N, pvar))
for (k in 1:pvar) {
  for (i in 1:N) {
    for (j in 1:N) {
      Q_all[i,j, k] = kernalLap(Xfun[i, which(label_right[,2] == c_reg[k])], 
                                Xfun[j, which(label_right[,2] == c_reg[k])], para = h)
    }
  }
}

Q_m4 <- Tensor_f(Q_all)
######################################
cat("Results with Laplace kernel are as follow:","\n")
source("CommonProGL.R")





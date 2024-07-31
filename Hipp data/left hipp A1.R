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

Rcpp::sourceCpp("Kernelmatrix.cpp")
source("CommonMethods.R")
source("LoadData.R")
registerDoMC(20)


####################################

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


source("CommonProMK.R")

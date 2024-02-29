library(gss)
library(convexjlr)
library("doMC")
library("foreach")
library(xtable)
source("Manifold_simulation.R")
source("Klein_simulation.R")

###########  #############
# Simulation needs to be parallel computation, 
# this paper uses 10 cores for parallel computation, 
# the reader can also according to their own needs 
# and hardware conditions to set the number of cores that need to be parallel computation
########################


# warning! It takes about a day to run.
#################### sigma = 0.1^(1/2) ################
set.seed(2021)
registerDoMC(10)

result = foreach(i.rep = 1:200, .combine = cbind) %dopar%
{ 
  
  sim_m1 = simulation_M(n.fit = 200, n.test = 100, h=0.1, 
                        scale=c(1/4), sigma=0.1^(1/2))
  sim_m2 = simulation_M(n.fit = 300, n.test = 100, h=0.1, 
                        scale=c(1/4), sigma=0.1^(1/2))
  sim_m3 = simulation_M(n.fit = 400, n.test = 100, h=0.1, 
                        scale=c(1/4), sigma=0.1^(1/2))

  sim_k1 = simulation_K(n.fit = 200, n.test = 100, h=1, 
                        scale=c(1/4), sigma=0.1^(1/2))
  sim_k2 = simulation_K(n.fit = 300, n.test = 100, h=1, 
                        scale=c(1/4), sigma=0.1^(1/2))
  sim_k3 = simulation_K(n.fit = 400, n.test = 100, h=1, 
                        scale=c(1/4), sigma=0.1^(1/2))
  
  rbind(sim_m1, sim_m2, sim_m3,
        sim_k1, sim_k2, sim_k3)
}


cmat <- matrix(c(1:48), 16, 3)

row_means = round(rowMeans(result, na.rm=T), 4)

### three-dimensional vector
M_meantimes_1 = row_means[t(cmat[1,])]
M_meantimes_2 = row_means[t(cmat[2,])]
M_predict_error_m1 = row_means[t(cmat[3,])]
M_fit_error_m1 = row_means[t(cmat[4,])]
M_predict_error_m2 = row_means[t(cmat[5,])]
M_fit_error_m2 = row_means[t(cmat[6,])]

M_meantimes_1_no = row_means[cmat[7,]]
M_meantimes_2_no = row_means[cmat[8,]]
M_predict_error_m1_no = row_means[cmat[9,]]
M_fit_error_m1_no = row_means[cmat[10,]]
M_predict_error_m2_no = row_means[cmat[11,]]
M_fit_error_m2_no = row_means[cmat[12,]]

M_fit_error_m3 = row_means[cmat[13,]]
M_predict_error_m3 = row_means[cmat[14,]]
M_fit_error_m4 = row_means[cmat[15,]]
M_predict_error_m4 = row_means[cmat[16,]]

K_meantimes_1 = row_means[t(cmat[1,] + 48)]
K_meantimes_2 = row_means[t(cmat[2,] + 48)]
K_predict_error_m1 = row_means[t(cmat[3,] + 48)]
K_fit_error_m1 = row_means[t(cmat[4,] + 48)]
K_predict_error_m2 = row_means[t(cmat[5,] + 48)]
K_fit_error_m2 = row_means[t(cmat[6,] + 48)]

K_meantimes_1_no = row_means[cmat[7,] + 48]
K_meantimes_2_no = row_means[cmat[8,] + 48]
K_predict_error_m1_no = row_means[cmat[9,] + 48]
K_fit_error_m1_no = row_means[cmat[10,] + 48]
K_predict_error_m2_no = row_means[cmat[11,] + 48]
K_fit_error_m2_no = row_means[cmat[12,] + 48]


K_fit_error_m3 = row_means[cmat[13,] + 48]
K_predict_error_m3 = row_means[cmat[14,] + 48]
K_fit_error_m4 = row_means[cmat[15,] + 48]
K_predict_error_m4 = row_means[cmat[16,] + 48]


########
sim_sd = round(apply(result, 1, sd, na.rm=T), 4)

M_predict_sd_m1 = sim_sd[t(cmat[3,])]
M_fit_sd_m1 = sim_sd[t(cmat[4,])]
M_predict_sd_m2 = sim_sd[t(cmat[5,])]
M_fit_sd_m2 = sim_sd[t(cmat[6,])]

M_predict_sd_m1_no = sim_sd[cmat[9,]]
M_fit_sd_m1_no = sim_sd[cmat[10,]]
M_predict_sd_m2_no = sim_sd[cmat[11,]]
M_fit_sd_m2_no = sim_sd[cmat[12,]]

M_fit_sd_m3 = sim_sd[cmat[13,]]
M_predict_sd_m3 = sim_sd[cmat[14,]]
M_fit_sd_m4 = sim_sd[cmat[15,]]
M_predict_sd_m4 = sim_sd[cmat[16,]]

K_predict_sd_m1 = sim_sd[t(cmat[3,]+48)]
K_fit_sd_m1 = sim_sd[t(cmat[4,]+48)]
K_predict_sd_m2 = sim_sd[t(cmat[5,]+48)]
K_fit_sd_m2 = sim_sd[t(cmat[6,]+48)]

K_predict_sd_m1_no = sim_sd[cmat[9,]+48]
K_fit_sd_m1_no = sim_sd[cmat[10,]+48]
K_predict_sd_m2_no = sim_sd[cmat[11,]+48]
K_fit_sd_m2_no = sim_sd[cmat[12,]+48]

K_fit_sd_m3 = sim_sd[cmat[13,]+48]
K_predict_sd_m3 = sim_sd[cmat[14,]+48]
K_fit_sd_m4 = sim_sd[cmat[15,]+48]
K_predict_sd_m4 = sim_sd[cmat[16,]+48]



data.f1 <- data.frame(Model = matrix(c(rep('M1',3), rep('M2',3)), 6, 1),
                      
                      Sample.size = matrix(c(rep(c(200,300,400),2)), 6, 1),
 
                      MSE.predict_M2_4 = matrix(c(c(M_predict_error_m2[1:3]), c(K_predict_error_m2[1:3])), 6, 1),
                      
                      MSE.predict_M2_sd_4 = matrix(c(c(M_predict_sd_m2[1:3]), c(K_predict_sd_m2[1:3])), 6, 1),
                      
                      MSE.predict_M2_no = matrix(c(c(M_predict_error_m2_no), c(K_predict_error_m2_no)), 6, 1),
                      
                      MSE.predict_M2_no_sd = matrix(c(c(M_predict_sd_m2_no), c(K_predict_sd_m2_no)), 6, 1),
                      
                      MSE.predict_M3 = matrix(c(c(M_predict_error_m3), c(K_predict_error_m3)), 6, 1),
                      
                      MSE.predict_M3_sd = matrix(c(c(M_predict_sd_m3), c(K_predict_sd_m3)), 6, 1),
                      
                      MSE.predict_M4 = matrix(c(c(M_predict_error_m4), c(K_predict_error_m4)), 6, 1),
                      
                      MSE.predict_M4_sd = matrix(c(c(M_predict_sd_m4), c(K_predict_sd_m4)), 6, 1),
                      
                      stringsAsFactors=FALSE)

xtable(data.f1,digits=4)

# warning! It takes about a day to run.
#################### sigma = 0.3^(1/2) ################
set.seed(2021)
registerDoMC(10)
result_sigma_3 = foreach(i.rep = 1:200, .combine = cbind) %dopar%
  { 
    
    sim_m1 = simulation_M(n.fit = 200, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.3^(1/2))
    sim_m2 = simulation_M(n.fit = 300, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.3^(1/2))
    sim_m3 = simulation_M(n.fit = 400, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.3^(1/2))
    
    sim_k1 = simulation_K(n.fit = 200, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.3^(1/2))
    sim_k2 = simulation_K(n.fit = 300, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.3^(1/2))
    sim_k3 = simulation_K(n.fit = 400, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.3^(1/2))

    
    rbind(sim_m1, sim_m2, sim_m3,
          sim_k1, sim_k2, sim_k3)
  }


cmat <- matrix(c(1:48), 16, 3)

row_means = round(rowMeans(result_sigma_3, na.rm=T), 4)

### three-dimensional vector
M_meantimes_1 = row_means[t(cmat[1,])]
M_meantimes_2 = row_means[t(cmat[2,])]
M_predict_error_m1 = row_means[t(cmat[3,])]
M_fit_error_m1 = row_means[t(cmat[4,])]
M_predict_error_m2 = row_means[t(cmat[5,])]
M_fit_error_m2 = row_means[t(cmat[6,])]

M_meantimes_1_no = row_means[cmat[7,]]
M_meantimes_2_no = row_means[cmat[8,]]
M_predict_error_m1_no = row_means[cmat[9,]]
M_fit_error_m1_no = row_means[cmat[10,]]
M_predict_error_m2_no = row_means[cmat[11,]]
M_fit_error_m2_no = row_means[cmat[12,]]

M_fit_error_m3 = row_means[cmat[13,]]
M_predict_error_m3 = row_means[cmat[14,]]
M_fit_error_m4 = row_means[cmat[15,]]
M_predict_error_m4 = row_means[cmat[16,]]




K_meantimes_1 = row_means[t(cmat[1,] + 48)]
K_meantimes_2 = row_means[t(cmat[2,] + 48)]
K_predict_error_m1 = row_means[t(cmat[3,] + 48)]
K_fit_error_m1 = row_means[t(cmat[4,] + 48)]
K_predict_error_m2 = row_means[t(cmat[5,] + 48)]
K_fit_error_m2 = row_means[t(cmat[6,] + 48)]

K_meantimes_1_no = row_means[cmat[7,] + 48]
K_meantimes_2_no = row_means[cmat[8,] + 48]
K_predict_error_m1_no = row_means[cmat[9,] + 48]
K_fit_error_m1_no = row_means[cmat[10,] + 48]
K_predict_error_m2_no = row_means[cmat[11,] + 48]
K_fit_error_m2_no = row_means[cmat[12,] + 48]


K_fit_error_m3 = row_means[cmat[13,] + 48]
K_predict_error_m3 = row_means[cmat[14,] + 48]
K_fit_error_m4 = row_means[cmat[15,] + 48]
K_predict_error_m4 = row_means[cmat[16,] + 48]


sim_sd = round(apply(result_sigma_3, 1, sd, na.rm=T), 4)

M_predict_sd_m1 = sim_sd[t(cmat[3,])]
M_fit_sd_m1 = sim_sd[t(cmat[4,])]
M_predict_sd_m2 = sim_sd[t(cmat[5,])]
M_fit_sd_m2 = sim_sd[t(cmat[6,])]

M_predict_sd_m1_no = sim_sd[cmat[9,]]
M_fit_sd_m1_no = sim_sd[cmat[10,]]
M_predict_sd_m2_no = sim_sd[cmat[11,]]
M_fit_sd_m2_no = sim_sd[cmat[12,]]

M_fit_sd_m3 = sim_sd[cmat[13,]]
M_predict_sd_m3 = sim_sd[cmat[14,]]
M_fit_sd_m4 = sim_sd[cmat[15,]]
M_predict_sd_m4 = sim_sd[cmat[16,]]

K_predict_sd_m1 = sim_sd[t(cmat[3,]+48)]
K_fit_sd_m1 = sim_sd[t(cmat[4,]+48)]
K_predict_sd_m2 = sim_sd[t(cmat[5,]+48)]
K_fit_sd_m2 = sim_sd[t(cmat[6,]+48)]

K_predict_sd_m1_no = sim_sd[cmat[9,]+48]
K_fit_sd_m1_no = sim_sd[cmat[10,]+48]
K_predict_sd_m2_no = sim_sd[cmat[11,]+48]
K_fit_sd_m2_no = sim_sd[cmat[12,]+48]

K_fit_sd_m3 = sim_sd[cmat[13,]+48]
K_predict_sd_m3 = sim_sd[cmat[14,]+48]
K_fit_sd_m4 = sim_sd[cmat[15,]+48]
K_predict_sd_m4 = sim_sd[cmat[16,]+48]


data.f1 <- data.frame(Model = matrix(c(rep('M1',3), rep('M2',3)), 6, 1),
                      
                      Sample.size = matrix(c(rep(c(200,300,400),2)), 6, 1),
                               
                      MSE.predict_M2_4 = matrix(c(c(M_predict_error_m2[1:3]), c(K_predict_error_m2[1:3])), 6, 1),
                      
                      MSE.predict_M2_sd_4 = matrix(c(c(M_predict_sd_m2[1:3]), c(K_predict_sd_m2[1:3])), 6, 1),
                      
                      MSE.predict_M2_no = matrix(c(c(M_predict_error_m2_no), c(K_predict_error_m2_no)), 6, 1),
                      
                      MSE.predict_M2_no_sd = matrix(c(c(M_predict_sd_m2_no), c(K_predict_sd_m2_no)), 6, 1),
                      
                      MSE.predict_M3 = matrix(c(c(M_predict_error_m3), c(K_predict_error_m3)), 6, 1),
                      
                      MSE.predict_M3_sd = matrix(c(c(M_predict_sd_m3), c(K_predict_sd_m3)), 6, 1),
                      
                      MSE.predict_M4 = matrix(c(c(M_predict_error_m4), c(K_predict_error_m4)), 6, 1),
                      
                      MSE.predict_M4_sd = matrix(c(c(M_predict_sd_m4), c(K_predict_sd_m4)), 6, 1),
                      
                      stringsAsFactors=FALSE)

xtable(data.f1,digits=4)







# warning! It takes about a day to run.
#################### sigma = 0.5^(1/2) ################
set.seed(2021)
registerDoMC(10)
result_sigma_5 = foreach(i.rep = 1:200, .combine = cbind) %dopar%
  { 
    
    sim_m1 = simulation_M(n.fit = 200, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.5^(1/2))
    sim_m2 = simulation_M(n.fit = 300, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.5^(1/2))
    sim_m3 = simulation_M(n.fit = 400, n.test = 100, h=0.1, 
                          scale=c(1/4), sigma=0.5^(1/2))
    
    sim_k1 = simulation_K(n.fit = 200, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.5^(1/2))
    sim_k2 = simulation_K(n.fit = 300, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.5^(1/2))
    sim_k3 = simulation_K(n.fit = 400, n.test = 100, h=1, 
                          scale=c(1/4), sigma=0.5^(1/2))

    rbind(sim_m1, sim_m2, sim_m3,
          sim_k1, sim_k2, sim_k3)
  }


cmat <- matrix(c(1:48), 16, 3)

row_means = round(rowMeans(result_sigma_5, na.rm=T), 4)

### three-dimensional vector
M_meantimes_1 = row_means[t(cmat[1,])]
M_meantimes_2 = row_means[t(cmat[2,])]
M_predict_error_m1 = row_means[t(cmat[3,])]
M_fit_error_m1 = row_means[t(cmat[4,])]
M_predict_error_m2 = row_means[t(cmat[5,])]
M_fit_error_m2 = row_means[t(cmat[6,])]

M_meantimes_1_no = row_means[cmat[7,]]
M_meantimes_2_no = row_means[cmat[8,]]
M_predict_error_m1_no = row_means[cmat[9,]]
M_fit_error_m1_no = row_means[cmat[10,]]
M_predict_error_m2_no = row_means[cmat[11,]]
M_fit_error_m2_no = row_means[cmat[12,]]

M_fit_error_m3 = row_means[cmat[13,]]
M_predict_error_m3 = row_means[cmat[14,]]
M_fit_error_m4 = row_means[cmat[15,]]
M_predict_error_m4 = row_means[cmat[16,]]




K_meantimes_1 = row_means[t(cmat[1,] + 48)]
K_meantimes_2 = row_means[t(cmat[2,] + 48)]
K_predict_error_m1 = row_means[t(cmat[3,] + 48)]
K_fit_error_m1 = row_means[t(cmat[4,] + 48)]
K_predict_error_m2 = row_means[t(cmat[5,] + 48)]
K_fit_error_m2 = row_means[t(cmat[6,] + 48)]

K_meantimes_1_no = row_means[cmat[7,] + 48]
K_meantimes_2_no = row_means[cmat[8,] + 48]
K_predict_error_m1_no = row_means[cmat[9,] + 48]
K_fit_error_m1_no = row_means[cmat[10,] + 48]
K_predict_error_m2_no = row_means[cmat[11,] + 48]
K_fit_error_m2_no = row_means[cmat[12,] + 48]


K_fit_error_m3 = row_means[cmat[13,] + 48]
K_predict_error_m3 = row_means[cmat[14,] + 48]
K_fit_error_m4 = row_means[cmat[15,] + 48]
K_predict_error_m4 = row_means[cmat[16,] + 48]


sim_sd = round(apply(result_sigma_5, 1, sd, na.rm=T), 4)

M_predict_sd_m1 = sim_sd[t(cmat[3,])]
M_fit_sd_m1 = sim_sd[t(cmat[4,])]
M_predict_sd_m2 = sim_sd[t(cmat[5,])]
M_fit_sd_m2 = sim_sd[t(cmat[6,])]

M_predict_sd_m1_no = sim_sd[cmat[9,]]
M_fit_sd_m1_no = sim_sd[cmat[10,]]
M_predict_sd_m2_no = sim_sd[cmat[11,]]
M_fit_sd_m2_no = sim_sd[cmat[12,]]

M_fit_sd_m3 = sim_sd[cmat[13,]]
M_predict_sd_m3 = sim_sd[cmat[14,]]
M_fit_sd_m4 = sim_sd[cmat[15,]]
M_predict_sd_m4 = sim_sd[cmat[16,]]

K_predict_sd_m1 = sim_sd[t(cmat[3,]+48)]
K_fit_sd_m1 = sim_sd[t(cmat[4,]+48)]
K_predict_sd_m2 = sim_sd[t(cmat[5,]+48)]
K_fit_sd_m2 = sim_sd[t(cmat[6,]+48)]

K_predict_sd_m1_no = sim_sd[cmat[9,]+48]
K_fit_sd_m1_no = sim_sd[cmat[10,]+48]
K_predict_sd_m2_no = sim_sd[cmat[11,]+48]
K_fit_sd_m2_no = sim_sd[cmat[12,]+48]

K_fit_sd_m3 = sim_sd[cmat[13,]+48]
K_predict_sd_m3 = sim_sd[cmat[14,]+48]
K_fit_sd_m4 = sim_sd[cmat[15,]+48]
K_predict_sd_m4 = sim_sd[cmat[16,]+48]


data.f1 <- data.frame(Model = matrix(c(rep('M1',3), rep('M2',3)), 6, 1),
                      
                      Sample.size = matrix(c(rep(c(200,300,400),2)), 6, 1),
                                       
                      MSE.predict_M2_4 = matrix(c(c(M_predict_error_m2[1:3]), c(K_predict_error_m2[1:3])), 6, 1),
                      
                      MSE.predict_M2_sd_4 = matrix(c(c(M_predict_sd_m2[1:3]), c(K_predict_sd_m2[1:3])), 6, 1),
                      
                      MSE.predict_M2_no = matrix(c(c(M_predict_error_m2_no), c(K_predict_error_m2_no)), 6, 1),
                      
                      MSE.predict_M2_no_sd = matrix(c(c(M_predict_sd_m2_no), c(K_predict_sd_m2_no)), 6, 1),
                      
                      MSE.predict_M3 = matrix(c(c(M_predict_error_m3), c(K_predict_error_m3)), 6, 1),
                      
                      MSE.predict_M3_sd = matrix(c(c(M_predict_sd_m3), c(K_predict_sd_m3)), 6, 1),
                      
                      MSE.predict_M4 = matrix(c(c(M_predict_error_m4), c(K_predict_error_m4)), 6, 1),
                      
                      MSE.predict_M4_sd = matrix(c(c(M_predict_sd_m4), c(K_predict_sd_m4)), 6, 1),
                      
                      stringsAsFactors=FALSE)

xtable(data.f1,digits=4)











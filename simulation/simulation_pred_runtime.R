library(gss)
library("doMC")
library("foreach")
source("Manifold_simulation.R")
source("Klein_simulation.R")

###########  #############
# Simulation needs to be parallel computation, 
# this paper uses 20 cores for parallel computation, 
# the reader can also according to their own needs 
# and hardware conditions to set the number of cores that need to be parallel computation
########################


# warning! It takes about one week to run.

#################################################################
set.seed(2022)
registerDoMC(20)

result300 = foreach(i.rep = 1:200, .combine = cbind) %dopar%
  { 
    
     sim_m1 = simulation_M(n.fit = 300, n.test = 100, h=0.1, scale=c(1/4), sigma=0.05)
     sim_m2 = simulation_M(n.fit = 300, n.test = 100, h=0.1, scale=c(1/4), sigma=0.10)
     sim_m3 = simulation_M(n.fit = 300, n.test = 100, h=0.1, scale=c(1/4), sigma=0.15)
    
     sim_k1 = simulation_K(n.fit = 300, n.test = 100, h=1, scale=c(1/4), sigma=0.05)
     sim_k2 = simulation_K(n.fit = 300, n.test = 100, h=1, scale=c(1/4), sigma=0.10)
     sim_k3 = simulation_K(n.fit = 300, n.test = 100, h=1, scale=c(1/4), sigma=0.15)
    
    rbind(sim_m1, sim_m2, sim_m3, sim_k1, sim_k2, sim_k3)
  }

save(result300,file="simu_runtime300.Rdata")

result = result300
matrix(round(apply(result, 1, function(x){mean(x)}), 4),12,)



############################
set.seed(2022)
result200 = foreach(i.rep = 1:200, .combine = cbind) %dopar%
  { 
    
    sim_m1 = simulation_M(n.fit = 200, n.test = 100, h=0.1, scale=c(1/4), sigma=0.05)
    sim_m2 = simulation_M(n.fit = 200, n.test = 100, h=0.1, scale=c(1/4), sigma=0.10)
    sim_m3 = simulation_M(n.fit = 200, n.test = 100, h=0.1, scale=c(1/4), sigma=0.15)
    
    sim_k1 = simulation_K(n.fit = 200, n.test = 100, h=1, scale=c(1/4), sigma=0.05)
    sim_k2 = simulation_K(n.fit = 200, n.test = 100, h=1, scale=c(1/4), sigma=0.10)
    sim_k3 = simulation_K(n.fit = 200, n.test = 100, h=1, scale=c(1/4), sigma=0.15)
    
    rbind(sim_m1, sim_m2, sim_m3, sim_k1, sim_k2, sim_k3)
  }

save(result200,file="simu_runtime200.Rdata")


result = result200
matrix(round(apply(result, 1, function(x){mean(x)}), 4),12,)


####################################################
set.seed(2022)
result400 = foreach(i.rep = 1:200, .combine = cbind) %dopar%
  { 
    
    sim_m1 = simulation_M(n.fit = 400, n.test = 100, h=0.1, scale=c(1/4), sigma=0.05)
    sim_m2 = simulation_M(n.fit = 400, n.test = 100, h=0.1, scale=c(1/4), sigma=0.10)
    sim_m3 = simulation_M(n.fit = 400, n.test = 100, h=0.1, scale=c(1/4), sigma=0.15)
    
    sim_k1 = simulation_K(n.fit = 400, n.test = 100, h=1, scale=c(1/4), sigma=0.05)
    sim_k2 = simulation_K(n.fit = 400, n.test = 100, h=1, scale=c(1/4), sigma=0.10)
    sim_k3 = simulation_K(n.fit = 400, n.test = 100, h=1, scale=c(1/4), sigma=0.15)
   
    rbind(sim_m1, sim_m2, sim_m3, sim_k1, sim_k2, sim_k3)
  }


save(result400,file="simu_runtime400.Rdata")

result = result400

matrix(round(apply(result, 1, function(x){mean(x,na.rm = TRUE)}), 4),12,)
#################################### 
load("simu_runtime200.Rdata")
load("simu_runtime300.Rdata")
load("simu_runtime400.Rdata")


result = result200
est200 = matrix(round(apply(result, 1, function(x){mean(x,na.rm = TRUE)}), 4),12,)
sd200 = matrix(round(apply(result, 1, function(x){sd(x,na.rm = TRUE)}), 4),12,)

result = result300
est300 = matrix(round(apply(result, 1, function(x){mean(x,na.rm = TRUE)}), 4),12,)
sd300 = matrix(round(apply(result, 1, function(x){sd(x,na.rm = TRUE)}), 4),12,)

result = result400
result[11, ] = result[11, ]*60
tmp = result[12, ] 
result[12, tmp < 0.5] = result[12, tmp < 0.5] * 60

result[23, ] = result[23, ]*60
tmp = result[24, ] 
result[24, tmp < 0.5] = result[24, tmp < 0.5] * 60

result[35, ] = result[35, ]*60
tmp = result[36, ] 
result[36, tmp < 0.5] = result[36, tmp < 0.5] * 60

result[47, ] = result[47, ]*60
tmp = result[48, ] 
result[48, tmp < 0.5] = result[48, tmp < 0.5] * 60

result[59, ] = result[59, ]*60
tmp = result[60, ] 
result[60, tmp < 0.5] = result[60, tmp < 0.5] * 60

result[71, ] = result[71, ]*60
tmp = result[72, ] 
result[72, tmp < 0.5] = result[72, tmp < 0.5] * 60

est400 = matrix(round(apply(result, 1, function(x){mean(x,na.rm = TRUE)}), 4),12,)
sd400 = matrix(round(apply(result, 1, function(x){sd(x,na.rm = TRUE)}), 4),12,)


resu = NULL


i=1
tmp1 = cbind(est200[, i], est300[, i], est400[, i])
tmp2 = cbind(sd200[, i], sd300[, i], sd400[, i])

resu = rbind(resu, paste(0.05, "& SO(3) &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[3, 1], "(", tmp2[3, 1],")&", tmp1[5, 1], "(", tmp2[5, 1],")&",
                         tmp1[7, 1], "(", tmp2[7, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                         tmp1[3, 2], "(", tmp2[3, 2],")&", tmp1[5, 2], "(", tmp2[5, 2],")&",
                         tmp1[7, 2], "(", tmp2[7, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",
                          tmp1[3, 3], "(", tmp2[3, 3],")&", tmp1[5, 3], "(", tmp2[5, 3],")&",
                          tmp1[7, 3], "(", tmp2[7, 3],")\\",sep=""))



i=4
tmp1 = cbind(est200[, i], est300[, i], est400[, i])
tmp2 = cbind(sd200[, i], sd300[, i], sd400[, i])

resu = rbind(resu, paste("& KB &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[3, 1], "(", tmp2[3, 1],")&", tmp1[5, 1], "(", tmp2[5, 1],")&",
                         tmp1[7, 1], "(", tmp2[7, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[3, 2], "(", tmp2[3, 2],")&", tmp1[5, 2], "(", tmp2[5, 2],")&",
                          tmp1[7, 2], "(", tmp2[7, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",
                          tmp1[3, 3], "(", tmp2[3, 3],")&", tmp1[5, 3], "(", tmp2[5, 3],")&",
                          tmp1[7, 3], "(", tmp2[7, 3],")\\",sep=""))



i=2
tmp1 = cbind(est200[, i], est300[, i], est400[, i])
tmp2 = cbind(sd200[, i], sd300[, i], sd400[, i])

resu = rbind(resu, paste(0.10, "& SO(3) &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[3, 1], "(", tmp2[3, 1],")&", tmp1[5, 1], "(", tmp2[5, 1],")&",
                         tmp1[7, 1], "(", tmp2[7, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[3, 2], "(", tmp2[3, 2],")&", tmp1[5, 2], "(", tmp2[5, 2],")&",
                          tmp1[7, 2], "(", tmp2[7, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",
                          tmp1[3, 3], "(", tmp2[3, 3],")&", tmp1[5, 3], "(", tmp2[5, 3],")&",
                          tmp1[7, 3], "(", tmp2[7, 3],")\\",sep=""))



i=5
tmp1 = cbind(est200[, i], est300[, i], est400[, i])
tmp2 = cbind(sd200[, i], sd300[, i], sd400[, i])

resu = rbind(resu, paste("& KB &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[3, 1], "(", tmp2[3, 1],")&", tmp1[5, 1], "(", tmp2[5, 1],")&",
                         tmp1[7, 1], "(", tmp2[7, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[3, 2], "(", tmp2[3, 2],")&", tmp1[5, 2], "(", tmp2[5, 2],")&",
                          tmp1[7, 2], "(", tmp2[7, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",
                          tmp1[3, 3], "(", tmp2[3, 3],")&", tmp1[5, 3], "(", tmp2[5, 3],")&",
                          tmp1[7, 3], "(", tmp2[7, 3],")\\",sep=""))


i=3
tmp1 = cbind(est200[, i], est300[, i], est400[, i])
tmp2 = cbind(sd200[, i], sd300[, i], sd400[, i])

resu = rbind(resu, paste(0.15, "& SO(3) &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[3, 1], "(", tmp2[3, 1],")&", tmp1[5, 1], "(", tmp2[5, 1],")&",
                         tmp1[7, 1], "(", tmp2[7, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[3, 2], "(", tmp2[3, 2],")&", tmp1[5, 2], "(", tmp2[5, 2],")&",
                          tmp1[7, 2], "(", tmp2[7, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",
                          tmp1[3, 3], "(", tmp2[3, 3],")&", tmp1[5, 3], "(", tmp2[5, 3],")&",
                          tmp1[7, 3], "(", tmp2[7, 3],")\\",sep=""))



i=6
tmp1 = cbind(est200[, i], est300[, i], est400[, i])
tmp2 = cbind(sd200[, i], sd300[, i], sd400[, i])

resu = rbind(resu, paste("& KB &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[3, 1], "(", tmp2[3, 1],")&", tmp1[5, 1], "(", tmp2[5, 1],")&",
                         tmp1[7, 1], "(", tmp2[7, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[3, 2], "(", tmp2[3, 2],")&", tmp1[5, 2], "(", tmp2[5, 2],")&",
                          tmp1[7, 2], "(", tmp2[7, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",
                          tmp1[3, 3], "(", tmp2[3, 3],")&", tmp1[5, 3], "(", tmp2[5, 3],")&",
                          tmp1[7, 3], "(", tmp2[7, 3],")\\",sep=""))

resu


write.table(resu, file="simu.txt", col.names=F, row.names=F)

############################# run times #######################

resu = NULL


i=1
tmp1 = cbind(est200[9:12, i], est300[9:12, i], est400[9:12, i])
tmp2 = cbind(sd200[9:12, i], sd300[9:12, i], sd400[9:12, i])

resu = rbind(resu, paste(0.05, "& SO(3) &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[2, 1], "(", tmp2[2, 1],")&", tmp1[3, 1], "(", tmp2[3, 1],")&",
                         tmp1[4, 1], "(", tmp2[4, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[2, 2], "(", tmp2[2, 2],")&", tmp1[3, 2], "(", tmp2[3, 2],")&",
                          tmp1[4, 2], "(", tmp2[4, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",tmp1[2, 3], "(", tmp2[2, 3],")&", tmp1[3, 3], "(", tmp2[3, 3],")&",tmp1[4, 3], "(", tmp2[4, 3],")\\",sep=""))



i=4
tmp1 = cbind(est200[9:12, i], est300[9:12, i], est400[9:12, i])
tmp2 = cbind(sd200[9:12, i], sd300[9:12, i], sd400[9:12, i])

resu = rbind(resu, paste( "& KB &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[2, 1], "(", tmp2[2, 1],")&", tmp1[3, 1], "(", tmp2[3, 1],")&",
                         tmp1[4, 1], "(", tmp2[4, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[2, 2], "(", tmp2[2, 2],")&", tmp1[3, 2], "(", tmp2[3, 2],")&",
                          tmp1[4, 2], "(", tmp2[4, 2],")\\",sep=""))


resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",tmp1[2, 3], "(", tmp2[2, 3],")&", tmp1[3, 3], "(", tmp2[3, 3],")&",tmp1[4, 3], "(", tmp2[4, 3],")\\",sep=""))



i=2
tmp1 = cbind(est200[9:12, i], est300[9:12, i], est400[9:12, i])
tmp2 = cbind(sd200[9:12, i], sd300[9:12, i], sd400[9:12, i])


resu = rbind(resu, paste(0.10, "& SO(3) &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                          tmp1[2, 1], "(", tmp2[2, 1],")&", tmp1[3, 1], "(", tmp2[3, 1],")&",
                          tmp1[4, 1], "(", tmp2[4, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[2, 2], "(", tmp2[2, 2],")&", tmp1[3, 2], "(", tmp2[3, 2],")&",
                          tmp1[4, 2], "(", tmp2[4, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",tmp1[2, 3], "(", tmp2[2, 3],")&", tmp1[3, 3], "(", tmp2[3, 3],")&",tmp1[4, 3], "(", tmp2[4, 3],")\\",sep=""))

i=5
tmp1 = cbind(est200[9:12, i], est300[9:12, i], est400[9:12, i])
tmp2 = cbind(sd200[9:12, i], sd300[9:12, i], sd400[9:12, i])



resu = rbind(resu, paste("& KB &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[2, 1], "(", tmp2[2, 1],")&", tmp1[3, 1], "(", tmp2[3, 1],")&",
                         tmp1[4, 1], "(", tmp2[4, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[2, 2], "(", tmp2[2, 2],")&", tmp1[3, 2], "(", tmp2[3, 2],")&",
                          tmp1[4, 2], "(", tmp2[4, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",tmp1[2, 3], "(", tmp2[2, 3],")&", tmp1[3, 3], "(", tmp2[3, 3],")&",tmp1[4, 3], "(", tmp2[4, 3],")\\",sep=""))

i=3
tmp1 = cbind(est200[9:12, i], est300[9:12, i], est400[9:12, i])
tmp2 = cbind(sd200[9:12, i], sd300[9:12, i], sd400[9:12, i])


resu = rbind(resu, paste(0.15, "& SO(3) &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[2, 1], "(", tmp2[2, 1],")&", tmp1[3, 1], "(", tmp2[3, 1],")&",
                         tmp1[4, 1], "(", tmp2[4, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[2, 2], "(", tmp2[2, 2],")&", tmp1[3, 2], "(", tmp2[3, 2],")&",
                          tmp1[4, 2], "(", tmp2[4, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",tmp1[2, 3], "(", tmp2[2, 3],")&", tmp1[3, 3], "(", tmp2[3, 3],")&",tmp1[4, 3], "(", tmp2[4, 3],")\\",sep=""))


i=6
tmp1 = cbind(est200[9:12, i], est300[9:12, i], est400[9:12, i])
tmp2 = cbind(sd200[9:12, i], sd300[9:12, i], sd400[9:12, i])


resu = rbind(resu, paste("& KB &", 200, "&", tmp1[1, 1], "(", tmp2[1, 1],")&",
                         tmp1[2, 1], "(", tmp2[2, 1],")&", tmp1[3, 1], "(", tmp2[3, 1],")&",
                         tmp1[4, 1], "(", tmp2[4, 1],")\\",sep=""))

resu = rbind(resu, paste( "& &", 300, "&", tmp1[1, 2], "(", tmp2[1, 2],")&",
                          tmp1[2, 2], "(", tmp2[2, 2],")&", tmp1[3, 2], "(", tmp2[3, 2],")&",
                          tmp1[4, 2], "(", tmp2[4, 2],")\\",sep=""))

resu = rbind(resu, paste( "& &", 400, "&", tmp1[1, 3], "(", tmp2[1, 3],")&",tmp1[2, 3], "(", tmp2[2, 3],")&", tmp1[3, 3], "(", tmp2[3, 3],")&",tmp1[4, 3], "(", tmp2[4, 3],")\\",sep=""))

resu


write.table(resu, file="simu-time.txt", col.names=F, row.names=F)



############################# For comparison with Lu or Kong ###############

#################### sigma = 0.05, 0.10, 0.15 ################
set.seed(2022)
registerDoMC(20)

result100 = foreach(i.rep = 1:200, .combine = cbind) %dopar%
  { 
    
    sim_m1 = simulation_M(n.fit = 100, n.test = 100, h=0.1, scale=c(1/4), sigma=0.05)
    sim_m2 = simulation_M(n.fit = 100, n.test = 100, h=0.1, scale=c(1/4), sigma=0.10)
    sim_m3 = simulation_M(n.fit = 100, n.test = 100, h=0.1, scale=c(1/4), sigma=0.15)
    
    sim_k1 = simulation_K(n.fit = 100, n.test = 100, h=1, scale=c(1/4), sigma=0.05)
    sim_k2 = simulation_K(n.fit = 100, n.test = 100, h=1, scale=c(1/4), sigma=0.10)
    sim_k3 = simulation_K(n.fit = 100, n.test = 100, h=1, scale=c(1/4), sigma=0.15)
    
    rbind(sim_m1, sim_m2, sim_m3, sim_k1, sim_k2, sim_k3)
  }

save(result100,file="simu_runtime100.Rdata")

result = result100
round(apply(result, 1, function(x){mean(x)}), 4)

load("simu_runtime100.Rdata")

result = result100
est100 = matrix(round(apply(result, 1, function(x){mean(x,na.rm = TRUE)}), 4),12,)
sd100 = matrix(round(apply(result, 1, function(x){sd(x,na.rm = TRUE)}), 4),12,)

est100 = est100[c(1, 3, 5, 7), ]
dimnames(est100)[[2]] <- c("SO(3)", "Klein")
dimnames(est100)[[1]] <- c("Our(1/4)", "Our(all)", "Gaussian", "Laplacian")


sd100 = sd100[c(1, 3, 5, 7), ]
dimnames(sd100)[[2]] <- c("SO(3)", "Klein")
dimnames(sd100)[[1]] <- c("Our(1/4)", "Our(all)", "Gaussian", "Laplacian")

est100



resu = NULL

tmp1 = est100[ , 1]
tmp2 = sd100[ , 1]

resu = rbind(resu, paste(0.05, "& SO(3) &", 100, "&", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")&",
                         tmp1[4], "(", tmp2[4],")\\",sep=""))

tmp1 = est100[ , 4]
tmp2 = sd100[ , 4]

resu = rbind(resu, paste(  "   & Klein &",  "   &", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")&",
                         tmp1[4], "(", tmp2[4],")\\",sep=""))


tmp1 = est100[ , 2]
tmp2 = sd100[ , 2]

resu = rbind(resu, paste(0.10, "& SO(3) &", 100, "&", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")&",
                         tmp1[4], "(", tmp2[4],")\\",sep=""))

tmp1 = est100[ , 5]
tmp2 = sd100[ , 5]

resu = rbind(resu, paste(  "   & Klein &",  "   &", tmp1[1], "(", tmp2[1],")&",
                           tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")&",
                           tmp1[4], "(", tmp2[4],")\\",sep=""))


tmp1 = est100[ , 3]
tmp2 = sd100[ , 3]

resu = rbind(resu, paste(0.15, "& SO(3) &", 100, "&", tmp1[1], "(", tmp2[1],")&",
                         tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")&",
                         tmp1[4], "(", tmp2[4],")\\",sep=""))

tmp1 = est100[ , 6]
tmp2 = sd100[ , 6]

resu = rbind(resu, paste(  "   & Klein &",  "   &", tmp1[1], "(", tmp2[1],")&",
                           tmp1[2], "(", tmp2[2],")&", tmp1[3], "(", tmp2[3],")&",
                           tmp1[4], "(", tmp2[4],")\\",sep=""))


resu


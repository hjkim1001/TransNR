path <- "C:/Users/USER/Desktop/Transferlearning/Function/"
source(paste0(path,"Functions_rankestimation_simul.R")) # load r function



##### Simulation model in Section 3.4 with n0=120, p=80, q=50, r=8, h=10

### Simulating 200 target data
# rep: the number of replicates
# n=n0, r=rank of B
# corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
# cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
targestset <- generator_target_normal(rep=100,n=120,p=80,q=50,r=8, corx=0.5,cory=0)



### Simulating 200 source datasets
## rank source= rank of contrast
# corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
# cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
# generating 200 replicates for a source data with nk=300, transferring level=2, rank_source=5
sources <-generator_source_rankestimationsimul(rep=200, n=150,p=80,q=50,B=targestset$B,
                                                                              corx=0.5,cory=0,
                                                                              h=2,numsource=1, ranksource =5)

### Fitting algorithms to the first five replicates

## Ah-Trans-NR
# Ylist: contains the simulated Y
# Xlist: contains the simulated X
# auxYlist_list: contains the simulated auxYlist, wheer auxYlist is the list containing source response matrices
# auxXlist_list:  contains the simulated auxXlist, wheer auxYlist is the list containing source design matrices
# repstart, repend: simulation results are obtaeind from repstart to repend, 
#  e.g., when two values are set as 1 and 100, simulation results are obtained from the first 100 simulated data

# lamseq_w: candiates of tuning parameters for the firs step
# lamseq_delta: candiates of tuning parameters for the bias-correction step

# maxiter_inital: the maxium of the number of iterations in the ADMM algorithm for computing initial estimate
# tol_inital: tolerance error in the ADMM algorithm for computing initial estimate
# maxiter_biascorrection: the maxium of the number of iterations in the ADMM algorithm for computing contrast
# tol_inital: tolerance error in the ADMM algorithm for computing contrast
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

lamseq_w <-  c(0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1, 1.2, 1.5, 2, 2.5, 3, 5, 7, 10, 12, 15,18,20,25, 30, 40, 50, 60, 70,100)
lamseq_delta <- c(lamseq_w,  100, 200, 500, 1000, 2000, 5000, 7000, 10000, 15000, 20000)
Repfit_TransNR_BIC(Ylist=targestset$Ylist,Xlist=targestset$Xlist,
                   auxYlist_list = sources$auxYlist_list, auxXlist_list=sources$auxXlist_list,
                   lamseq_w=lamseq_w, lamseq_delta=lamseq_delta, numsim=5,
                   maxiter_inital=300, maxiter_biascorrection=300,
                   tol_inital=1e-04, tol_biascorrection=1e-04) -> res_transNR

Repfit_TransSCAD_BIC(Ylist=targestset$Ylist, Xlist=targestset$Xlist,
                     auxYlist_list=sources$auxYlist_list, auxXlist_list=sources$auxXlist_list, 
                     lamseq_w=lamseq_w, lamseq_delta=lamseq_delta, 
                     numsim=5, maxiter_inital = 300, maxiter_biascorrection = 300, 
                     tol_inital = 1e-04, tol_biascorrection = 1e-04) -> res_transSCAD

##### Results 
### PCR 
## Trans- SCAD, Trans-NR

rankvec_transSCAD <- unlist(lapply(res_transSCAD, function(x){ x$optrank_BIC }))
rankvec_transNR <- unlist(lapply(res_transNR, function(x){ x$optrank_BIC }))
c( sum(rankvec_transSCAD==8), sum(rankvec_transNR==8 ) )/5 
### rank 
# first row: average of ranks for Trans-SCAD, standard deviation of ranks for Trans-SCAD 
# second row: average of ranks for Trans-NR, standard deviation of ranks for Trans-NR 

round( rbind(c( mean(rankvec_transSCAD), sd(rankvec_transSCAD)),
             c(  mean(rankvec_transNR), sd(rankvec_transNR)) ), 2)

### Estimation error  
# first row: average of EEs for Trans-SCAD, standard deviation of EEs for Trans-SCAD 
# second row: average of EEs for Trans-NR, standard deviation of EEs for Trans-NR 

Errorvec_transSCAD <- unlist(lapply(res_transSCAD, function(x){ sqrt( mean( (x$B[-1,]- targestset$B)^2 )) }))
Errorvec_transNR <- unlist(lapply(res_transNR, function(x){ sqrt( mean( (x$B[-1,]- targestset$B)^2 )) }))

round( rbind(c( mean(Errorvec_transSCAD), sd(Errorvec_transSCAD)),
              c(  mean(Errorvec_transNR), sd(Errorvec_transNR)) ), 3)

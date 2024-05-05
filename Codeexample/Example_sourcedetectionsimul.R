path <- "C:/Users/USER/Desktop/Transferlearning/Function/"
source(paste0(path,"Functions_naiveapproaches.R")) # load r function
source(paste0(path,"Functions_MSDtrans.R")) # load r function
source(paste0(path,"Functions_FSDtrans.R")) # load r function
source(paste0(path,"Functions_nr-basedmethods_simul.R")) # load r function

###### Simulation setup: Simulation model in Section 3.1 with n0=50, p=30, q=20, r= 5, h=5, 
# nk=50 for k=1,...,K(=8).


### Simulating 100 target data
# rep: the number of replicates
# n=n0, r=rank of B
# corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
# cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
targestset <- generator_target_simulsourcedetection(rep=100,n=50,p=30,q=20,r=5, corx=0.5,cory=0)

### Simulating 100 source datasets
## rank source= rank of contrast
# generating 100 replicates for two source data with nk=50, transferring level=30, rank_source=5
Informativesets_numsourcesfour_h30 <- generator_source_simulsourcedetection(rep=100, n=50,p=30,q=20,
                                                                            rank_source=5,B=targestset$B,
                                                                            corx=0.5,cory=0,
                                                                            h=30,numsource=2)

# generating 100 repicates for three source data with nk=50, transferring level=20, rank_source=5
Informativesets_numsourcesfour_h20 <- generator_source_simulsourcedetection(rep=100, n=50,p=30,q=20,
                                                                            rank_source=5,B=targestset$B,
                                                                            corx=0.5,cory=0,
                                                                            h=20,numsource=3)

# generating 100 replicates for three source data with nk=50, transferring level=5, rank_source=5
Informativesets_numsourcesfour_h5 <-generator_source_simulsourcedetection(rep=100, n=50,p=30,q=20,
                                                                          rank_source=5,B=targestset$B,
                                                                          corx=0.5,cory=0,
                                                                          h=5,numsource=3)

# Sourceset_X[[i]]: a list containig eight source covariate matrices
# Sourceset_Y[[i]]: a list containig eight source response matrices

Sourceset_X <- Sourceset_Y <- list()
for(ss in 1:100){
  Sourceset_X_ss <- list()
  Sourceset_Y_ss <- list()
  ############## informative sources #############
  for(aa in 1:3){
    Sourceset_X_ss[[aa]] <- Informativesets_numsourcesfour_h5$auxXlist_list[[ss]][[aa]]
    Sourceset_Y_ss[[aa]] <- Informativesets_numsourcesfour_h5$auxYlist_list[[ss]][[aa]]
  }
  ###################### h=20 ###############
  for(aa in 4:6){
    Sourceset_X_ss[[aa]] <- Informativesets_numsourcesfour_h20$auxXlist_list[[ss]][[aa-3]]
    Sourceset_Y_ss[[aa]] <- Informativesets_numsourcesfour_h20$auxYlist_list[[ss]][[aa-3]]
  }
  ############### h=30 ######################
  Sourceset_X_ss[[7]] <- Informativesets_numsourcesfour_h30$auxXlist_list[[ss]][[1]]
  Sourceset_Y_ss[[7]] <- Informativesets_numsourcesfour_h30$auxYlist_list[[ss]][[1]]
  
  Sourceset_X_ss[[8]] <- Informativesets_numsourcesfour_h30$auxXlist_list[[ss]][[2]]
  Sourceset_Y_ss[[8]] <- Informativesets_numsourcesfour_h30$auxYlist_list[[ss]][[2]]
  
  Sourceset_X[[ss]] <- Sourceset_X_ss
  Sourceset_Y[[ss]] <- Sourceset_Y_ss
  
}



### Obtaining results for the five methods from the first simulated data to the third simulated data
## five methods: NR, Pooled-NR, [K]-Trans, MSD-Trans, FSD-Trans
# Ylist: contains the simulated Y
# Xlist: contains the simulated X
# auxYlist_list: contains the simulated auxYlist, wheer auxYlist is the list containing source response matrices
# auxXlist_list:  contains the simulated auxXlist, wheer auxYlist is the list containing source design matrices
# repstart, repend: simulation results are obtaeind from repstart to repend, 
#  e.g., when two values are set as 1 and 100, simulation results are obtained from the first 100 simulated data


testrep1to3 <- rep_forsimul(Ylist=targestset$Ylist,
                            Xlist=targestset$Xlist,
                            auxYlist_list=Sourceset_Y,
                            auxXlist_list=Sourceset_X, repstart=1, repend=3, Btrue=targestset$B)

# Summary of source detection results
Summary_ft_detection(testrep1to3, numinf=3) #numinf: the number of informative sources

# Summary of estimation results
Summary_ft_est(testrep1to3)

########################### Hetergoeneous design setting ###########################################
# ranks of contrasts, indexes of informative sources are randomly generated.
tgarget_n050p30r5q20 <- generator_target_simulsourcedetection(rep=100, n=50, p=30, q=20, r= 5, corx=0.5,cory=0)
Sourceset_h5 <- generator_source_simulsourcedetection_het(rep=100, nvec=rep(50,5),hvec = c(25, 5, 5, 25, 5), 
                                                           p=30,q=20,ranksourcevec=c(11,6,9,3,4),
                                                           B=tgarget_n050p30r5q20$B,cory=0,numsource=5)
testrep1to3_het <- rep_forsimul(Ylist=tgarget_n050p30r5q20$Ylist,
                            Xlist=tgarget_n050p30r5q20$Xlist,
                            auxYlist_list=Sourceset_h5$auxYlist_list,
                            auxXlist_list=Sourceset_h5$auxXlist_list, repstart=1, repend=3, 
                            Btrue=tgarget_n050p30r5q20$B)

# Summary of source detection results
## truesupport: informative set
Summary_ft_detection_het(testrep1to3_het, truesupport=c(2,3,5))
# Summary of estimation results
Summary_ft_est(testrep1to3_het)

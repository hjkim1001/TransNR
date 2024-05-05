path <- "C:/Users/USER/Desktop/Transferlearning/Function/"
source(paste0(path,"Functions_naiveapproaches.R")) # load r function
source(paste0(path,"Functions_MSDtrans.R")) # load r function
source(paste0(path,"Functions_FSDtrans.R")) # load r function
source(paste0(path,"Functions_transSCAD.R")) # load r function


# Load data
path <- "C:/Users/USER/Desktop/Transferlearning/data"
dataname <- "CCLEdataset.RDATA"
load(paste(path,dataname,sep="/"))

path <- "C:/Users/USER/Desktop/Transferlearning/data"
dataname <- "CCLEdataset_trainingtest.RDATA"
load(paste(path,dataname,sep="/"))

########### Common parameters ###################################################################
# Y: response matrix in the target data
# X: matrix of covariates in the target data
# B: A in algorithm in Section S2.1 (if it is set as NULL, null matrix will be used as initial value)
# L: C in algorithm in Section S2.1 (if it is set as NULL, null matrix will be used as initial value)
# eta: rho in algorithm in Section S2.1 (recommend =1)
# tol: tolerance error (default=1e-04)
# maxiter: the number of maximum of iterations (recommend maxiter=100)

###########################################################################################

########### Common returned values ############################
# B: Bhat
# What: What
# Deltahat: Deltahat
# lamopt: lambde selected by CV in NR and Pooled-NR
# lamoptw: lambda_W selected by [K]-Trans, MSD-Trans, FSD-Trans
# lamoptdelta: lambda_delta selected by [K]-Trans, MSD-Trans, FSD-Trans
# detectedsources: detected sources in FSD and MSD
################################################################

### Five-fold CV for NR
lamLassoseq <- c(0.001, 0.01, 0.05, 0.1,0.15,0.2, 0.3,0.35, 0.4, 0.5, 0.6, 0.7,  0.8, 0.9, 1, 1,1, 1.2, 1.5, 1.8, 2, 2.2 ,2.5)
# lamseq: candidates of lambda in algorithm in Section S2.1
cv.nuclear(Y=CCLEdataset_trainingtest[[1]]$Ytraining, 
           X=CCLEdataset_trainingtest[[1]]$Xtraining, 
           B=NULL,
           L=NULL,
           eta=1,lamseq = lamLassoseq,
           tol=1e-04,maxiter=100) -> cvLasso 

p <- ncol(CCLEdataset_trainingtest[[1]]$Xtraining)
q <- ncol(CCLEdataset_trainingtest[[1]]$Ytrainin)
dLasso <- svd(matrix(cvLasso$B[-1,],p,q))$d # corresponding singular value
predLasso <- cbind(1,CCLEdataset_trainingtest[[1]]$Xtest) %*% cvLasso$B
Lassoerr <- mean( (predLasso-CCLEdataset_trainingtest[[1]]$Ytest)^2 ,na.rm=TRUE) # MSPE

### Five-fold CV for Pooled-NR
lamseq_w <- lamLassoseq
# auxYlist: a list contaring resposne matrix for source datasets
# auxXlist: a list contaring coviaraites matrix for source datasets
# lamseq: candidates of lambda in algorithm in Section S2.1
cv.pooledNR(Y=CCLEdataset_trainingtest[[1]]$Ytraining, 
            X=CCLEdataset_trainingtest[[1]]$Xtraining, 
            auxYlist = CCLEdataset$auxYlist,
            auxXlist = CCLEdataset$auxXlist,  
            B=NULL,
            L=NULL,eta=1,lamseq = lamLassoseq,
            tol=1e-04,maxiter=100) -> PooledNRres

dPooledNR <- svd(matrix(PooledNRres$B[-1,],p,q))$d # corresponding singular value
predPooledNR <- cbind(1,CCLEdataset_trainingtest[[1]]$Xtest) %*% PooledNRres$B
PooledNRerr <- mean( (predPooledNR-CCLEdataset_trainingtest[[1]]$Ytest)^2 ,na.rm=TRUE) # MSPE

### Five-fold CV for [K]-Trans
lam_delta <- c(lamseq_w,3)
# auxYlist: a list contaring resposne matrix for source datasets
# auxXlist: a list contaring coviaraites matrix for source datasets
# lamseq_w: candidates of lambda_W in the step 1 of Algorithm 1 in the main paper
# lamseq_delta: candidates of lambda_W in the step 1 of Algorithm 1 in the main paper
cv.twostep(Y=CCLEdataset_trainingtest[[1]]$Ytraining, 
           X=CCLEdataset_trainingtest[[1]]$Xtraining, 
           B=NULL, L=NULL,eta=1,
           lamseq_w = lamseq_w,
           lamseq_delta=lam_delta,
           tol=1e-04,maxiter=100,
           auxYlist=CCLEdataset$auxYlist, 
           auxXlist=CCLEdataset$auxXlist) -> TwosteptransRes

dtwostep <- svd(matrix(TwosteptransRes$B[-1,],p,q))$d # corresponding singular value
predTwostep <- cbind(1,CCLEdataset_trainingtest[[1]]$Xtest) %*% TwosteptransRes$B
Twosteperr <- mean( (predTwostep-CCLEdataset_trainingtest[[1]]$Ytest)^2 ,na.rm=TRUE) # MSPE



### Running MSD-Trans 
# auxYlist: a list contaring resposne matrix for source datasets
# auxXlist: a list contaring coviaraites matrix for source datasets
# lamseq_w: candidates of lambda_W in the step 1 of Algorithm 1 in the main paper
# lamseq_delta: candidates of lambda_W in the step 1 of Algorithm 1 in the main paper
# nfold: the number of folds to tune penalty parameters
# nfold_choiceforC: the number of folds to select
# nfold_selectionstep_pe: the number of folds to estimate predition error when performing source detection
# C: the candidates of C

cl <- makeCluster(3)
registerDoParallel(cl) # for fsd, msd
MSDtrans(Y=CCLEdataset_trainingtest[[1]]$Ytraining, 
         X=CCLEdataset_trainingtest[[1]]$Xtraining, 
         B=NULL, L=NULL,eta=1,
         lamseq_w = lamseq_w,
         lamseq_delta=lam_delta,
         tol=1e-04,maxiter=100,
         auxYlist=CCLEdataset$auxYlist, 
         auxXlist=CCLEdataset$auxXlist,
         nfold=5,nfold_choiceforC=3,
         nfold_selectionstep_pe=3,
         C=c(0.001, 0.01,  0.05,  0.1)) -> MSDRes

dMSD <- svd(matrix(MSDRes$B[-1,],p,q))$d # corresponding singular value
predMSD<- cbind(1,CCLEdataset_trainingtest[[1]]$Xtest) %*% MSDRes$B
MSDerr <- mean( (predTwostep-CCLEdataset_trainingtest[[1]]$Ytest)^2 ,na.rm=TRUE) # MSPE


#### Running FSD-Trans-NR
FSDtrans(Y=CCLEdataset_trainingtest[[1]]$Ytraining, 
         X=CCLEdataset_trainingtest[[1]]$Xtraining, 
         B=NULL, L=NULL,eta=1,
         lamseq_w = lamseq_w,
         lamseq_delta=lam_delta,
         tol=1e-04,maxiter=100,
         auxYlist=CCLEdataset$auxYlist, 
         auxXlist=CCLEdataset$auxXlist,
         nfold=5,nfold_choiceforC=3,
         nfold_selectionstep_pe=3,
         C=c(0.001, 0.01,  0.05,  0.1)) -> FSDRes

dFSD <- svd(matrix(FSDRes$B[-1,],p,q))$d # corresponding singular value
predFSD <- cbind(1,CCLEdataset_trainingtest[[1]]$Xtest) %*% FSDRes$B
FSDerr <- mean( (predFSD-CCLEdataset_trainingtest[[1]]$Ytest)^2 ,na.rm=TRUE) # MSPE
stopCluster(cl)
stopImplicitCluster()

### Running FSD-Trans-SCAD
lamseq_w_SCAD <-  c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1, 1.2, 1.5, 2, 2.5)
lamseq_delta_SCAD <- c(lamseq_w_SCAD,3)
sourceindex <- FSDRes$detectedsources
auxYlist_SCAD <- auxXlsit_SCAD <- list()
for(ss in 1:length(sourceindex)){
  auxYlist_SCAD[[ss]] <- CCLEdataset$auxYlist[[sourceindex[ss]]]
  auxXlsit_SCAD[[ss]] <- CCLEdataset$auxXlist[[sourceindex[ss]]]
}
FSDSCADres <- BIC.TransSCADinitcv(Y=CCLEdataset_trainingtest[[1]]$Ytraining, 
                            X=CCLEdataset_trainingtest[[1]]$Xtraining,
                            auxYlist=auxYlist_SCAD,
                            auxXlist=auxXlsit_SCAD,
                            lamseq_w=lamseq_w_SCAD, 
                            lamseq_delta=lamseq_delta_SCAD,
                            eta=1, a=3.7,
                            B=NULL,L=NULL,Delta=NULL,H=NULL,Pi=NULL,
                            maxiter_inital=300, maxiter_biascorrection=300,
                            tol_inital=1e-04, tol_biascorrection=1e-04,
                            standardize=T)

dFSDSCAD <- svd(matrix(FSDSCADres$B[-1,],p,q))$d # corresponding singular value
predFSDSCAD <- cbind(1,CCLEdataset_trainingtest[[1]]$Xtest) %*% FSDSCADres$B
FSDSCADerr <- mean( (predFSDSCAD-CCLEdataset_trainingtest[[1]]$Ytest)^2 ,na.rm=TRUE) # MSPE


### Results 
## prediction error
round(c(Lassoerr, PooledNRerr, Twosteperr, MSDerr, FSDerr, FSDSCADerr),3)
computerank <- function(x,tol=1e-03){
  sum(x>1e-02)
}
c(computerank(dLasso), computerank(dPooledNR),
  computerank(dtwostep), computerank(dMSD), computerank(dFSD),
  computerank(dFSDSCAD))

##### Multiple run
nrep <- 3
errorlist <- list() # to save error
ranklist <- list() # to save rank
Blist <- list() # to save estimates
lamLassoseq <- c(0.001, 0.01, 0.05, 0.1,0.15,0.2, 0.3,0.35, 0.4, 0.5, 0.6, 0.7,  0.8, 0.9, 1, 1,1, 1.2, 1.5, 1.8, 2, 2.2 ,2.5)
lamseq_w <- lamLassoseq
lamseq_delta <- c(lamLassoseq,3)
p <- ncol(CCLEdataset_trainingtest[[1]]$Xtraining)
q <- ncol(CCLEdataset_trainingtest[[1]]$Ytrainin)
cl <- makeCluster(3)
registerDoParallel(cl) # for fsd, msd

lamseq_w_SCAD <-  c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1, 1.2, 1.5, 2, 2.5)
lamseq_delta_SCAD <- c(lamseq_w_SCAD,3)
for(i in 1:nrep){
  cv.nuclear(Y=CCLEdataset_trainingtest[[i]]$Ytraining, 
             X=CCLEdataset_trainingtest[[i]]$Xtraining, 
             B=NULL,
             L=NULL,
             eta=1,lamseq = lamLassoseq,
             tol=1e-04,maxiter=100) -> cvLasso 
  
  cv.pooledNR(Y=CCLEdataset_trainingtest[[i]]$Ytraining, 
              X=CCLEdataset_trainingtest[[i]]$Xtraining, 
              auxYlist = CCLEdataset$auxYlist,
              auxXlist = CCLEdataset$auxXlist,  
              B=NULL,
              L=NULL,eta=1,lamseq = lamLassoseq,
              tol=1e-04,maxiter=100) -> PooledNRres
  
  cv.twostep(Y=CCLEdataset_trainingtest[[i]]$Ytraining, 
             X=CCLEdataset_trainingtest[[i]]$Xtraining, 
             B=NULL, L=NULL,eta=1,
             lamseq_w = lamseq_w,
             lamseq_delta=lam_delta,
             tol=1e-04,maxiter=100,
             auxYlist=CCLEdataset$auxYlist, 
             auxXlist=CCLEdataset$auxXlist) -> TwosteptransRes
  
  MSDtrans(Y=CCLEdataset_trainingtest[[i]]$Ytraining, 
           X=CCLEdataset_trainingtest[[i]]$Xtraining, 
           B=NULL, L=NULL,eta=1,
           lamseq_w = lamseq_w,
           lamseq_delta=lam_delta,
           tol=1e-04,maxiter=100,
           auxYlist=CCLEdataset$auxYlist, 
           auxXlist=CCLEdataset$auxXlist,
           nfold=5,nfold_choiceforC=3,
           nfold_selectionstep_pe=3,
           C=c(0.001, 0.01,  0.05,  0.1)) -> MSDRes
  
  FSDtrans(Y=CCLEdataset_trainingtest[[i]]$Ytraining, 
           X=CCLEdataset_trainingtest[[i]]$Xtraining, 
           B=NULL, L=NULL,eta=1,
           lamseq_w = lamseq_w,
           lamseq_delta=lam_delta,
           tol=1e-04,maxiter=100,
           auxYlist=CCLEdataset$auxYlist, 
           auxXlist=CCLEdataset$auxXlist,
           nfold=5,nfold_choiceforC=3,
           nfold_selectionstep_pe=3,
           C=c(0.001, 0.01,  0.05,  0.1)) -> FSDRes
  
  
  predLasso <- cbind(1,CCLEdataset_trainingtest[[i]]$Xtest) %*% cvLasso$B
  Lassoerr <- mean( (predLasso-CCLEdataset_trainingtest[[i]]$Ytest)^2 ,na.rm=TRUE) 
  
  predPooledNR <- cbind(1,CCLEdataset_trainingtest[[i]]$Xtest) %*% PooledNRres$B
  PooledNRerr <- mean( (predPooledNR-CCLEdataset_trainingtest[[i]]$Ytest)^2 ,na.rm=TRUE) 
  
  predTwostep <- cbind(1,CCLEdataset_trainingtest[[i]]$Xtest) %*% TwosteptransRes$B
  Twosteperr <- mean( (predTwostep-CCLEdataset_trainingtest[[i]]$Ytest)^2 ,na.rm=TRUE) 
  
  predMSD<- cbind(1,CCLEdataset_trainingtest[[i]]$Xtest) %*% MSDRes$B
  MSDerr <- mean( (predMSD-CCLEdataset_trainingtest[[i]]$Ytest)^2 ,na.rm=TRUE) 
  
  predFSD <- cbind(1,CCLEdataset_trainingtest[[i]]$Xtest) %*% FSDRes$B
  FSDerr <- mean( (predFSD-CCLEdataset_trainingtest[[i]]$Ytest)^2 ,na.rm=TRUE) 
  
  
  
  sourceindex <- FSDRes$detectedsources
  if(length(sourceindex)==0){
    FSDSCADerr <- Lassoerr
    Bscad <- cvLasso$B
    rankSCAD <- sum(svd(matrix(Bscad[-1,],p,q))$d > 1e-02)
  }else{
    auxYlist_SCAD <- auxXlist_SCAD <- list()
    for(ss in 1:length(sourceindex)){
      auxYlist_SCAD[[ss]] <- CCLEdataset$auxYlist[[sourceindex[ss]]]
      auxXlist_SCAD[[ss]] <- CCLEdataset$auxXlist[[sourceindex[ss]]]
    }
    FSDSCADres <-BIC.TransSCAD(Y=CCLEdataset_trainingtest[[i]]$Ytraining, 
                                X=CCLEdataset_trainingtest[[i]]$Xtraining,
                                auxYlist=auxYlist_SCAD,
                                auxXlist=auxXlist_SCAD,
                                lamseq_w=lamseq_w_SCAD, 
                                lamseq_delta=lamseq_delta_SCAD,
                                eta=1, a=3.7,
                                B=NULL,L=NULL,Delta=NULL,H=NULL,Pi=NULL,
                                maxiter_inital=300, maxiter_biascorrection=300,
                                tol_inital=1e-04, tol_biascorrection=1e-04,
                                standardize=T)
    predFSDSCAD <- cbind(1,CCLEdataset_trainingtest[[i]]$Xtest) %*% FSDSCADres$B
    FSDSCADerr <- mean( (predFSDSCAD-CCLEdataset_trainingtest[[i]]$Ytest)^2 ,na.rm=TRUE) # MSPE
    Bscad <- FSDSCADres$B
  }
  
  
  error_i <- c(Lassoerr, PooledNRerr, Twosteperr, MSDerr, FSDerr, FSDSCADerr)
  names(error_i) <- c("NR", "Pooled-NR", "[K]-Trans", "MSD-Trans", "FSD-Trans", "FSD-Trans-SCAD")
  errorlist[[i]] <- error_i # save test errors
  
  Blist_i <- list(cvLasso$B, PooledNRres$B, TwosteptransRes$B,
                  MSDRes$B, FSDRes$B, Bscad)
  
  names(Blist_i) <- c("NR", "Pooled-NR", "[K]-Trans", "MSD-Trans", "FSD-Trans", "FSD-Trans-SCAD")
  Blist[[i]]  <- Blist_i # save Bhat 
  
  ranklist_i <- unlist(lapply(Blist_i, function(x){
    sum(svd(matrix(x[-1,],p,q))$d > 1e-03)
  }))
  ranklist_i[2] <- sum(svd(Blist_i[[2]])$d > 0.005)
  names(ranklist_i) <- c("NR", "Pooled-NR", "[K]-Trans", "MSD-Trans", "FSD-Trans", "FSD-Trans-SCAD")
  ranklist[[i]] <- ranklist_i
  
  
  
}

#### mean prediction error
apply(do.call('rbind',errorlist),2,mean)
#### standard deviation of prediction errors
apply(do.call('rbind',errorlist),2,sd)

#### mean ranks
apply(do.call('rbind',ranklist),2,mean)
#### standard deviation of ranks
apply(do.call('rbind',ranklist),2,sd)

#### Wilcoxon singed rank tests
FSDvsNR <- wilcox.test(do.call('rbind',errorlist)[,5],
                       do.call('rbind',errorlist)[,1], paired = TRUE, alternative = "two.sided")
FSDvsNR$p.value

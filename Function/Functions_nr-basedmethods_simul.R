rep_forsimul <- function(Ylist, Xlist, auxYlist_list, auxXlist_list, repstart, repend, Btrue){
  # Ylist: contains the simulated Y
  # Xlist: contains the simulated X
  # auxYlist_list: contains the simulated auxYlist, wheer auxYlist is the list containing source response matrices
  # auxXlist_list:  contains the simulated auxXlist, wheer auxYlist is the list containing source design matrices
  # repstart, repend: simulation results are obtaeind from repstart to repend, 
  #  e.g., when two values are set as 1 and 100, simulation results are obtained from the first 100 simulated data
  
  # Btrue: true coefficient matrix for a simulation model
  
  library(parallel)
  library(doParallel)
  library(foreach)
  cl <- makeCluster(3)
  registerDoParallel(cl)
  
  library(splitTools)
  
  lamLassoseq <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.8, 0.9, 1, 1.2, 1.5, 2, 2.5, 3, 5, 7, 10, 12, 15) # lambds for NR
  lam_wNaiveseq <- lamLassoseq # lamdbds for the first step
  lam_deltaNaiveseq <- c(lamLassoseq,18,20,25) # lambda for debiasing step
  maxit <- 100
  C <- c(0.001, 0.01, 0.05, 0.1) # candidates for C, which is the perameter in source detection
  
  B.update <- function(L,U,eta, covXY, gramridgeInv){
    B_pred <- covXY - U + eta*L
    return( (gramridgeInv %*% B_pred) )
  }
  
  L.update <- function(B,U,eta,lambda,p,q){
    Lpred <- B + U/eta
    svdres <- corpcor::fast.svd(Lpred)
    d <- svdres$d
    dthres <- pmax(d-(lambda/eta),0)
    return( list(L = (svdres$u %*% diag(dthres) %*% t(svdres$v)),d=dthres ))
  }
  
  ADMM.nuclear <- function(Y,X,B,L,eta,lambda,maxiter,tol=1e-04,
                           standardize=TRUE){
    n <- nrow(Y)
    p <- ncol(X)
    q <- ncol(Y)
    
    if(standardize){
      Ymean <- apply(Y,2,mean)
      Ycenter <- Y - matrix(Ymean,n,q,byrow=T)
      Xmean <- apply(X,2,mean)
      Xcenter <- X - matrix(Xmean,n,p,byrow=T)
      Xcenternorm <- apply(Xcenter,2,function(x){sqrt(sum(x^2)/length(x))})
      Xstd <- Xcenter / matrix(Xcenternorm,n,p,byrow=T)
      X <- Xstd
      Y <- Ycenter
    }
    
    
    Onematrix <- matrix(1,p,q)
    gramridge <- diag(eta,p,p) + (1/n)*t(X)%*%X
    gramridgeInv <- chol2inv(chol(gramridge))
    covXY <- (1/n)*t(X)%*%Y
    
    iter <- 1
    isstop <- FALSE
    obj <- NULL
    errB <- errL <- errU <- NULL
    
    if(is.null(B)){
      B <- matrix(0,p,q)
    }
    if(is.null(L)){
      L <- matrix(0,p,q)
    }
    U <- B-L
    
    while(iter <= maxiter & !(isstop)){
      Bold <- B; Lold <- L; Uold <- U
      B <- B.update(L,U,eta, covXY, gramridgeInv)
      Lres <- L.update(B,U,eta,lambda,p,q)
      L <- Lres$L
      U <- U+eta*(B-L)
      res <- Y - X %*% B
      obj[iter] <- (1/(2*n))*sum(res^2) + lambda * sum(Lres$d )
      errB[iter] <-  sum((Bold -B)^2)
      errL[iter] <-  sum((Lold -L)^2)
      errU[iter] <-  sum((Uold -U)^2)
      isstop <- (max(errB[iter],errL[iter],errU[iter]) <= tol)
      iter <- iter + 1
    }
    if(standardize){
      B <- matrix(1/Xcenternorm,p,q) * B # rescale
      L <- matrix(1/Xcenternorm,p,q) * L
      Bint <- matrix(Ymean,1,q) - matrix(Xmean,1,p) %*% B
      Lint <- matrix(Ymean,1,q) - matrix(Xmean,1,p) %*% L
      B <- rbind(Bint,B)
      L <- rbind(Lint,L)
      
    }
    
    return(list(B=B,L=L,obj=obj,errB=errB,errL=errL,errU=errU,iter=iter-1,d=Lres$d))
    
  }
  
  cv.nuclear <- function(Y,X,B,L,eta,lamseq,maxiter=100,tol=1e-04,nfold=5){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    ######################### Divide index ##################################
    set.seed(100)
    n <- nrow(Y)
    ind <- 1:n
    foldset <- list()
    samplesizes <- rep(n %/% nfold, nfold)
    if( (n %% 5) > 0 ){
      samplesizes[1:(n %% 5)] <- samplesizes[1:(n %% 5)] + 1
    }
    
    for( k in 1:(nfold-1) ){
      foldset[[k]] <- sample(ind, samplesizes[[k]])
      ind <- setdiff(ind, foldset[[k]])
    }
    foldset[[nfold]] <- ind
    
    Xtra <- Ytra <- list()
    Xtest <- Ytest <- list()
    
    p <- ncol(X)
    q <- ncol(Y)
    Ymean <- apply(Y,2,mean)
    Ycenter <- Y-matrix(Ymean,byrow=T,n,q)
    Xmean <- apply(X,2,mean)
    Xcenter <- X-matrix(Xmean,n,p,byrow=T)
    Xcenternorm <- apply(Xcenter, 2, function(x){sqrt(sum(x^2)/n)})
    Xstd <- Xcenter / matrix(Xcenternorm,n,p, byrow=T)
    for(k in 1:nfold){
      Xtra[[k]] <- X[-foldset[[k]],]
      Ytra[[k]] <- Y[-foldset[[k]],]
      Xtest[[k]] <- X[foldset[[k]],]
      Ytest[[k]] <- Y[foldset[[k]],]
    }
    ##########################################################################
    
    
    ######################## K-fold validation ##############################
    lamfold <- expand.grid(lamseq, 1:nfold)
    lamfold <- as.matrix(lamfold)
    nlamfold <- nrow(lamfold)
    cverr <- lapply(1:nlamfold, function(x){
      lambda <- lamfold[x,1]
      fold <- lamfold[x,2]
      Ytrax <- Ytra[[fold]]
      Xtrax <- Xtra[[fold]]
      Ytestx <- Ytest[[fold]]
      Xtestx <- Xtest[[fold]]
      
      resx <- ADMM.nuclear(Y=Ytrax, X=Xtrax, B=NULL, L=NULL,
                           eta=eta, lambda = lambda,
                           maxiter=maxiter, tol=tol, standardize = TRUE)
      Bhat <- resx$L
      predY <- cbind(1,Xtestx) %*% Bhat
      avgerr <- mean( (Ytestx -predY)^2 )
      return( avgerr )
    })
    cverr <- unlist(cverr)
    cverrmean <- tapply(cverr, lamfold[,1], mean)
    lamopt <- lamseq[which.min(cverrmean)]
    errors <- cverr[lamfold[,1] == lamopt]
    ##########################################################################
    
    ############################# Fit model ##################################
    finalfit <- ADMM.nuclear(Y=Ycenter, X=Xstd, B=B, L=L, eta=eta,
                             lambda=lamopt, maxiter=maxiter, tol=tol,
                             standardize = FALSE)
    Bhat <- finalfit$L
    Slope <- 1/matrix(Xcenternorm,p,q) * Bhat
    Intercept <- matrix(Ymean, 1, q) -  matrix(Xmean,1,p) %*% Slope
    Bhat <- rbind(Intercept, Slope)
    ###########################################################################
    return( list(B=Bhat, d=finalfit$d, foldset=foldset, lamopt=lamopt,
                 opterr=min(cverrmean), errrorsd= sd(errors)))
  }
  
  cv.twostep <- function(Y,X,auxYlist,auxXlist,B,L,eta,
                         lamseq_w, lamseq_delta,
                         maxiter=100,tol=1e-04,nfold=5){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    K <- length(auxYlist)
    q <- ncol(Y)
    p <- ncol(X)
    n <- nrow(X)
    
    auxYcenterlist <- list()
    for(k in 1:K){
      auxY_k <- auxYlist[[k]]
      auxY_kmean <- apply(auxY_k,2,mean)
      auxY_cneter <- auxY_k - matrix(auxY_kmean, nrow(auxY_k), q,byrow=T)
      auxYcenterlist[[k]] <- auxY_cneter
    }
    auxXcenterlist <- list()
    for(k in 1:K){
      auxX_k <- auxXlist[[k]]
      nk <- nrow(auxX_k)
      auxX_kmean <- apply(auxX_k,2,mean)
      auxX_center <- auxX_k - matrix(auxX_kmean, nk, p,byrow=T)
      auxXcenterlist[[k]] <- auxX_center
    }
    
    auxX <- do.call("rbind", auxXcenterlist)
    auxY <- do.call("rbind", auxYcenterlist)
    auxY <- as.matrix(auxY)
    auxX <- as.matrix(auxX)
    
    Xmean <- apply(X,2,mean)
    Xcenter <- X-matrix(Xmean,n,p,byrow=T)
    Ymean <- apply(Y,2,mean)
    Ycenter <- Y - matrix(Ymean,n,q,byrow=T)
    Xcenter <- as.matrix(Xcenter)
    Ycenter <- as.matrix(Ycenter)
    ##########################################################################
    
    
    
    
    
    #################### Standardization aggregate samples ###################
    allY <- rbind(Ycenter,auxY)
    allX <- rbind(Xcenter,auxX)
    
    nall <- nrow(allY)
    allYmean <- apply(allY,2,mean)
    allYcenter <- allY-matrix(allYmean,byrow=T,nall,q)
    allXmean <- apply(allX,2,mean)
    allXcenter <- allX-matrix(allXmean,nall,p,byrow=T)
    
    ##########################################################################
    cv.nuclear(Y=allY, X=allX, B=B,L=L, eta=eta, lamseq = lamseq_w,
               maxiter=maxiter, tol=tol, nfold=nfold) -> cvW
    What <- cvW$B[-1,]
    Yres <- Ycenter - Xcenter %*% What
    
    cv.nuclear(Y=Yres,X=Xcenter, B=B, L=L, eta=eta, lamseq=lamseq_delta,
               maxiter=maxiter, tol=tol, nfold=nfold) -> rescorrection
    
    Deltahat <- rescorrection$B[-1,]
    intercept <-  matrix(Ymean, 1, q) -  matrix(Xmean,1,p) %*% (Deltahat+What)
    Slope <- Deltahat + What
    Bhat <- rbind(intercept, Slope)
    
    ###########################################################################
    return( list(B=Bhat, Deltahat=Deltahat, What=What,
                 lamoptw = cvW$lamopt, lamoptdelta= rescorrection$lamopt))
  }
  
  
  
  
  res_list <- list()
  for(sim in repstart:repend){
    
    Ytarget_sim <- Ylist[[sim]]
    Xtarget_sim <- Xlist[[sim]]
    p <- ncol(Xtarget_sim)
    q <- ncol(Ytarget_sim)
    auxY_sim <- auxYlist_list[[sim]]
    auxX_sim <- auxXlist_list[[sim]]
    auxvec <- 1:length(auxY_sim)
    
    Yauxmat <- do.call("rbind", auxY_sim)
    Xauxmat <- do.call("rbind", auxX_sim)
    
    ###################################### NR ##################################
    cv.nuclear(Y=Ytarget_sim, X=Xtarget_sim, B=NULL,
               L=NULL,eta=1,lamseq = lamLassoseq,
               tol=1e-04,maxiter=maxit) -> cvLasso
    
    dLasso <- svd(matrix(cvLasso$B[-1,],p,q))$d # singluar values
    Lassoerr <- mean( (cvLasso$B[-1,] -Btrue)^2 ) # RMSE
    
    ################################################################################ FSD and MSD ########################################################################################
    traind_i <- 1:nrow(Ytarget_sim)
    ntra_i <- length(traind_i)
    fold_selection <- 3 # the number of folds to be used to determine the optimal C
    set.seed(100)
    folds_i <- create_folds(traind_i,k=fold_selection, type = "basic")
    errlist_screening <- list()
    errlist_forward <- list()
    Blistlist <- list()
    Bset <- list()
    seset <- list()
    
    ############################### Code to find optimal C in FSD and MSD ####################################################################
    for(k1 in 1:fold_selection){
      Ytrai_tra <- Ytarget_sim[folds_i[[k1]],]
      Ytrai_val <- Ytarget_sim[-folds_i[[k1]],]
      
      Xtrai_tra <- Xtarget_sim[folds_i[[k1]],]
      Xtrai_val <- Xtarget_sim[-folds_i[[k1]],]
      
      set.seed(100)
      ni <- nrow(Ytarget_sim)
      indi <- (1:ni)
      Ytraitraset <- Xtraitraset <- list()
      Ytraitestset <- Xtraitestset <- list()
      traind_k1 <- 1:nrow(Ytrai_tra)
      folds_k1 <- create_folds(traind_k1,k=fold_selection, type = "basic")
      for(ss in 1:fold_selection){
        Ytraitestset[[ss]] <- Ytrai_tra[-folds_k1[[ss]],]
        Ytraitraset[[ss]] <-  Ytrai_tra[folds_k1[[ss]],]
        
        Xtraitestset[[ss]] <- Xtrai_tra[-folds_k1[[ss]],]
        Xtraitraset[[ss]] <-  Xtrai_tra[folds_k1[[ss]],]
        
      }
      
      p <- ncol(Xtraitraset[[ss]])
      q <- ncol(Ytraitraset[[ss]])
      
      foldbyaux <- expand.grid(1:fold_selection, 1:length(auxvec))
      foldbyaux <- as.matrix(foldbyaux)
      nauxfold <- nrow(foldbyaux)
      ######################### To compute validatieon errors when a source is combined with the target
      errorvec <- foreach::foreach(ss=1:nauxfold) %dopar%{
        
        fold <- foldbyaux[ss,1]
        g <- foldbyaux[ss,2]
        
        Ytrai_ss <- rbind(Ytraitraset[[fold]], auxY_sim[[ g ]] )
        Xtrai_ss <- rbind(Xtraitraset[[fold]], auxX_sim[[ g ]] )
        
        
        cv.nuclear(Y=Ytrai_ss, X=Xtrai_ss, B=NULL,
                   L=NULL,eta=1,lamseq = lamLassoseq,
                   tol=1e-04,maxiter=maxit) -> cvLasso_ss
        
        mean((Ytraitestset[[fold]] - cbind(1,Xtraitestset[[fold]]) %*% cvLasso_ss$B)^2)
      }
      ##############################################################################################
      
      
      ################################### Valdiation errors when using only target #################
      lassoerr_vec <- foreach::foreach(ss=1:3) %dopar%{
        Lasso_val <-    cv.nuclear(Y=Ytraitraset[[ss]], X=Xtraitraset[[ss]], B=NULL,
                                   L=NULL,eta=1,lamseq = lamLassoseq, tol=1e-04,maxiter=maxit)
        mean((Ytraitestset[[ss]]  - cbind(1,Xtraitestset[[ss]])  %*%  Lasso_val$B )^2)
      }
      lassoerr_vec <- unlist(lassoerr_vec)
      ###############################################################################################
      
      ################################## To save results for FSD ###################################
      sellist_list <- list()
      errlist <- list()
      ###############################################################################################
      errorvec <- unlist(errorvec)
      errorvecmean <- tapply(errorvec, foldbyaux[,2], mean)
      
      
      ############################## MSD-step ########################################################
      cv.nuclear(Y=Ytrai_tra, X=Xtrai_tra, B=NULL,
                 L=NULL,eta=1,lamseq = lamLassoseq,
                 tol=1e-04,maxiter=maxit) -> cvLasso_tra
      predLassoi <- cbind(1,Xtrai_val) %*% cvLasso_tra$B
      Lassoerr_tra <- mean( (Ytrai_val-predLassoi)^2 )
      
      # selectionres_tra[i]: selection results  when C is set as C[i]
      selectionres_tra <- foreach(gg= 1:length(C)) %dopar%{
        library(foreach)
        auxYlist1 <- auxXlist1 <- list()
        auxvec1 <- which( errorvecmean < ((1+C[gg])*mean(lassoerr_vec)) )
        if(length(auxvec1)==0){
          cv1err <- Lassoerr_tra
        }else{
          for(g in 1:length(auxvec1)){
            auxYlist1[[g]] <- auxY_sim[[ auxvec1[g] ]]
            auxXlist1[[g]] <- auxX_sim[[ auxvec1[g] ]]
          }
          cv1 <- cv.twostep(Y=Ytrai_tra, X=Xtrai_tra, B=NULL,
                            L=NULL,eta=1,
                            lamseq_w = lam_wNaiveseq,
                            lamseq_delta=lam_deltaNaiveseq,
                            tol=1e-04,maxiter=maxit,
                            auxYlist=auxYlist1, auxXlist=auxXlist1)
          cv1_pred <- cbind(1,Xtrai_val) %*% cv1$B
          cv1err <-   mean( (Ytrai_val -    cv1_pred)^2 )
        }
        
        list(err=cv1err)
      }
      #############################################################################################
      
      
      ################################## FSD-step ##################################################
      errorvecmean_init <- errorvecmean
      for(sss in 1:length(C)){
        errorvecmean <- errorvecmean_init
        detection <- any(errorvecmean < ((1+C[sss])*mean(lassoerr_vec)))
        sources <- NULL
        sources_cands <- auxvec
        auxYcand <- auxXcand <- NULL
        if(detection){
          #auxvec_old <- sources_cands
          sources <- auxvec[which.min(errorvecmean)]
          sources_cands <- auxvec[-which.min(errorvecmean)]
          #auxvec_new <- auxvec_new[-which.min(errorvecmean)]
          auxYcand <- auxY_sim[[which.min(errorvecmean)]]
          auxXcand <- auxX_sim[[which.min(errorvecmean)]]
          for(kkk in 1:(length(auxvec)-1)){
            
            
            cri_err <- errorvecmean[which.min(errorvecmean)]
            foldbyaux2 <- expand.grid(1:fold_selection, sources_cands  )
            foldbyaux2 <- as.matrix(foldbyaux2)
            nauxfold2 <- nrow(foldbyaux2)
            errorvec <- foreach::foreach(ss=1:nauxfold2) %dopar%{
              
              fold <- foldbyaux2[ss,1]
              g <- foldbyaux2[ss,2]
              auxYcand_g <- rbind(auxYcand, auxY_sim[[ g ]])
              auxXcand_g <- rbind(auxXcand, auxX_sim[[ g ]])
              
              Ytrai_ss <- rbind(Ytraitraset[[fold]], auxYcand_g )
              Xtrai_ss <- rbind(Xtraitraset[[fold]], auxXcand_g )
              
              val_g <- cv.nuclear(Y=Ytrai_ss, X=Xtrai_ss, B=NULL,L=NULL,eta=1,lamseq = lamLassoseq,
                                  tol=1e-04,maxiter=maxit)
              
              mean((Ytraitestset[[fold]]  - cbind(1,Xtraitestset[[fold]])  %*%  val_g$B )^2)
            }
            errorvec <- unlist(errorvec)
            errorvecmean <- tapply(errorvec, foldbyaux2[,2], mean)
            if(all( errorvecmean > ( (1+C[sss]) * cri_err)  )){
              break
            }else{
              sources <- c(sources, sources_cands[which.min(errorvecmean)])
              auxYcand <- rbind(auxYcand, auxY_sim[[ sources_cands[which.min(errorvecmean)] ]])
              auxXcand <- rbind(auxXcand, auxX_sim[[ sources_cands[which.min(errorvecmean)] ]])
              sources_cands <- sources_cands[-which.min(errorvecmean)]
            }
            
          }
          
          selist <- sources
          auxvec1 <- (1:length(auxvec))[which(auxvec %in% selist )]
          auxYlist1 <- auxXlist1 <- list()
          if(length(auxvec1)==0){
            cv1err <- Lassoerr_tra
          }else{
            for(g in 1:length(auxvec1)){
              auxYlist1[[g]] <- auxY_sim[[ auxvec1[g] ]]
              auxXlist1[[g]] <- auxX_sim[[ auxvec1[g] ]]
            }
          }
          
          cv1 <- cv.twostep(Y=Ytrai_tra, X=Xtrai_tra, B=NULL,
                            L=NULL,eta=1,
                            lamseq_w = lam_wNaiveseq,
                            lamseq_delta=lam_deltaNaiveseq,
                            tol=1e-04,maxiter=maxit,
                            auxYlist=auxYlist1, auxXlist=auxXlist1)
          cv1_pred <- cbind(1,Xtrai_val) %*% cv1$B
          cv1err <-   mean( (Ytrai_val -    cv1_pred)^2 )
          errlist[sss] <- cv1err
          sellist_list[[sss]] <- selist
          Blistlist[[sss]] <- cv1$B
          
          
        }else{
          errlist[sss] <- Lassoerr_tra
          sellist_list[[sss]] <- NULL
        }
        
      }
      ########################################################################################################
      errlist_screening[[k1]] <- unlist(selectionres_tra)
      errlist_forward[[k1]] <- errlist
      
      seset[[k1]] <- sellist_list
      Bset[[k1]] <- Blistlist
      
      
      
    }
    
    err_screening_mean <- apply(do.call('rbind', errlist_screening),2,mean) # err_screening_mean[i]: average of CVerros for MSD when C is set as C[i]
    err_forward_mean <- apply(do.call('rbind', lapply(errlist_forward, unlist)),2,mean) # err_ofrward_mean[i]: average of CVerros for FSD when C is set as C[i]
    
    ################################ Copmuting B using optimal Cs ####################################################
    
    #################### Dividing data to perform source selections baed on the optimal C ###########################
    Ytraitraset <- Xtraitraset <- list()
    Ytraitestset <- Xtraitestset <- list()
    for(ss in 1:fold_selection){
      Ytraitestset[[ss]] <- Ytarget_sim[-folds_i[[ss]],]
      Ytraitraset[[ss]] <-  Ytarget_sim[folds_i[[ss]],]
      
      Xtraitestset[[ss]] <- Xtarget_sim[-folds_i[[ss]],]
      Xtraitraset[[ss]] <-  Xtarget_sim[folds_i[[ss]],]
      
    }
    #################################################################################################################
    foldbyaux <- expand.grid(1:fold_selection, 1:length(auxvec))
    foldbyaux <- as.matrix(foldbyaux)
    nauxfold <- nrow(foldbyaux)
    
    ################# Marginal source detection #########################################
    errorvec <- foreach::foreach(ss=1:nauxfold) %dopar%{
      
      fold <- foldbyaux[ss,1]
      g <- foldbyaux[ss,2]
      
      Ytrai_ss <- rbind(Ytraitraset[[fold]], auxY_sim[[ g ]] )
      Xtrai_ss <- rbind(Xtraitraset[[fold]], auxX_sim[[ g ]] )
      
      cv.nuclear(Y=Ytrai_ss, X=Xtrai_ss, B=NULL,
                 L=NULL,eta=1,lamseq = lamLassoseq,
                 tol=1e-04,maxiter=maxit) -> cvLasso_ss
      
      mean((Ytraitestset[[fold]] - cbind(1,Xtraitestset[[fold]]) %*% cvLasso_ss$B)^2)
    }
    errorvec <- unlist(errorvec)
    errorvecmean <- tapply(errorvec, foldbyaux[,2], mean)
    
    lassoerr_vec <- foreach::foreach(ss=1:3) %dopar%{
      Lasso_val <-    cv.nuclear(Y=Ytraitraset[[ss]], X=Xtraitraset[[ss]], B=NULL,
                                 L=NULL,eta=1,lamseq = lamLassoseq,
                                 tol=1e-04,maxiter=maxit)
      mean((Ytraitestset[[ss]]  - cbind(1,Xtraitestset[[ss]])  %*%  Lasso_val$B )^2)
      
    }
    lassoerr_vec <- unlist(lassoerr_vec)
    
    opt_C <- C[which.min( unlist(err_screening_mean) )]
    auxvec1_opt <- which( errorvecmean < ((1+opt_C)*mean(lassoerr_vec)) )
    auxvec1_opt_screening <- auxvec1_opt
    if(length(auxvec1_opt)==0){
      cv1err_opt <- Lassoerr
      d_cv1_opt <- dLasso
      What_1_opt <- cvLasso$B
      Delthat_1_opt <- matrix(0,p,q)
      #ell1_cv1err_opt <- ell1_Lassoerr
    }else{
      auxvec1 <- which( errorvecmean < ((1+opt_C)*mean(lassoerr_vec)) )
      auxYlist1 <- auxXlist1 <- list()
      for(g in 1:length(auxvec1)){
        auxYlist1[[g]] <- auxY_sim[[ auxvec1[g] ]]
        auxXlist1[[g]] <- auxX_sim[[ auxvec1[g] ]]
      }
      cv1_opt <- cv.twostep(Y=Ytarget_sim, X=Xtarget_sim, B=NULL,
                            L=NULL,eta=1,
                            lamseq_w = lam_wNaiveseq,
                            lamseq_delta=lam_deltaNaiveseq,
                            tol=1e-04,maxiter=maxit,
                            auxYlist=auxYlist1, auxXlist=auxXlist1)
      d_cv1_opt <- svd(cv1_opt$B[-1,])$d
      What_1_opt <- cv1_opt$W
      Delthat_1_opt <- cv1_opt$Delta
      cv1err_opt <-   mean( (cv1_opt$B[-1,] -Btrue)^2 )
    }
    ##################################################################################################
    
    ################################ FSD ########################################################
    optC_sequence <-  C[which.min( unlist(err_forward_mean) )]
    detection <- any(errorvecmean < ((1+optC_sequence)*mean(lassoerr_vec)))
    sources_cands <- auxvec
    auxYcand <- auxXcand <- NULL
    selist_opt <- NULL
    if(detection){
      auxvec_old <- sources_cands
      sources <- auxvec[which.min(errorvecmean)]
      sources_cands <- auxvec[-which.min(errorvecmean)]
      #auxvec_new <- auxvec_new[-which.min(errorvecmean)]
      
      auxYcand <- auxY_sim[[which.min(errorvecmean)]]
      auxXcand <- auxX_sim[[which.min(errorvecmean)]]
      selist_opt <- sources
      
      for(kkk in 1:(length(auxvec)-1)){
        
        cri_err <- errorvecmean[which.min(errorvecmean)]
        foldbyaux2 <- expand.grid(1:fold_selection, sources_cands  )
        foldbyaux2 <- as.matrix(foldbyaux2)
        nauxfold2 <- nrow(foldbyaux2)
        errorvec <- foreach::foreach(ss=1:nauxfold2) %dopar%{
          
          fold <- foldbyaux2[ss,1]
          g <- foldbyaux2[ss,2]
          auxYcand_g <- rbind(auxYcand, auxY_sim[[ g ]])
          auxXcand_g <- rbind(auxXcand, auxX_sim[[ g ]])
          
          Ytrai_ss <- rbind(Ytraitraset[[fold]], auxYcand_g )
          Xtrai_ss <- rbind(Xtraitraset[[fold]], auxXcand_g )
          
          val_g <- cv.nuclear(Y=Ytrai_ss, X=Xtrai_ss, B=NULL,L=NULL,eta=1,lamseq = lamLassoseq,
                              tol=1e-04,maxiter=maxit)
          
          mean((Ytraitestset[[fold]]  - cbind(1,Xtraitestset[[fold]])  %*%  val_g$B )^2)
        }
        errorvec <- unlist(errorvec)
        errorvecmean <- tapply(errorvec, foldbyaux2[,2], mean)
        if(all( errorvecmean > ( (1+optC_sequence) * cri_err)  ) ){
          break
        }else{
          sources <- c(sources, sources_cands[which.min(errorvecmean)])
          auxYcand <- rbind(auxYcand, auxY_sim[[ sources_cands[which.min(errorvecmean)] ]])
          auxXcand <- rbind(auxXcand, auxX_sim[[ sources_cands[which.min(errorvecmean)] ]])
          sources_cands <- sources_cands[-which.min(errorvecmean)]
          selist_opt <- sources
        }
        
      }
    }
    
    auxvec1_opt <- (1:length(auxvec))[auxvec %in% selist_opt ]
    if(length(auxvec1_opt)==0){
      cv1err_opt_seq <- Lassoerr
      d_cv1_opt <- dLasso
      What_1_opt <- cvLasso$B
      Delthat_1_opt <- matrix(0,p,q)
      #cv1_pred_opt_seq <- predLassoi
      #ell1_cv1err_opt <- ell1_Lassoerr
    }else{
      auxYlist1 <- auxXlist1 <- list()
      for(g in 1:length(auxvec1_opt)){
        auxYlist1[[g]] <- auxY_sim[[ auxvec1_opt[g] ]]
        auxXlist1[[g]] <- auxX_sim[[ auxvec1_opt[g] ]]
      }
      cv1_opt <- cv.twostep(Y=Ytarget_sim, X=Xtarget_sim, B=NULL,
                            L=NULL,eta=1,
                            lamseq_w = lam_wNaiveseq,
                            lamseq_delta=lam_deltaNaiveseq,
                            tol=1e-04,maxiter=maxit,
                            auxYlist=auxYlist1, auxXlist=auxXlist1)
      d_cv1_opt <- svd(cv1_opt$B[-1,])$d
      What_1_opt <- cv1_opt$W
      Delthat_1_opt <- cv1_opt$Delta
      cv1err_opt_seq <-    mean( (cv1_opt$B[-1,] -Btrue)^2 )
    }
    
    #####################################################################################################
    
    ###################################################################################################################################################################
    
    cvNaive <- cv.twostep(Y=Ytarget_sim, X=Xtarget_sim, B=NULL,
                          L=NULL,eta=1,
                          lamseq_w = lam_wNaiveseq,
                          lamseq_delta=lam_deltaNaiveseq,
                          tol=1e-04,maxiter=maxit,
                          auxYlist=auxY_sim, auxXlist=auxX_sim) # cross-validation for two-step transfer learning
    
    
    d_naitrans <- svd(cvNaive$B[-1,])$d
    What_naive <- cvNaive$W
    Delthat_naive <- cvNaive$Delta
    NaiveTransfer_err <-  mean( (cvNaive$B[-1,] -Btrue)^2 ) # RMSE of two-step transfer learning
    
    AllYtrai <- rbind(Ytarget_sim, do.call('rbind',auxY_sim) )
    AllXtrai <- rbind(Xtarget_sim, do.call('rbind',auxX_sim) )
    cvAll <- cv.nuclear(Y=AllYtrai, X=AllXtrai, B=NULL,
                        L=NULL,eta=1,lamseq = lamLassoseq,
                        tol=1e-04,maxiter=maxit)
    ALL_cer <-  mean( (cvAll$B[-1,] -Btrue)^2 ) # RMSE of pooled-NR
    d_all <- svd(cvAll$B[-1,])$d
    
    target_fh_cverr_i <- c(Lassoerr,
                           ALL_cer, NaiveTransfer_err, cv1err_opt, cv1err_opt_seq )
    
    
    
    # error: estimation error
    # dectected_sources_foward: sources detected through FSD
    # detected_sources_screening: sources detected through MSD
    # optC_sequence: C selected by CV when FSD is considered 
    # opt_C_screening: C selected by CV when MSD is considered 
    
    res_list[[sim]] <- list(error=target_fh_cverr_i,
                            dectected_sources_foward= auxvec1_opt,
                            detected_sources_screening = auxvec1_opt_screening, 
                            optC_sequence=optC_sequence,
                            opt_C_screening=opt_C)
    
  }
  res_list
}

Summary_ft_detection <- function(results, numinf=3){
  len_res <- length(results)
  Num_cor <- list()
  truesupport=c(1)
  if(numinf==3){
    truesupport = c(1,2,3)
  }
  for(aaa in 1:len_res){
    detect_aaa_forward <- results[[aaa]]$dectected_sources_foward
    detect_aaa_screening <- results[[aaa]]$detected_sources_screening
    Num_cor_aaa <- NULL
    if(length(detect_aaa_forward) == length(truesupport)){
      Num_cor_aaa[1] <- ifelse(  all(detect_aaa_forward==truesupport), TRUE, FALSE )
    }else{
      Num_cor_aaa[1] <- FALSE
    }
    
    if(length(detect_aaa_screening) == length(truesupport)){
      Num_cor_aaa[2] <- ifelse(  all(detect_aaa_screening==truesupport), TRUE, FALSE )
    }else{
      Num_cor_aaa[2] <- FALSE
    }
    Num_cor[[aaa]] <- Num_cor_aaa
  }
  
  
  mean_pci = round(apply( do.call('rbind', Num_cor), 2, mean),3)
  mean_card = c(
    mean(unlist(lapply(results, function(x){length(x$dectected_sources_foward)})))  ,
    mean(unlist(lapply(results, function(x){length(x$detected_sources_screening)}))  ) )
  mean_pci = paste0(mean_pci*100,"%")
  sd_card = c(
    sd(unlist(lapply(results, function(x){length(x$dectected_sources_foward)})))  ,
    sd(unlist(lapply(results, function(x){length(x$detected_sources_screening)}))  ) )
  sd_card = round(sd_card,3)
  sd_card = paste0("(", sd_card, ")")
  mean_sd_card <- paste0(mean_card, sd_card)
  
  res <- c(mean_pci[1], mean_sd_card[1], mean_pci[2], mean_sd_card[2])
  names(res) <- c("PCI:FSD", "Card:FSD", "PCI:MSD", "Card:MSD")
  #res[c(1,3)] <- paste0(res[c(1,3)],"&")
  #res[2] <- paste0(res[2], "&&")
  res
}

Summary_ft_est <- function(results){
  len_res <- length(results)
  mean_ee <- round( apply(do.call('rbind', lapply(results, function(x){sqrt(x$error)})),2,mean),3)
  sd_ee <- round( apply(do.call('rbind', lapply(results, function(x){sqrt(x$error)})),2,sd),3)
  sd_ee <- paste0("(", sd_ee, ")")
  mean_sd_ee <- paste0(mean_ee, sd_ee)
  #mean_sd_ee[1:4] <- paste0(mean_sd_ee[1:4], "&")
  names(mean_sd_ee) <- c("NR", "Pooled-NR", "[K]-Trans", "MSD", "FSD")
  mean_sd_ee
  
}


generator_target_simulsourcedetection <- function(rep,n,p,q,r, corx=0.5,cory=0){
  # rep: the number of replicates
  # n: sample size, p: the number of covariates, q: the number of responses
  # r: rank of B
  # corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
  # cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
  set.seed(100)
  B1 <- runif(p*r,-1,1)
  B1 <- matrix(B1,p,r)
  B2 <- runif(r*q,-1,1)
  B2 <- matrix(B2,r,q)
  B <- B1 %*% B2
  
  
  Sigma_e <- matrix(0,q,q)
  for(i1 in 1:q){
    for(i2 in 1:q){
      Sigma_e[i1,i2] <- cory^{abs(i1-i2)}
    }
  }
  diag(Sigma_e) <- 1
  
  
  Sigma_x <- matrix(0,p,p)
  for(i1 in 1:p){
    for(i2 in 1:p){
      Sigma_x[i1,i2] <- corx^{abs(i1-i2)}
    }
  }
  
  set.seed(100)
  Xlist <- Ylist <- list()
  for(sim in 1:rep){
    E <- MASS::mvrnorm(n=n,mu=rep(0,q),Sigma=Sigma_e)
    X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma_x)
    Y <- X %*% B + E
    Ylist[[sim]] <- Y
    Xlist[[sim]] <- X
  }
  return(list(Ylist=Ylist,Xlist=Xlist,B=B))
  
}



generator_source_simulsourcedetection <-  function(rep, n,p,q, h,numsource, 
                                                   rank_source, B, corx=0.5,cory=0){
  
  # rep: the number of replicates
  # n: sample size, p: the number of covariates, q: the number of responses
  # h: transferring level
  # numsource: the number of source dataset with h transferring level
  # rank_source: rank of contrast
  # B: target parameter
  # corx: [Sigma_x]^{|i-j|}=corx^{|i-j|}
  # cory: [Sigma_epsilon]^{|i-j|}=Sigma_epsilon^{|i-j|}
  
  Sigma_e <- matrix(0,q,q)
  for(i1 in 1:q){
    for(i2 in 1:q){
      Sigma_e[i1,i2] <- cory^{abs(i1-i2)}
    }
  }
  diag(Sigma_e) <- 1
  Sigma_x <- matrix(0,p,p)
  for(i1 in 1:p){
    for(i2 in 1:p){
      Sigma_x[i1,i2] <- corx^{abs(i1-i2)}
    }
  }
  
  set.seed(100)
  auxYlist_list <- auxXlist_list <- list()
  
  
  
  Delta_list <- list()
  for(kk in 1:numsource){
    A <- runif(p*q,-0.5,0.5)
    A <- matrix(A,p,q)
    svd_B <- svd(A)
    dd1 <- rep(0,min(p,q)-rank_source)
    dd2 <- h/sum(1:rank_source)*(1:rank_source)
    dd <- c(rev(dd2), dd1)
    Delta_list[[kk]] <- svd_B$u %*% diag(dd) %*% t(svd_B$v)
  }
  
  
  for(sim in 1:rep){
    auxYlist_sim <- list()
    auxXlist_sim <- list()
    for(k in 1:numsource){
      E <- MASS::mvrnorm(n=n,mu=rep(0,q),Sigma=Sigma_e)
      X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma_x)
      Y <- X %*% (B+Delta_list[[k]]) + E
      auxYlist_sim[[k]] <- Y
      auxXlist_sim[[k]] <- X
    }
    
    
    auxYlist_list[[sim]] <- auxYlist_sim
    auxXlist_list[[sim]] <- auxXlist_sim
    
  }
  
  return(list(Delta_list=Delta_list, auxYlist_list=auxYlist_list, auxXlist_list=auxXlist_list, numsource=numsource,
              d=dd))
  
}




############################# Heoterogeneous design #########################################


generator_source_simulsourcedetection_het <- function(rep, nvec, p, q, hvec, numsource, ranksourcevec, B,
                                                      cory=0){
   Sigma_e <- matrix(0, q, q)
    for (i1 in 1:q) {
        for (i2 in 1:q) {
            Sigma_e[i1, i2] <- cory^{
                abs(i1 - i2)
            }
        }
    }
    Sigmaxlist <- list()

    for(ss in 1:numsource){
        Sigma_x <- matrix(0, p, p)
        for (i1 in 1:p) {
          for (i2 in 1:p) {
            Sigma_x[i1, i2] <- (ss/(numsource+2))^{abs(i1 - i2) }
          }
        }
        Sigmaxlist[[ss]] <- Sigma_x
    }
    set.seed(100)
    auxYlist_list <- auxXlist_list <- list()
    Delta_list <- list()
    for (kk in 1:numsource) {
        A <- runif(p * q, -0.5, 0.5)
        A <- matrix(A, p, q)
        svd_B <- svd(A)
        dd1 <- rep(0, min(p, q) - ranksourcevec[kk])
        dd2 <- hvec[kk]/sum(1:ranksourcevec[kk]) * (1:ranksourcevec[kk])
        dd <- c(rev(dd2), dd1)
        Delta_list[[kk]] <- svd_B$u %*% diag(dd) %*% t(svd_B$v)
    }
    for(sim in 1:rep){
        auxYlist_sim <- list()
        auxXlist_sim <- list()
        for (k in 1:numsource) {
            E <- MASS::mvrnorm(n = nvec[k], mu = rep(0, q), Sigma = Sigma_e)
            X <- MASS::mvrnorm(n = nvec[k], mu = rep(0, p), Sigma = Sigmaxlist[[k]])
            Y <- X %*% (B + Delta_list[[k]]) + E
            auxYlist_sim[[k]] <- Y
            auxXlist_sim[[k]] <- X
        }
        auxYlist_list[[sim]] <- auxYlist_sim
        auxXlist_list[[sim]] <- auxXlist_sim
    }
    return(list(Delta_list = Delta_list, auxYlist_list = auxYlist_list,
        auxXlist_list = auxXlist_list, numsource = numsource,
        d = dd))

}



Summary_ft_detection_het <- function (results, truesupport) {
  len_res <- length(results)
  Num_cor <- list()
  #truesupport = 1:numinf
  for (aaa in 1:len_res) {
    detect_aaa_forward <- results[[aaa]]$dectected_sources_foward
    detect_aaa_screening <- results[[aaa]]$detected_sources_screening
    Num_cor_aaa <- NULL
    if (length(detect_aaa_forward) == length(truesupport)) {
      Num_cor_aaa[1] <- ifelse(all(detect_aaa_forward ==
                                     truesupport), TRUE, FALSE)
    }
    else {
      Num_cor_aaa[1] <- FALSE
    }
    if (length(detect_aaa_screening) == length(truesupport)) {
      Num_cor_aaa[2] <- ifelse(all(detect_aaa_screening ==
                                     truesupport), TRUE, FALSE)
    }
    else {
      Num_cor_aaa[2] <- FALSE
    }
    Num_cor[[aaa]] <- Num_cor_aaa
  }
  mean_pci = round(apply(do.call("rbind", Num_cor), 2, mean),
                   3)
  mean_card = c(mean(unlist(lapply(results, function(x) {
    length(x$dectected_sources_foward)
  }))), mean(unlist(lapply(results, function(x) {
    length(x$detected_sources_screening)
  }))))
  mean_pci = paste0(mean_pci * 100, "%")
  sd_card = c(sd(unlist(lapply(results, function(x) {
    length(x$dectected_sources_foward)
  }))), sd(unlist(lapply(results, function(x) {
    length(x$detected_sources_screening)
  }))))
  sd_card = round(sd_card, 3)
  sd_card = paste0("(", sd_card, ")")
  mean_sd_card <- paste0(mean_card, sd_card)
  res <- c(mean_pci[1], mean_sd_card[1], mean_pci[2], mean_sd_card[2])
  names(res) <- c("PCI:FSD", "Card:FSD", "PCI:MSD", "Card:MSD")
  res
}


generator_target_normal<- function(rep,n,p,q,r, corx=0.5,cory=0){
  set.seed(100)
  B1 <- rnorm(p*r)
  B2 <- rnorm(r*q)
  B <- matrix(B1,p,r) %*% matrix(B2,r,q) * 0.2


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
  rx=40
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

generator_source_rankestimationsimul <- function(rep, n,p,q,h, numsource, ranksource, B, corx=0.5,cory=0){
  
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
  #dd <- runif(min(p,q),0,1)
  #dd <- sort(dd, decreasing = TRUE)
  
  
  
  A <- runif(p*q,-0.5,0.5)
  A <- matrix(A,p,q)
  svd_A <- svd(A)
  
  dd1 <- rep(0,min(p,q)-ranksource)
  dd2 <- h/sum(1:ranksource)*(1:ranksource)
  dd <- c(rev(dd2), dd1)
  Delta <- svd_A$u %*% diag(dd) %*% t(svd_A$v)
  
  for(sim in 1:rep){
    auxYlist_sim <- list()
    auxXlist_sim <- list()
    for(k in 1:numsource){
      E <- MASS::mvrnorm(n=n,mu=rep(0,q),Sigma=Sigma_e)
      X <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma_x)
      Y <- X %*% (B-Delta) + E
      auxYlist_sim[[k]] <- Y
      auxXlist_sim[[k]] <- X
    }
    
    
    auxYlist_list[[sim]] <- auxYlist_sim
    auxXlist_list[[sim]] <- auxXlist_sim
    
  }
  
  return(list(Delta=Delta, auxYlist_list=auxYlist_list, auxXlist_list=auxXlist_list, numsource=numsource,
              d=dd))
  
}


################## Common parameters ########################
# Ylist: contains the simulated Y
# Xlist: contains the simulated X
# auxYlist_list: contains the simulated auxYlist, wheer auxYlist is the list containing source response matrices
# auxXlist_list:  contains the simulated auxXlist, wheer auxYlist is the list containing source design matrices
# numsim: the number of replicates, i.e., length(Ylist)
# lamseq_w: candiates of tuning parameters for the firs step
# lamseq_delta: candiates of tuning parameters for the bias-correction step

# maxiter_inital: the maxium of the number of iterations in the ADMM algorithm for computing initial estimate
# tol_inital: tolerance error in the ADMM algorithm for computing initial estimate
# maxiter_biascorrection: the maxium of the number of iterations in the ADMM algorithm for computing contrast
# tol_inital: tolerance error in the ADMM algorithm for computing contrast
#############################################################

########### Computing BICpath for Trans-NR to multiple simulated example ##################
Repfit_TransNR_BIC <- function(Ylist,Xlist, auxYlist_list, auxXlist_list, lamseq_w, lamseq_delta, numsim,
                               maxiter_inital=100, maxiter_biascorrection=300, 
                               tol_inital=1e-04, tol_biascorrection=1e-04){
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
  
  
  
  path.biascorrectionNR <- function(Y,X,auxYlist,auxXlist,B,L,eta,
                                    lamseq_w, lamseq_delta,
                                    maxiter=100,tol=1e-04){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    K <- length(auxYlist)
    q <- ncol(Y)
    p <- ncol(X)
    n <- nrow(X)
    ################## Centering each data #################################
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
    
    
    
    
    
    #################### Standardization aggregated samples ###################
    allY <- rbind(Ycenter,auxY)
    allX <- rbind(Xcenter,auxX)
    
    nall <- nrow(allY)
    allYmean <- apply(allY,2,mean)
    allYcenter <- allY-matrix(allYmean,byrow=T,nall,q)
    allXmean <- apply(allX,2,mean)
    allXcenter <- allX-matrix(allXmean,nall,p,byrow=T)
    
    ##########################################################################
    lamseq_pair <- expand.grid(lamseq_w, lamseq_delta)
    lamseq_pair <- as.matrix(lamseq_pair)
    nallpairs <- nrow(lamseq_pair)
    pathres <- lapply(1:nallpairs, function(x){
      lambda_w <- lamseq_pair[x,1]
      lambda_delta <- lamseq_pair[x,2]
      
      ADMM.nuclear(Y=allY, X=allX, B=B,L=L, eta=eta, lambda = lambda_w,
                   maxiter=maxiter, tol=tol) -> fit_W
      
      Yres <- Ycenter - Xcenter %*% fit_W$L[-1,]
      
      ADMM.nuclear(Y=Yres, X=Xcenter, B=B,L=L, eta=eta, lambda=lambda_delta, maxiter=maxiter_inital,
                   tol=tol_inital, standardize=TRUE) -> fit_delta
      
      intercept <-  matrix(Ymean, 1, q) -  matrix(Xmean,1,p) %*% ( fit_delta$L[-1,]+fit_W$L[-1,])
      Slope <-  fit_delta$L[-1,]+fit_W$L[-1,]
      Bhat <- rbind(intercept, Slope)
      Bhat
    })
    
    
    
    ###########################################################################
    return( pathres )
  }
  
  n0 <- nrow(Ylist[[1]])
  q <- ncol(Ylist[[1]])
  p <- sum(svd(Xlist[[1]])$d> 1e-02) # rank of x
  foreach::foreach(ss = 1:numsim) %dopar%{
    path.biascorrectionNR(Y=Ylist[[ss]], X=Xlist[[ss]],
                          auxYlist=auxYlist_list[[ss]],
                          auxXlist = auxXlist_list[[ss]],
                          B=NULL, L=NULL,eta=1, lamseq_w = lamseq_w,lamseq=lamseq_delta,
                          maxiter=maxiter_biascorrection,
                          tol=tol_biascorrection) -> pathres
    
    residual_list_ss <-   lapply(pathres, function(x){ ( (Ylist[[ss]] - cbind(1, Xlist[[ss]])  %*% x ) )  } )
    det_ss <- unlist(lapply(residual_list_ss, function(x){ sum(x^2)})     ) # sum of squared errors
    rankseq <- unlist(lapply(pathres, function(x){ sum(corpcor::fast.svd(x[-1,])$d > 1e-02) }))
    df <- rankseq*(p+q) - rankseq^2
    logdet <- log(det_ss)
    penaltyterm_1 <- ( log(n0*q)/(n0*q) ) * df
    penaltyterm_2 <- log(log(n0))/(n0*q)*log(p)  * df
    
    
    BICs <- (logdet + penaltyterm_2)
    # #remove candidates corresponding to zero estimate and full-rank estimate
    BICs_nozero <- BICs[rankseq!=0 & rankseq!=min(p,q)]
    optwhich_nozero <- which.min(BICs_nozero)
    optwhich_BIC <- which(BICs== BICs_nozero[optwhich_nozero])[1]
    
    
    optrank_BIC <- rankseq[optwhich_BIC]
    optB_BIC <- pathres[[optwhich_BIC]]
    
    list(optrank_BIC=optrank_BIC, B=optB_BIC, optwhich_BIC=optwhich_BIC,
         BICs=BICs, rankseq=rankseq)
  }
  
}


Repfit_TransSCAD_BIC <- function(Ylist, Xlist,
                                 auxYlist_list, auxXlist_list,
                                 lamseq_w, lamseq_delta,numsim,
                                 maxiter_inital=100, maxiter_biascorrection=300,
                                 tol_inital=1e-04, tol_biascorrection=1e-04, nfold=5){

  ######################## Functions for the one-step algorithm ############################
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
  ########################################################################################


  ADMM.biascorrection.scad_onestep <- function(Y,X,What,a=3.7, Delta,H,Pi,eta,lambda,maxiter,tol=1e-04,
                                               standardize=TRUE){

    # eta: rho in the manuscript
    mu <- lambda/2

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

      What <- matrix(Xcenternorm,p,q) * What
    }

    d_What <- corpcor::fast.svd(What)$d
    if(length(d_What) < min(p,q)){
      d_What = c(d_What, rep(0, min(p,q)- length(d_What)))
    }
    # drivative of scad penalty
    shrinkageterms_pi <- ifelse(d_What <= mu, mu, ifelse(d_What < (a*mu), (a*mu-d_What)/(a-1), 0)  )

    gramridge <- 2*diag(eta,p,p) + (1/n)*t(X)%*%X
    gramridgeinv <- chol2inv(chol(gramridge))
    covXY <- (1/n)*t(X)%*%Y

    iter <- 1
    isstop <- FALSE
    obj <- NULL
    errDelta <- errH <- errPi <- errZ1 <- errZ2 <- NULL

    if(is.null(Delta)){
      Delta <- matrix(0,p,q)
    }
    if(is.null(H)){
      H <- matrix(0,p,q)
    }
    if(is.null(Pi)){
      Pi <- matrix(0,p,q)
    }
    Z1 <- Delta-H
    Z2 <- Delta+What-Pi



    while(iter <= maxiter & !(isstop)){
      Deltaold <- Delta; Hold <- H; Piold <- Pi
      Z1old <- Z1; Z2old <- Z2

      ########## Delta-step #########################
      Deltapred <- covXY - (Z1 + Z2) + eta * (H+Pi-What)
      Delta <- gramridgeinv %*% Deltapred

      ########## H-step #########################
      Hpred <- Delta + Z1/eta
      svd_H <- corpcor::fast.svd(Hpred)
      dthres <- pmax(svd_H$d-lambda/eta,0)
      H <- svd_H$u %*% diag(dthres) %*% t(svd_H$v)

      ############### Pi-step ######################
      Pipred <- Delta + What + Z2/eta
      svd_Pi <- corpcor::fast.svd(Pipred)
      len_d <- length(svd_Pi$d)
      dthres_pi <- pmax(svd_Pi$d-shrinkageterms_pi[1:len_d]/eta, 0)
      Pi <- svd_Pi$u %*% diag(dthres_pi) %*% t(svd_Pi$v)

      Z1 <- Z1 + eta*(Delta - H)
      Z2 <- Z2 + eta*(Delta+What - Pi)


      errDelta[iter] <-  sum((Deltaold -Delta)^2)
      errH[iter] <-  sum((Hold -H)^2)
      errPi[iter] <-  sum((Piold -Pi)^2)
      errZ1[iter] <-  sum((Z1old -Z1)^2)
      errZ2[iter] <-  sum((Z2old -Z2)^2)
      isstop <- (max(errDelta[iter],errH[iter],errPi[iter],errZ1[iter],errZ2[iter]) <= tol)
      iter <- iter + 1
    }

    B <- Pi
    Delta <- H
    if(standardize){
      B <- matrix(1/Xcenternorm,p,q) * B # rescale
      Bint <- matrix(Ymean,1,q) - matrix(Xmean,1,p) %*% B
      B <- rbind(Bint,B)
    }

    return(list(B=B, Delta=Delta, obj=obj,errDelta=errDelta,errH=errH,errPi=errPi,
                errZ1=errZ1, errZ2=errZ2,
                iter=iter-1,d=dthres_pi))

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
  if( (n %% nfold) > 0 ){
    samplesizes[1:(n %% nfold)] <- samplesizes[1:(n %% nfold)] + 1
  }

  for( k in 1:(nfold-1) ){
    foldset[[k]] <- sample(ind, samplesizes[[k]])
    ind <- setdiff(ind, foldset[[k]])
  }
  foldset[[nfold]] <- ind
  ##############################################################################
  Xtra <- Ytra <- list()
  Xtest <- Ytest <- list()
  ################### Standardization step ################################
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
    if(length(foldset[[k]]) == 1 ){
      Xtest[[k]] <- matrix(X[foldset[[k]],],nrow=1)
      Ytest[[k]] <- matrix(Y[foldset[[k]],],nrow=1)
    }else{
      Xtest[[k]] <- X[foldset[[k]],]
      Ytest[[k]] <- Y[foldset[[k]],]
    }

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


  path.scadbiascorrection_onestep <- function(Y,X,auxYlist,auxXlist,a=3.7, B,L,Delta,H,Pi,eta,
                                              lamseq_w, lamseq_delta,
                                              maxiter_inital=100, maxiter_biascorrection=300,
                                              tol_inital=1e-04, tol_biascorrection=1e-04,
                                              standardize=TRUE){

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

    allY <- rbind(Ycenter,auxY)
    allX <- rbind(Xcenter,auxX)

    nall <- nrow(allY)
    allYmean <- apply(allY,2,mean)
    allYcenter <- allY-matrix(allYmean,byrow=T,nall,q)
    allXmean <- apply(allX,2,mean)
    allXcenter <- allX-matrix(allXmean,nall,p,byrow=T)

    intialcvress <- cv.nuclear(Y=allY, X=allX, B=B,L=L, eta=eta, lamseq = lamseq_w,
                   maxiter=maxiter_inital, tol=tol_inital, nfold=nfold)



    pathres <- lapply(1:length(lamseq_delta), function(x){

      Yres <- Ycenter - Xcenter %*% intialcvress$B[-1,]

      ADMM.biascorrection.scad_onestep(Y=Yres, X=Xcenter,What=intialcvress$B[-1,], a=a, Delta=Delta, H=H, Pi=Pi, eta=eta,
                                       lambda=lamseq_delta[x], maxiter=maxiter_biascorrection,
                                       tol=tol_biascorrection, standardize=TRUE) -> fit_delta
      fit_delta$B
    })


    return(list(pathres=pathres, W=intialcvress$B))
  }


  n0 <- nrow(Ylist[[1]])
  q <- ncol(Ylist[[1]])
  p <- sum(svd(Xlist[[1]])$d> 1e-02)
  foreach::foreach(ss = 1:numsim) %dopar%{

    path.scadbiascorrection_onestep(Y=Ylist[[ss]], X=Xlist[[ss]],
                                    auxYlist=auxYlist_list[[ss]],
                                    auxXlist = auxXlist_list[[ss]],a=3.7, B=NULL, L=NULL, Delta=NULL,H=NULL,Pi=NULL,eta=1,
                                    lamseq_w = lamseq_w,
                                    lamseq=lamseq_delta,
                                    maxiter_inital=300,
                                    maxiter_biascorrection=300,
                                    tol_inital=1e-04,
                                    tol_biascorrection=1e-04) -> pathres

    residual_list_ss <-   lapply(pathres$pathres, function(x){ ( (Ylist[[ss]] - cbind(1, Xlist[[ss]])  %*% x ) )  } )
    det_ss <- unlist(lapply(residual_list_ss, function(x){ sum(x^2) })     )
    rankseq <- unlist(lapply(pathres$pathres, function(x){ sum(corpcor::fast.svd(x[-1,])$d > 1e-02) }))
    df <- rankseq*(p+q) - rankseq^2
    logdet <- log(det_ss) # logarithm of sum of squred errors
    penaltyterm <- log(log(n0))/(n0*q)*log(p)  * df

    BICs <- (logdet + penaltyterm)
    #remove candidates corresponding to zero estimate and full-rank estimate
    BICs_nozero <- BICs[rankseq!=0 & rankseq!=min(p,q)]
    optwhich_nozero <- which.min(BICs_nozero)
    optwhich_BIC <- which(BICs== BICs_nozero[optwhich_nozero])[1]


    optrank_BIC <- rankseq[optwhich_BIC]
    optB_BIC <- lapply(pathres$pathres,function(x){x})[[optwhich_BIC]]
    #optW_BIC <- lapply(pathres,function(x){x$W})[[optwhich_BIC]]
    list(optrank_BIC=optrank_BIC, B=optB_BIC, optwhich_BIC=optwhich_BIC,
         BICs=BICs, rankseq=rankseq, W=pathres$W)

  }

}

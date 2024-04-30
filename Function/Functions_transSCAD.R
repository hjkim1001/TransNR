BIC.TransSCAD <- function(Y, X, auxYlist, auxXlist, lamseq_w, lamseq_delta, eta,
    a = 3.7, B = NULL, L = NULL, Delta = NULL, H = NULL, Pi = NULL,
    maxiter_inital = 100, maxiter_biascorrection = 300, tol_inital = 1e-04,
    tol_biascorrection = 1e-04, standardize = TRUE, nfold=5)
{
    B.update <- function(L, U, eta, covXY, gramridgeInv) {
        B_pred <- covXY - U + eta * L
        return((gramridgeInv %*% B_pred))
    }
    L.update <- function(B, U, eta, lambda, p, q) {
        Lpred <- B + U/eta
        svdres <- corpcor::fast.svd(Lpred)
        d <- svdres$d
        dthres <- pmax(d - (lambda/eta), 0)
        return(list(L = (svdres$u %*% diag(dthres) %*% t(svdres$v)),
            d = dthres))
    }
    ADMM.nuclear <- function(Y, X, B, L, eta, lambda, maxiter,
        tol = 1e-04, standardize = TRUE) {
        n <- nrow(Y)
        p <- ncol(X)
        q <- ncol(Y)
        if (standardize) {
            Ymean <- apply(Y, 2, mean)
            Ycenter <- Y - matrix(Ymean, n, q, byrow = T)
            Xmean <- apply(X, 2, mean)
            Xcenter <- X - matrix(Xmean, n, p, byrow = T)
            Xcenternorm <- apply(Xcenter, 2, function(x) {
                sqrt(sum(x^2)/length(x))
            })
            Xstd <- Xcenter/matrix(Xcenternorm, n, p, byrow = T)
            X <- Xstd
            Y <- Ycenter
        }
        Onematrix <- matrix(1, p, q)
        gramridge <- diag(eta, p, p) + (1/n) * t(X) %*% X
        gramridgeInv <- chol2inv(chol(gramridge))
        covXY <- (1/n) * t(X) %*% Y
        iter <- 1
        isstop <- FALSE
        obj <- NULL
        errB <- errL <- errU <- NULL
        if (is.null(B)) {
            B <- matrix(0, p, q)
        }
        if (is.null(L)) {
            L <- matrix(0, p, q)
        }
        U <- B - L
        while (iter <= maxiter & !(isstop)) {
            Bold <- B
            Lold <- L
            Uold <- U
            B <- B.update(L, U, eta, covXY, gramridgeInv)
            Lres <- L.update(B, U, eta, lambda, p, q)
            L <- Lres$L
            U <- U + eta * (B - L)
            res <- Y - X %*% B
            obj[iter] <- (1/(2 * n)) * sum(res^2) + lambda *
                sum(Lres$d)
            errB[iter] <- sum((Bold - B)^2)
            errL[iter] <- sum((Lold - L)^2)
            errU[iter] <- sum((Uold - U)^2)
            isstop <- (max(errB[iter], errL[iter], errU[iter]) <=
                tol)
            iter <- iter + 1
        }
        if (standardize) {
            B <- matrix(1/Xcenternorm, p, q) * B
            L <- matrix(1/Xcenternorm, p, q) * L
            Bint <- matrix(Ymean, 1, q) - matrix(Xmean, 1, p) %*%
                B
            Lint <- matrix(Ymean, 1, q) - matrix(Xmean, 1, p) %*%
                L
            B <- rbind(Bint, B)
            L <- rbind(Lint, L)
        }
        return(list(B = B, L = L, obj = obj, errB = errB, errL = errL,
            errU = errU, iter = iter - 1, d = Lres$d))
    }
    ADMM.biascorrection.scad_onestep <- function(Y, X, What,
        a = 3.7, Delta, H, Pi, eta, lambda, maxiter, tol = 1e-04,
        standardize = TRUE) {
        mu <- lambda/2
        n <- nrow(Y)
        p <- ncol(X)
        q <- ncol(Y)
        if (standardize) {
            Ymean <- apply(Y, 2, mean)
            Ycenter <- Y - matrix(Ymean, n, q, byrow = T)
            Xmean <- apply(X, 2, mean)
            Xcenter <- X - matrix(Xmean, n, p, byrow = T)
            Xcenternorm <- apply(Xcenter, 2, function(x) {
                sqrt(sum(x^2)/length(x))
            })
            Xstd <- Xcenter/matrix(Xcenternorm, n, p, byrow = T)
            X <- Xstd
            Y <- Ycenter
            What <- matrix(Xcenternorm, p, q) * What
        }
        d_What <- corpcor::fast.svd(What)$d
        if (length(d_What) < min(p, q)) {
            d_What = c(d_What, rep(0, min(p, q) - length(d_What)))
        }
        mu <- 0
        shrinkageterms_pi <- ifelse(d_What <= mu, mu, ifelse(d_What <
            (a * mu), (a * mu - d_What)/(a - 1), 0))
        gramridge <- 2 * diag(eta, p, p) + (1/n) * t(X) %*% X
        gramridgeinv <- chol2inv(chol(gramridge))
        covXY <- (1/n) * t(X) %*% Y
        iter <- 1
        isstop <- FALSE
        obj <- NULL
        errDelta <- errH <- errPi <- errZ1 <- errZ2 <- NULL
        if (is.null(Delta)) {
            Delta <- matrix(0, p, q)
        }
        if (is.null(H)) {
            H <- matrix(0, p, q)
        }
        if (is.null(Pi)) {
            Pi <- matrix(0, p, q)
        }
        Z1 <- Delta - H
        Z2 <- Delta + What - Pi
        while (iter <= maxiter & !(isstop)) {
            Deltaold <- Delta
            Hold <- H
            Piold <- Pi
            Z1old <- Z1
            Z2old <- Z2
            Deltapred <- covXY - (Z1 + Z2) + eta * (H + Pi -
                What)
            Delta <- gramridgeinv %*% Deltapred
            Hpred <- Delta + Z1/eta
            svd_H <- corpcor::fast.svd(Hpred)
            dthres <- pmax(svd_H$d - lambda/eta, 0)
            H <- svd_H$u %*% diag(dthres) %*% t(svd_H$v)
            Pipred <- Delta + What + Z2/eta
            svd_Pi <- corpcor::fast.svd(Pipred)
            len_d <- length(svd_Pi$d)
            dthres_pi <- pmax(svd_Pi$d - shrinkageterms_pi[1:len_d]/eta,
                0)
            Pi <- svd_Pi$u %*% diag(dthres_pi) %*% t(svd_Pi$v)
            Z1 <- Z1 + eta * (Delta - H)
            Z2 <- Z2 + eta * (Delta + What - Pi)
            errDelta[iter] <- sum((Deltaold - Delta)^2)
            errH[iter] <- sum((Hold - H)^2)
            errPi[iter] <- sum((Piold - Pi)^2)
            errZ1[iter] <- sum((Z1old - Z1)^2)
            errZ2[iter] <- sum((Z2old - Z2)^2)
            isstop <- (max(errDelta[iter], errH[iter], errPi[iter],
                errZ1[iter], errZ2[iter]) <= tol)
            iter <- iter + 1
        }
        B <- Pi
        Delta <- H
        if (standardize) {
            B <- matrix(1/Xcenternorm, p, q) * B
            Bint <- matrix(Ymean, 1, q) - matrix(Xmean, 1, p) %*%
                B
            B <- rbind(Bint, B)
            Delta <- matrix(1/Xcenternorm, p, q) * Delta
            Deltaint <- matrix(Ymean, 1, q) - matrix(Xmean, 1, p) %*% Delta
            Delta <- rbind(Deltaint, Delta)
        }
        return(list(B = B, Delta = Delta, obj = obj, errDelta = errDelta,
            errH = errH, errPi = errPi, errZ1 = errZ1, errZ2 = errZ2,
            iter = iter - 1, d = dthres_pi))
    }
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    K <- length(auxYlist)
    q <- ncol(Y)
    p <- ncol(X)
    n <- nrow(X)
    n0 <- n
    auxYcenterlist <- list()
    for (k in 1:K) {
        auxY_k <- auxYlist[[k]]
        auxY_kmean <- apply(auxY_k, 2, mean)
        auxY_cneter <- auxY_k - matrix(auxY_kmean, nrow(auxY_k),
            q, byrow = T)
        auxYcenterlist[[k]] <- auxY_cneter
    }
    auxXcenterlist <- list()
    for (k in 1:K) {
        auxX_k <- auxXlist[[k]]
        nk <- nrow(auxX_k)
        auxX_kmean <- apply(auxX_k, 2, mean)
        auxX_center <- auxX_k - matrix(auxX_kmean, nk, p, byrow = T)
        auxXcenterlist[[k]] <- auxX_center
    }
    auxX <- do.call("rbind", auxXcenterlist)
    auxY <- do.call("rbind", auxYcenterlist)
    auxY <- as.matrix(auxY)
    auxX <- as.matrix(auxX)
    Xmean <- apply(X, 2, mean)
    Xcenter <- X - matrix(Xmean, n, p, byrow = T)
    Ymean <- apply(Y, 2, mean)
    Ycenter <- Y - matrix(Ymean, n, q, byrow = T)
    Xcenter <- as.matrix(Xcenter)
    Ycenter <- as.matrix(Ycenter)
    allY <- rbind(Ycenter, auxY)
    allX <- rbind(Xcenter, auxX)
    nall <- nrow(allY)
    allYmean <- apply(allY, 2, mean)
    allYcenter <- allY - matrix(allYmean, byrow = T, nall, q)
    allXmean <- apply(allX, 2, mean)
    allXcenter <- allX - matrix(allXmean, nall, p, byrow = T)


    intialcvress <- cv.nuclear(Y=allY, X=allX, B=B,L=L, eta=eta, lamseq = lamseq_w,
                   maxiter=maxiter_inital, tol=tol_inital, nfold=nfold)



    pathres <- lapply(1:length(lamseq_delta), function(x){

      Yres <- Ycenter - Xcenter %*% intialcvress$B[-1,]

      ADMM.biascorrection.scad_onestep(Y=Yres, X=Xcenter,What=intialcvress$B[-1,], a=a, Delta=Delta, H=H, Pi=Pi, eta=eta,
                                       lambda=lamseq_delta[x], maxiter=maxiter_biascorrection,
                                       tol=tol_biascorrection, standardize=TRUE) -> fit_delta
      intercept <-  matrix(Ymean, 1, q) -  matrix(Xmean,1,p) %*% ( fit_delta$Delta[-1,]+intialcvress$B[-1,])
      Slope <-  fit_delta$Delta[-1,]+intialcvress$B[-1,]
      Bhat <- rbind(intercept, Slope)
      Bhat
    })
    residual_list <- lapply(pathres, function(x) {
        ((Y - cbind(1, X) %*% x))
    })
    sse <- unlist(lapply(residual_list, function(x) {
        sum(x^2)
    }))
    rankseq <- unlist(lapply(pathres, function(x) {
        sum(corpcor::fast.svd(x[-1, ])$d > 0.00001)
    }))
    df <- rankseq * (p + q) - rankseq^2
    logsse <- log(sse)
    penaltyterm <- log(log(n0))/(n0 * q) * log(p) * df
    BICs <- (logsse + penaltyterm)
    BICs_nozero <- BICs#[rankseq != 0 & rankseq != min(p, q)]
    optwhich_nozero <- which.min(BICs_nozero)
    optwhich_BIC <- which(BICs == BICs_nozero[optwhich_nozero])[1]
    optrank_BIC <- rankseq[optwhich_BIC]
    optB_BIC <- pathres[[optwhich_BIC]]
    return(list(optrank_BIC = optrank_BIC, B = optB_BIC, optwhich_BIC = optwhich_BIC,
        BICs = BICs, rankseq = rankseq))
}

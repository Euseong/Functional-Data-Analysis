library(fda)
library(fda.usc)
library(dplyr)

# To differentiate more than 4 times, functions fdata2fd and fdata.deriv in fda.usc package are modified.
fdata2fd <- function (fdataobj, type.basis = NULL, nbasis = NULL, nderiv = 0, 
                      lambda = NULL, norder=4,...) # norder(order of derivative) is added.
{
  if (is.fdata(fdataobj)) 
    DATA = fdataobj[["data"]]
  else stop("No fdata object")
  np = ncol(DATA)
  tt = fdataobj[["argvals"]]
  rtt <- fdataobj[["rangeval"]]
  if (is.null(lambda)) 
    lambda = 3e-08/diff(rtt)
  if (is.null(nbasis)) {
    nbasis = ifelse(floor(np/3) > floor((np - nderiv - 4)/2), 
                    floor((np - nderiv - 4)/2), floor(np/3))
  }
  as <- list()
  as[[1]] <- rtt
  names(as)[[1]] <- "rangeval"
  as[[2]] <- nbasis
  names(as)[[2]] <- "nbasis"
  C <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("DATA", "type.basis", "nbasis", 
               "nderiv"), names(mf), 0L)
  imetri <- m[2]
  if (imetri == 0) {
    type.basis1 = "bspline"
    a1 <- create.bspline.basis(norder=norder) # B-spline basis is created with norder.
    len.metric <- length(formals(a1))
    vv <- array(0, dim = c(len.metric))
  }
  else {
    a1 <- paste("create.", type.basis, ".basis", 
                sep = "")
    len.metric <- length(formals(a1))
    vv <- array(0, dim = c(len.metric))
  }
  ii <- imetri + 1
  if (C[ii] != "NULL()") {
    ind.m <- 3
    while (C[ii] != "NULL()" && ind.m <= len.metric) {
      aa <- any(names(C) == names(formals(a1))[ind.m])
      if (aa) {
        vv[ind.m] <- which(names(C) == names(formals(a1)[ind.m]))
        ii <- ii + 1
        as[[ind.m]] <- C[[vv[ind.m]]]
        names(as)[[ind.m]] <- names(formals(a1)[ind.m])
      }
      ind.m <- ind.m + 1
    }
  }
  b1.1 <- do.call(a1, as)
  class(DATA) = "matrix"
  fd1.1 <- Data2fd(argvals = tt, y = t(DATA), basisobj = b1.1, 
                   lambda = lambda, ...)
  if (nderiv > 0) 
    fd1.1 = deriv.fd(fd1.1, int2Lfd(nderiv))
  fd1.1
}

fdata.deriv <- function (fdataobj, nderiv = 1, method = "bspline", class.out = "fdata", 
                         nbasis = NULL, norder=4,...) # norder is added.
{
  if (!is.fdata(fdataobj)) 
    fdataobj = fdata(fdataobj)
  nas1 <- is.na(fdataobj)
  if (any(nas1)) 
    stop("fdataobj contain ", sum(nas1), " curves with some NA value \n")
  DATA <- fdataobj[["data"]]
  tt = fdataobj[["argvals"]]
  rtt = fdataobj[["rangeval"]]
  labels = fdataobj[["names"]]
  nr <- nrow(DATA)
  nc <- ncol(DATA)
  if (method == "diff") {
    res = matrix(NA, nrow = nr, ncol = nc)
    for (i in 1:nr) {
      a = diff(DATA[i, ], differences = nderiv)/(tt[2:nc] - 
                                                   tt[1:(nc - 1)])^nderiv
      ab = matrix(NA, ncol = nc, nrow = 2)
      ab[1, 2:nc] = a
      ab[2, 1:(nc - 1)] = a
      res[i, ] = colMeans(ab, na.rm = TRUE)
    }
    labels$main <- paste("d(", labels$main, ",", 
                         nderiv, ")", sep = "")
    res <- fdata(res, tt, rtt, names = labels)
  }
  else {
    if (any(method == c("fmm", "periodic", "natural", 
                        "monoH.FC"))) {
      res = matrix(NA, nrow = nrow(DATA), ncol = ncol(DATA))
      for (i in 1:nrow(DATA)) {
        f1 <- splinefun(x = tt, y = DATA[i, ], method = method)
        res[i, ] = f1(tt, deriv = nderiv)
      }
      labels$main <- paste("d(", labels$main, ",", 
                           nderiv, ")", sep = "")
      res <- fdata(res, tt, rtt, names = labels)
    }
    else {
      if (any(method == c("bspline", "exponential", 
                          "fourier", "monomial", "polynomial"))) {
        # convert fdata object(fda.usc package) to fd object(fda package) with norder.
        res = fdata2fd(fdataobj = fdataobj, type.basis = method,
                       nbasis = nbasis, nderiv = nderiv, norder = norder,...)
        if (class.out == "fdata") {
          ff <- eval.fd(tt, res)
          labels$ylab <- paste("d(", labels$ylab, 
                               ",", nderiv, ")", sep = "")
          res = fdata(t(ff), tt, rtt, names = labels)
        }
      }
    }
  }
  res
}

get.pw.const.beta <- function(fdata., m, method="bspline", n.basis) {
  # returns fitted constant weights for m-th order linear homogeneous differential equation.
  # using point-wise least squares estimation.
  # m : order of the differential equation.
  # method : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # n.basis : number of basis when taking derivative of given data curve used in fdata.deriv function.
  library(fda.usc)
  n <- nrow(fdata.$data) # number of observations.
  x.len <- ncol(fdata.$data) # number of sample points for each curve.
  n.order <- ifelse(m>=4, m+1, 4) # order of splines when taking derivative of given data curve(must be greater than m).
  
  # least squares estimation of data matrix differentiated 0 to m orders.
  y_Zb <- matrix(0, n*x.len, m+1)
  y_Zb[,1] <- as.vector(fdata.deriv(fdata., nderiv=m, norder=n.order, nbasis=n.basis, method=method)$data)
  for (i in 0:(m-1)) { # X
    y_Zb[,i+2] <- -as.vector(fdata.deriv(fdata.$data, nderiv=i, norder=n.order, nbasis=n.basis, method=method)$data)
  }
  y_Zb <- as.data.frame(y_Zb)
  names(y_Zb)[1] <- "y"
  
  pw.ls.fit <- lm(y ~ .-1, data=y_Zb)
  return(coef(pw.ls.fit))
}

get.RSQ <- function(fdata., coefs, method="bspline", n.basis) {
  # returns RSQ of each data curve for fitted weights(coefs).
  d <- length(coefs)
  x.len <- ncol(fdata.$data)
  X <- list(x=fdata.$data)
  for (i in 1:d) {
    if (i >= 4) {n.order <- i + 1} else n.order <- 4
    X[[paste("D", as.character(i), sep="")]] <- fdata.deriv(fdata., nderiv=i, norder=n.order, nbasis=n.basis, method=method)$data
    X[[i]] <- X[[i]] * coefs[i]
  }
  PSSE.0 <- apply(X[[d+1]]^2, 1, sum)
  Lx <- Reduce("+", X)
  PSSE.L <- apply(Lx^2, 1, sum)
  RSQ <- (PSSE.0 - PSSE.L) / PSSE.0
  return(RSQ)
}

get.psi <- function(coefs, x) {
  # compute the solution functions(psi) of given constant weights(coefs).
  # x : values of sample points.
  result <- list()
  solutions <- polyroot(c(coefs, 1))
  for (i in 1:length(solutions)) {
    if (Im(solutions[i]) == 0) { result[[paste("psi", as.character(i), sep=".")]] <- exp(Re(solutions[i])*x) }
    else if (Im(solutions[i]) > 0) {
      result[[paste("psi", as.character(i), sep=".")]] <- exp(Re(solutions[i])*x)*cos( abs(Im(solutions[i])) * x )
    }
    else {
      result[[paste("psi", as.character(i), sep=".")]] <- exp(Re(solutions[i])*x)*sin( abs(Im(solutions[i])) * x )
    }
  }
  return(result)
}

get.fpc.score <- function(fdata, fdata.comp) {
  # returns FPC scores data.frame.
  # fdata.comp : principal component basis(result of create.pc.basis).
  mean.curve <- as.vector(fdata.comp$mean$data)
  epsilons <- list()
  n <- nrow(fdata$data)
  ncomp <- length(fdata.comp$l) # number of FPCs.
  for (m in 1:ncomp) {
    epsilon.im <- c()
    for (i in 1:n) { # calculate the FPC score for each data curve.
      diff <- fdata$data[i,] - mean.curve
      fpc.curve <- fdata.comp$basis$data[m,]
      epsilon.im[i] <- sum(diff*fpc.curve)
    }
    colname <- paste("e", as.character(m), sep="")
    epsilons[[colname]] <- epsilon.im
  }
  result <- as.data.frame(epsilons)
  return(result)
}

RSQ.cv.error <- function(fdata., label, K, max.order, deriv.method="bspline", nbasis) {
  # returns minimum K-fold cross validation classification error via PDA RSQ
  # with its standard deviation and the order of differntial equation.
  # label : label vector aligned with fdata.'s row index.
  # K : fold of cross validation.
  # max.order : maximum order of differential equation. RSQs are obtained from first to max.order-th differential equations.
  # deriv.method : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # nbasis : number of basis when taking derivative of given data curve used in fdata.deriv function.
  n <- nrow(fdata.$data)
  folds <- (1:n)%%K + 1
  folds <- sample(folds, n)
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  cv.errors <- matrix(0, K, max.order)
  for (order in 1:max.order) { # RSQ except i fold data is calculated via foreach function.
    cv.errors[,order] <- foreach (i=1:K, .combine="c", .packages="fda.usc",
                                  .export=c("fdata2fd", "fdata.deriv", "get.RSQ", "get.pw.const.beta")) %dopar% {
                                    train <- (1:n)[folds != i]
                                    # for each fold, train data is separated by its classes or groups
                                    train.class1 <- fdata.
                                    train.class2 <- fdata.
                                    train.class1$data <- fdata.$data[train[which(label[train] == class1)],]
                                    train.class2$data <- fdata.$data[train[which(label[train] == class2)],]
                                    val.fdata <- fdata.
                                    val.fdata$data <- fdata.$data[folds == i,]
                                    
                                    # obtain constant weight for all classes respectively
                                    class1.beta <- get.pw.const.beta(train.class1, m=order, method=deriv.method, n.basis=nbasis)
                                    class2.beta <- get.pw.const.beta(train.class2, m=order, method=deriv.method, n.basis=nbasis)
                                    
                                    # calculate RSQs of each class for test data
                                    class1.RSQ <- get.RSQ(val.fdata, class1.beta, method=deriv.method, n.basis=nbasis)
                                    class2.RSQ <- get.RSQ(val.fdata, class2.beta, method=deriv.method, n.basis=nbasis)
                                    RSQ.df <- data.frame(class1=class1.RSQ, class2=class2.RSQ, y=label)
                                    
                                    # classify validation data as one class with maximum RSQ value
                                    logistic.fit <- glm(y~., data=RSQ.df[train,], family=binomial)
                                    pi.hat <- predict(logistic.fit, newdata=RSQ.df[folds == i,], type="response")
                                    pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1),
                                                             c(class1, class2)))[1:sum(folds == i)]
                                    return(mean(pred != label[folds == i]))
                                  }
    print(paste("order ", " complete", sep=as.character(order)))
  }
  result <- data.frame(mean=apply(cv.errors, 2, mean), standard.error=apply(cv.errors, 2, sd))
  min.error <- result[which.min(result$mean),]
  min.error$order <- which.min(result$mean)
  return(min.error)
}

pda.score.cv.error <- function(fdata., label, K, max.order, deriv.method="bspline", nbasis) {
  # returns minimum K-fold cross validation classification error via PDA scores
  # with its standard deviation and the order of differntial equation.
  # label : label vector aligned with fdata.'s row index.
  # K : fold of cross validation.
  # max.order : maximum order of differential equation. RSQs are obtained from first to max.order-th differential equations.
  # deriv.method : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # nbasis : number of basis when taking derivative of given data curve used in fdata.deriv function.
  n <- nrow(fdata.$data)
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  folds <- (1:n)%%K + 1
  folds <- sample(folds, n)
  cv.errors <- matrix(0, K, max.order)
  for (order in 1:max.order) {
    cv.errors[,order] <- foreach (i=1:K, .combine="c", .packages="fda.usc",
                                  .export=c("fdata2fd", "fdata.deriv", "get.psi", "get.pw.const.beta")) %dopar% {
                                    n.order <- ifelse(order>=4, order+1, 4)
                                    train <- (1:n)[folds != i]
                                    train.fdata <- fdata.
                                    train.fdata$data <- fdata.$data[train,]
                                    coefs <- get.pw.const.beta(train.fdata, m=order, method=deriv.method, n.basis=nbasis)
                                    psi.list <- get.psi(coefs, train.fdata$argvals)
                                    
                                    # for some eigen functions, huge values are indicated as Inf causing error when calculating PDA scores.
                                    # To continue the iteration in spite of Inf value, return Inf for i-th fold classification error.
                                    for (j in 1:order) {
                                      if (sum(is.infinite(psi.list[[j]]))) { # this is done only when if any eigen equation has Inf value.
                                        return(Inf)
                                      }
                                    }
                                    
                                    # calculate PDA scores using eigen functions(psi.list) obtained with train data.
                                    mean.curve <- apply(train.fdata$data, 2, mean)
                                    diff <- as.matrix(fdata.$data) - matrix(mean.curve, nrow=n, ncol=length(mean.curve), byrow=T)
                                    scores <- list()
                                    scores[["y"]] <- label
                                    for (j in 1:length(psi.list)) {
                                      scores[[as.character(j)]] <- diff %*% matrix(psi.list[[j]], ncol=1)
                                    }
                                    scores <- as.data.frame(scores)
                                    
                                    # classify validation data using calculated PDA scores.
                                    logistic.fit <- glm(y~., data=scores[train,], family=binomial)
                                    pi.hat <- predict(logistic.fit, newdata=scores[folds == i,], type="response")
                                    pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1),
                                                             c(class1, class2)))[1:sum(folds == i)]
                                    return(mean(pred != scores$y[folds == i]))
                                  }
    print(paste("order ", " complete", sep=as.character(order)))
  }
  result <- data.frame(mean=apply(cv.errors, 2, mean), standard.error=apply(cv.errors, 2, sd))
  min.error <- result[which.min(result$mean),]
  min.error$order <- which.min(result$mean)
  return(min.error)
}

fpc.score.cv.error <- function(fdata., label, K, max.fpc, lam=1) {
  # returns minimum K-fold cross validation classification error using first to max.fpc-th FPC scores with its standard deviation.
  # label : label vector aligned with fdata.'s row index.
  # K : fold of cross validation.
  # max.fpc : maximum number of funtional principal components.
  # lam : roughness penalty parameter used in create.pc.basis.
  n <- nrow(fdata.$data)
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  folds <- (1:n)%%K + 1
  folds <- sample(folds, n)
  cv.errors <- matrix(0, K, max.fpc)
  for (fpc in 1:max.fpc) {
    cv.errors[,fpc] <- foreach (i=1:K, .combine="c", .packages="fda.usc",
                                .export=c("fdata2fd", "fdata.deriv", "get.fpc.score")) %dopar% {
                                  train <- (1:n)[folds != i]
                                  train.fdata <- fdata.
                                  train.fdata$data <- fdata.$data[train,]
                                  
                                  # obtain pc basis with train data.
                                  train.fpca <- create.pc.basis(train.fdata, l=1:fpc, lambda=lam)
                                  
                                  scores <- get.fpc.score(fdata., train.fpca)
                                  scores$y <- label
                                  
                                  # classify validation data using FPC scores.
                                  logistic.fit <- glm(y~., data=scores[train,], family=binomial)
                                  pi.hat <- predict(logistic.fit, newdata=scores[-train,], type="response")
                                  pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1),
                                                           c(class1, class2)))[1:length(scores$y[-train])]
                                  return(mean(pred != scores$y[-train]))
                                }
    print(paste("max.fpc ", " complete", sep=as.character(fpc)))
  }
  result <- data.frame(mean=apply(cv.errors, 2, mean), standard.error=apply(cv.errors, 2, sd))
  return(result)
}

simul.data <- function(functions, N, mu, error=T, covariance=F) {
  # returns simulation data for given eigen functions.
  # functions : the list of eigen functions which simulation data is based on.
  # N : number of data curves to generate.
  # mu : mean of random weight which is multiplied by each eigen function.
  # error : weight or coefficient of random noise. random noise from N(0,1) is multiplied by this argument, and added to each data curve.
  # covariance : covariance matrix among data points of eigen functions. If False, all data points are assumed to be independent each other.
  #              If covariance matrix is given, errors are generated from multivariate normal distribution, MVN(mu, covariance).
  if (!is.list(functions)) {
    print("functions must be list")
    return()
  }
  n.functions <- length(functions) # number of given eigen functions.
  x.len <- length(functions[[1]]) # number of sample points.
  result <- matrix(0, N, x.len)
  
  if (covariance) {
    if (length(mu) < 2) {print("the number of means must be greater than 1"); return()}
    for (i in 1:n.functions) {
      for (j in 1:N) {
        result[j, ] <- result[j, ] + rnorm(1,mu,1)*functions[[i]] + mvrnorm(1, mu=rep(error, 100), Sigma=covariance)
      }
    }
  return(result)
  }else{
  for (i in 1:n.functions) {
    for (j in 1:N) {
      result[j, ] <- result[j, ] + rnorm(1,mu,1)*functions[[i]] + error*rnorm(x.len, 0, 1)
    }
  }
  return(result)
  }
}

x <- seq(-5, 5, len=200) # whole simulation data uses 200 sample points equally separated within interval (-5, 5).

# source of eigen functions.
{exp0.25 <- exp(0.25*x)
  exp0.5 <-exp(0.5*x)
  exp3 <-exp(3*x)
  exp5 <-exp(5*x)
  sin1 <- sin(1*pi*x)
  cos1 <- cos(1*pi*x)
  sin1.05 <- sin(1.05*pi*x)
  cos1.05 <- cos(1.05*pi*x)
  sin3 <- sin(3*pi*x)
  cos3 <- cos(3*pi*x)
  x.poly1 <- 2*x + x^2 + 3*x^3
  x.poly2 <- 2*x + x^2 - 3*x^3 + x^4
  x.poly3 <- 0.01*x + 0.4*x^2 - 0.003*x^3 + 0.0006*x^4}

# plot simulation data for each class.
{par(mfcol=c(2,2))
  plot.fdata(data1.fd, col=1:5, main=paste("E(a) =", as.character(mu), ", b =", as.character(err)))
  plot.fdata(data2.fd, col=1:5, main=NULL)}
par(mfrow=c(1,1))

library(doParallel)
nc <- detectCores()
registerDoParallel(nc-1)

# generate simulation data set.
nfold <- 8
err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1
{set.seed(100)
  data1 <- simul.data(list(cos1, sin1), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(cos3, sin3), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set1 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label1 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}


# data set2.
err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1
{set.seed(100)
  data1 <- simul.data(list(exp5), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(exp0.5), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set2 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label2 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}



# data set3.
err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1
{set.seed(300)
  data1 <- simul.data(list(exp0.25*cos1, exp0.25*sin1), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(exp0.5*cos3, exp0.5*sin3), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set3 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label3 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}


# data set4.
err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1
{set.seed(700)
  data1 <- simul.data(list(cos1, sin1), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(cos1.05, sin1.05), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set4 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label4 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}


# data set5.
err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1
{set.seed(100)
  data1 <- simul.data(list(x.poly3), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(3*cos3, 3*sin3), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set5 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label5 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}


# PDA after smoothing.
smooth.fdata <- function(fdata., basis, lambda, argvals.=NULL) {
  # returns smoothed functional data.
  # basis : spline basis used for smoothing the functional data.
  # lambda : roughness penalty parameter for smoothing the functional data.
  # argvals. : sample points of out put. If argvals=NULL, sample points is same with fdata..
  library(fda)
  library(fda.usc)
  l <- lambda
  n <- nrow(fdata.$data)
  if (is.null(argvals.)) {
    argvals. <- fdata.$argvals
  }
  basis.mat <- getbasismatrix(argvals., basis) # returns the matrix of values of each basis at sample points.
  fd0 <- fd(basisobj=basis)
  basis.fdPar <- fdPar(fd0, lambda=l)
  fdata.smoothed <- smooth.basis(argvals., t(fdata.$data), basis.fdPar)$fd # coefficients of each basis after smoothing.
  result <- matrix(0, n, length(argvals.))
  
  # convert the coefficients of each basis to the values at sample points.
  for (i in 1:n) {
    result[i,] <- apply(t(basis.mat) * fdata.smoothed$coefs[,i], 2, sum)
  }
  return(fdata(result, argvals=argvals.))
}


basis.length <- round(0.5*length(x)) # number of B-spline basis.
bs.basis <- create.bspline.basis(c(-5,5), basis.length)
bsfd0 <- fd(basisobj=bs.basis)
bsfdPar <- fdPar(bsfd0, lambda=0.5)

f.basis <- create.fourier.basis(c(-5,5), 11) # number of Fourier basis is set to 11.


# plot of smoothed data using B-spline basis.
# data1.fd and data2.fd can be obtained from data set1 to data set5.
sm.basis <- bs.basis
data1.smooth <- smooth.fdata(data1.fd, basis=sm.basis, lambda=0.1)
data2.smooth <- smooth.fdata(data2.fd, basis=sm.basis, lambda=0.1)

{par(mfcol=c(2,2))
  plot.fdata(data1.smooth, col=1:5, main=paste("E(a) =", as.character(mu), "B-spline smoothed"))
  plot.fdata(data2.smooth, col=1:5, main=NULL)}
par(mfrow=c(1,1))

# plot of smoothed data using Fourier basis.
sm.basis <- f.basis
data1.smooth <- smooth.fdata(data1.fd, basis=sm.basis, lambda=0.1)
data2.smooth <- smooth.fdata(data2.fd, basis=sm.basis, lambda=0.1)
{plot.fdata(data1.smooth, col=1:5, main=paste("E(a) =", as.character(mu), "Fourier smoothed"))
  plot.fdata(data2.smooth, col=1:5, main=NULL)}
par(mfrow=c(1,1))

# smoothed data set.
{set.seed(100)
  index <- sample(200, 200)
  data.set.smoothed <- fdata(rbind(data1.smooth$data, data2.smooth$data)[index,], argvals=x)
  label.sm <- as.factor(c(rep("1", 100), rep("2", 100))[index])}

## codes for repeated test
RSQ.cv.order <- function(fdata., label, K, max.order, deriv.method="bspline", s.basis=NULL, nbasis) {
  # returns the order of differntial equation minimizing K-fold cross validation classification error via PDA RSQ.
  # label : label vector aligned with fdata.'s row index.
  # K : fold of cross validation.
  # max.order : maximum order of differential equation. RSQs are obtained from first to max.order-th differential equations.
  # deriv.method : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # nbasis : number of basis when taking derivative of given data curve used in fdata.deriv function.
  n <- nrow(fdata.$data)
  folds <- (1:n)%%K + 1
  folds <- sample(folds, n)
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  cv.errors <- matrix(0, K, max.order)
  for (order in 1:max.order) { # RSQ except i fold data is calculated via foreach function.
    cv.errors[,order] <- foreach (i=1:K, .combine="c", .packages="fda.usc",
                                  .export=c("fdata2fd", "fdata.deriv", "get.RSQ", "get.pw.const.beta")) %dopar% {
                                    train <- (1:n)[folds != i]
                                    # for each fold, train data is separated by its classes or groups
                                    train.class1 <- fdata.
                                    train.class2 <- fdata.
                                    train.class1$data <- fdata.$data[train[which(label[train] == class1)],]
                                    train.class2$data <- fdata.$data[train[which(label[train] == class2)],]
                                    
                                    # obtain constant weight for all classes respectively
                                    class1.beta <- get.pw.const.beta(train.class1, m=order, method=deriv.method, n.basis=nbasis)
                                    class2.beta <- get.pw.const.beta(train.class2, m=order, method=deriv.method, n.basis=nbasis)
                                    
                                    # calculate RSQs of each class for test data
                                    class1.RSQ <- get.RSQ(fdata., class1.beta, method=deriv.method, n.basis=nbasis)
                                    class2.RSQ <- get.RSQ(fdata., class2.beta, method=deriv.method, n.basis=nbasis)
                                    RSQ.df <- data.frame(class1=class1.RSQ, class2=class2.RSQ, y=label)
                                    
                                    # classify validation data through logistic regression model
                                    logistic.fit <- glm(y~., data=RSQ.df[train,], family=binomial)
                                    pi.hat <- predict(logistic.fit, newdata=RSQ.df[folds == i,], type="response")
                                    pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1),
                                                             c(class1, class2)))[1:sum(folds == i)]
                                    return(mean(pred != label[folds == i]))
                                  }
  }
  result <- apply(cv.errors, 2, mean)
  best.order <- which.min(result)
  return(best.order)
}

pda.score.cv.order <- function(fdata., label, K, max.order, deriv.method="bspline", s.basis=NULL, nbasis) {
  # returns the order of differntial equation minimizing K-fold cross validation classification error via PDA scores.
  # label : label vector aligned with fdata.'s row index.
  # K : fold of cross validation.
  # max.order : maximum order of differential equation. RSQs are obtained from first to max.order-th differential equations.
  # deriv.method : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # nbasis : number of basis when taking derivative of given data curve used in fdata.deriv function.
  n <- nrow(fdata.$data)
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  folds <- (1:n)%%K + 1
  folds <- sample(folds, n)
  p <- length(fdata.$argvals)
  cv.errors <- matrix(0, K, max.order)
  for (order in 1:max.order) {
    cv.errors[,order] <- foreach (i=1:K, .combine="c", .packages="fda.usc",
                                  .export=c("fdata2fd", "fdata.deriv", "get.psi", "get.pw.const.beta")) %dopar% {
                                    n.order <- ifelse(order>=4, order+1, 4)
                                    train <- (1:n)[folds != i]
                                    train.fdata <- fdata.
                                    train.fdata$data <- fdata.$data[train,]
                                    coefs <- get.pw.const.beta(train.fdata, m=order, method=deriv.method, n.basis=nbasis)
                                    
                                    # for some eigen functions, huge values are indicated as Inf causing error when calculating PDA scores.
                                    # To continue the iteration in spite of Inf value, return Inf for i-th fold classification error.
                                    psi.list <- get.psi(coefs, train.fdata$argvals)
                                    if (sum(sapply(psi.list, is.infinite))) { # this is done only when if any eigen equation has Inf value.
                                      return(Inf)
                                    }
                                    
                                    # calculate PDA scores using eigen functions(psi.list) obtained with train data.
                                    mean.curve <- apply(train.fdata$data, 2, mean)#[x.i:x.f]
                                    diff <- as.matrix(fdata.$data) - matrix(mean.curve, nrow=n, ncol=length(mean.curve), byrow=T)
                                    scores <- list()
                                    scores[["y"]] <- label
                                    for (j in 1:length(psi.list)) {
                                      scores[[as.character(j)]] <- diff %*% matrix(psi.list[[j]], ncol=1)
                                    }
                                    scores <- as.data.frame(scores)
                                    
                                    # classify validation data using calculated PDA scores.
                                    logistic.fit <- glm(y~., data=scores[train,], family=binomial)
                                    pi.hat <- predict(logistic.fit, newdata=scores[folds == i,], type="response")
                                    pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1),
                                                             c(class1, class2)))[1:sum(folds == i)]
                                    return(mean(pred != scores$y[folds == i]))
                                  }
  }
  result <- apply(cv.errors, 2, mean)
  best.order <- which.min(result)
  return(best.order)
}


RSQ.error <- function(fdata., label, n.repeat, K., max.order., deriv.method.="bspline", nbasis.) {
  # returns mean of repeated classification error via PDA RSQ.
  # label : label vector aligned with fdata.'s row index.
  # n.repeat : number of repetition for calculating classification errors.
  # K. : fold of cross validation.
  # max.order. : maximum order of differential equation. RSQs are obtained from first to max.order-th differential equations.
  # deriv.method. : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # nbasis. : number of basis when taking derivative of given data curve used in fdata.deriv function.
  n <- nrow(fdata.$data)
  n.train <- round(0.8*n)
  n.test <- n - n.train
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  test.results <- foreach (i=1:n.repeat, .packages=c("fda.usc", "doParallel", "dplyr"),
                          .export=c("fdata2fd", "fdata.deriv", "get.RSQ", "get.pw.const.beta", "RSQ.cv.order", "foreach")) %dopar% {
                            set.seed(100*i)
                            train <- sample(1:n, n.train)
                            train.fdata <- fdata.
                            train.class1 <- fdata.
                            train.class2 <- fdata.
                            train.fdata$data <- fdata.$data[train,]
                            train.class1$data <- fdata.$data[train[which(label[train] == class1)],]
                            train.class2$data <- fdata.$data[train[which(label[train] == class2)],]
                            
                            # obtain best order of differential equation that minimizes K-fold CV error via PDA RSQ.
                            best.order <- RSQ.cv.order(train.fdata, label[train], K=K., max.order=max.order.,
                                                       deriv.method=deriv.method., nbasis=nbasis.)
                            
                            # classify test data using RSQ calculated with best.order through logistic regression model.
                            class1.beta <- get.pw.const.beta(train.class1, m=best.order, method=deriv.method., n.basis=nbasis.)
                            class2.beta <- get.pw.const.beta(train.class2, m=best.order, method=deriv.method., n.basis=nbasis.)
                            class1.RSQ <- get.RSQ(fdata., class1.beta, method=deriv.method., n.basis=nbasis.)
                            class2.RSQ <- get.RSQ(fdata., class2.beta, method=deriv.method., n.basis=nbasis.)
                            RSQ.df <- data.frame(class1=class1.RSQ, class2=class2.RSQ, y=label)
                            
                            logistic.fit <- glm(y~., data=RSQ.df[train,], family=binomial)
                            pi.hat <- predict(logistic.fit, newdata=RSQ.df[-train,], type="response")
                            pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1),
                                                     c(class1, class2)))[1:length(pi.hat)]
                            pi.hat.df <- data.frame(pi.hat=logistic.fit$fitted.values, y=logistic.fit$data$y) %>% arrange(pi.hat)
                            posterior <- c(0.5*mean(pi.hat.df$y == "2"))
                            for (j in 1:(nrow(pi.hat.df)-1)) {
                              posterior.j <- 0.5*mean(pi.hat.df$y[1:j] == "1") + 0.5*mean(pi.hat.df$y[-(1:j)] == "2")
                              posterior <- c(posterior, posterior.j)
                            }
                            posterior <- c(posterior, 0.5*mean(pi.hat.df$y == "1"))
                            
                            result.i <- data.frame(test.error=mean(pred != label[-train]), orders=best.order)
                            result.i$bayes.error <- 1 - max(posterior)
                            return(result.i)
                          }
  result <- as.data.frame(apply(t(sapply(test.results, cbind)), 2, unlist))
  return(result)
}

pda.score.error <- function(fdata., label, n.repeat, K., max.order., deriv.method.="bspline", nbasis.) {
  # returns the mean of repeated classification error via PDA scores.
  # label : label vector aligned with fdata.'s row index.
  # n.repeat : number of repetition for calculating classification errors.
  # K. : fold of cross validation.
  # max.order. : maximum order of differential equation. RSQs are obtained from first to max.order-th differential equations.
  # deriv.method. : spline method when taking derivative of given data curve(fdata.) used in fdata.deriv function.
  # nbasis. : number of basis when taking derivative of given data curve used in fdata.deriv function.
  n <- nrow(fdata.$data)
  n.train <- round(0.8*n)
  n.test <- n - n.train
  p <- length(fdata.$argvals)
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  test.errors <- numeric(n.repeat)
  test.results <- foreach (i=1:n.repeat, .packages=c("fda.usc", "doParallel", "dplyr"),
                           .export=c("fdata2fd", "fdata.deriv", "get.psi", "get.pw.const.beta", "pda.score.cv.order", "foreach")) %dopar% {
                             set.seed(100*i)
                             train <- sample(1:n, n.train)
                             train.fdata <- fdata.
                             train.fdata$data <- fdata.$data[train,]
                             
                             # obtain best order of differential equation that minimizes K-fold CV error via PDA scores.
                             best.order <- pda.score.cv.order(train.fdata, label[train], K=K., max.order=max.order.,
                                                              deriv.method=deriv.method., nbasis=nbasis.)
                             coefs <- get.pw.const.beta(train.fdata, m=best.order, method=deriv.method., n.basis=nbasis.)
                             psi.list <- get.psi(coefs, train.fdata$argvals)
                             
                             if (sum(sapply(psi.list, is.infinite))) {
                               return(Inf)
                             }
                             mean.curve <- apply(train.fdata$data, 2, mean)
                             diff <- as.matrix(fdata.$data) - matrix(mean.curve, nrow=n, ncol=length(mean.curve), byrow=T)
                             scores <- list()
                             scores[["y"]] <- label
                             for (j in 1:length(psi.list)) {
                               scores[[as.character(j)]] <- diff %*% matrix(psi.list[[j]], ncol=1)
                             }
                             scores <- as.data.frame(scores)
                             
                             # classify test data using scores calculated with best.order.
                             logistic.fit <- glm(y~., data=scores[train,], family=binomial)
                             pi.hat <- predict(logistic.fit, newdata=scores[-train,], type="response")
                             pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1), c(class1, class2)))[1:n.test]
                             pi.hat.df <- data.frame(pi.hat=logistic.fit$fitted.values, y=logistic.fit$data$y) %>% arrange(pi.hat)
                             
                             posterior <- c(0.5*mean(pi.hat.df$y == "2"))
                             for (j in 1:(nrow(pi.hat.df)-1)) {
                               posterior.j <- 0.5*mean(pi.hat.df$y[1:j] == "1") + 0.5*mean(pi.hat.df$y[-(1:j)] == "2")
                               posterior <- c(posterior, posterior.j)
                             }
                             posterior <- c(posterior, 0.5*mean(pi.hat.df$y == "1"))
                             
                             result.i <- data.frame(test.error=mean(pred != label[-train]), orders=best.order)
                             result.i$bayes.error <- 1 - max(posterior)
                             return(result.i)
                           }
  result <- as.data.frame(apply(t(sapply(test.results, cbind)), 2, unlist))
  return(result)
}

fpc.score.error <- function(fdata., label, n.repeat, max.fpc, lam=1) {
  # returns the mean of repeated classification error using first to max.fpc-th FPC scores.
  # label : label vector aligned with fdata.'s row index.
  # n.repeat : number of repetition for calculating classification errors.
  # max.fpc : maximum number of funtional principal components.
  # lam : roughness penalty parameter used in create.pc.basis.
  n <- nrow(fdata.$data)
  n.train <- round(0.8*n)
  n.test <- n - n.train
  class1 <- levels(label)[1]
  class2 <- levels(label)[2]
  test.errors <- matrix(0, n.repeat, max.fpc)
  bayes.errors <- matrix(0, n.repeat, max.fpc)
  for (fpc in 1:max.fpc) {
    test.results <- foreach (i=1:n.repeat, .packages="fda.usc",
                                  .export=c("fdata2fd", "fdata.deriv", "get.fpc.score")) %dopar% {
                                    set.seed(100*i)
                                    train <- sample(1:n, n.train)
                                    train.fdata <- fdata.
                                    train.fdata$data <- fdata.$data[train,]
                                    train.fpca <- create.pc.basis(train.fdata, l=1:fpc, lambda=lam)
                                    scores <- get.fpc.score(fdata., train.fpca)
                                    scores$y <- label
                                    
                                    logistic.fit <- glm(y~., data=scores[train,], family=binomial)
                                    pi.hat <- predict(logistic.fit, newdata=scores[-train,], type="response")
                                    pred <- as.factor(append(ifelse(pi.hat > 0.5, class2, class1), c(class1, class2)))[1:n.test]
                                    pi.hat.df <- data.frame(pi.hat=logistic.fit$fitted.values, y=logistic.fit$data$y) %>% arrange(pi.hat)
                                    
                                    posterior <- c(0.5*mean(pi.hat.df$y == "2"))
                                    for (j in 1:(nrow(pi.hat.df)-1)) {
                                      posterior.j <- 0.5*mean(pi.hat.df$y[1:j] == "1") + 0.5*mean(pi.hat.df$y[-(1:j)] == "2")
                                      posterior <- c(posterior, posterior.j)
                                    }
                                    posterior <- c(posterior, 0.5*mean(pi.hat.df$y == "1"))
                                    
                                    result.i <- data.frame(test.error=mean(pred != label[-train]), bayes.error=1 - max(posterior))
                                    return(result.i)
                                  }
    test.results <- as.data.frame(apply(t(sapply(test.results, cbind)), 2, unlist))
    test.errors[,fpc] <- test.results$test.error
    bayes.errors[,fpc] <- test.results$bayes.error
    print(paste("max.fpc ", " complete", sep=as.character(fpc)))
  }
  result <- list()
  result$test.errors <- data.frame(mean=apply(test.errors, 2, mean), standard.error=apply(test.errors, 2, sd))
  result$bayes.error <- data.frame(mean=apply(bayes.errors, 2, mean), standard.error=apply(bayes.errors, 2, sd))
  return(result)
}

## simulation data with covariance structure
 # data set 1.
err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1
rho <- 0.1 # used +-0.1 or +-0.4
err.cov <- matrix(rho, 100, 100)
diag(err.cov) <- 1
err.cov <- err*err.cov

{set.seed(100)
  data1 <- simul.data(list(cos1, sin1), 100, rep(mu, 100), error=err, covariance=err.cov)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(cos3, sin3), 100, rep(mu, 100), error=err, covariance=err.cov)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set1 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label1 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}

# data set2.
{set.seed(100)
  data1 <- simul.data(list(exp5), 100, rep(mu, 100), error=err, covariance=err.cov)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(exp0.5), 100, rep(mu, 100), error=err, covariance=err.cov)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set2 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label2 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}



# data set3.
{set.seed(300)
  data1 <- simul.data(list(exp0.25*cos1, exp0.25*sin1), 100, rep(mu, 100), error=err, covariance=err.cov)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(exp0.5*cos3, exp0.5*sin3), 100, rep(mu, 100), error=err, covariance=err.cov)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set3 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label3 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}


# data set4.
{set.seed(700)
  data1 <- simul.data(list(cos1, sin1), 100, rep(mu, 100), error=err, covariance=err.cov)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(cos1.05, sin1.05), 100, rep(mu, 100), error=err, covariance=err.cov)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set4 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label4 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}


# data set5.
{set.seed(100)
  data1 <- simul.data(list(x.poly3), 100, rep(mu, 100), error=err, covariance=err.cov)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(3*cos3, 3*sin3), 100, rep(mu, 100), error=err, covariance=err.cov)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set5 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label5 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}

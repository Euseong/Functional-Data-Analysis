gene <- read.csv("yeast_with_label.csv")
gene <- na.omit(gene[, c(2, 7:24, 80)])

gene$Peak <- as.factor(ifelse(gene$Peak == "G1" | gene$Peak == "M/G1", "G1", "non_G1"))
gene$Peak <- as.factor(ifelse(gene$Peak == "G1", "G1", "non_G1"))

n <- nrow(gene)

library(fda)
bs.basis <- create.bspline.basis(c(0,119), 18)
bsfd0 <- fd(basisobj=bs.basis)
bsfdPar <- fdPar(bsfd0)

cbasis <- create.constant.basis(c(0,119))
cfd0 <- fd(0, cbasis)
cfdPar <- fdPar(cfd0)

mbasis <- create.monomial.basis(c(0,119), 18)
mfd0 <- fd(matrix(0,18,1), mbasis)
mfdPar <- fdPar(mfd0)

train.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[-test, c(2:19)])), bs.basis)$fd
test.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[test, c(2:19)])), bs.basis)$fd
train.g1.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[g1.train, c(2:19)])), bs.basis)$fd
train.g2m.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[g2m.train, c(2:19)])), bs.basis)$fd
train.mg1.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[mg1.train, c(2:19)])), bs.basis)$fd
train.s.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[s.train, c(2:19)])), bs.basis)$fd
train.sg2.fd <- smooth.basis(seq(0,119,len=18), t(as.matrix(gene[sg2.train, c(2:19)])), bs.basis)$fd

fd.Par <- bsfdPar
train.pda <- pda.fd(xfdlist=list(train.fd), bwtlist=list(fd.Par))
train.g1.pda <- pda.fd(xfdlist=list(train.g1.fd), bwtlist=list(fd.Par, fd.Par))
train.g2m.pda <- pda.fd(xfdlist=list(train.g2m.fd), bwtlist=list(fd.Par, fd.Par))
train.s.pda <- pda.fd(xfdlist=list(train.s.fd), bwtlist=list(fd.Par, fd.Par))
train.mg1.pda <- pda.fd(xfdlist=list(train.mg1.fd), bwtlist=list(fd.Par, fd.Par))
train.sg2.pda <- pda.fd(xfdlist=list(train.sg2.fd), bwtlist=list(fd.Par, fd.Par))

library(ggplot2)
g1.weights1 <- data.frame(g1=train.g1.pda$bwtlist[[1]]$fd$coefs)
ggplot(data=g1.weights1, aes(x=1:18, y=g1)) + geom_line() + ylim(-1, 1) +
  geom_line(aes(x=1:18, y=train.g2m.pda$bwtlist[[1]]$fd$coefs, color="G2/M")) +
  geom_line(aes(x=1:18, y=train.mg1.pda$bwtlist[[1]]$fd$coefs, color="M/G1")) +
  geom_line(aes(x=1:18, y=train.s.pda$bwtlist[[1]]$fd$coefs, color="S")) +
  geom_line(aes(x=1:18, y=train.sg2.pda$bwtlist[[1]]$fd$coefs, color="S/G2"))

g1.weights2 <- data.frame(g1=train.g1.pda$bwtlist[[2]]$fd$coefs)
ggplot(data=g1.weights2, aes(x=1:18, y=g1)) + geom_line() + ylim(-1, 1) +
  geom_line(aes(x=1:18, y=train.g2m.pda$bwtlist[[2]]$fd$coefs, color="G2/M")) +
  geom_line(aes(x=1:18, y=train.mg1.pda$bwtlist[[2]]$fd$coefs, color="M/G1")) +
  geom_line(aes(x=1:18, y=train.s.pda$bwtlist[[2]]$fd$coefs, color="S")) +
  geom_line(aes(x=1:18, y=train.sg2.pda$bwtlist[[2]]$fd$coefs, color="S/G2"))
# 모든 Peak에서 weight가 비슷한 양상을 보임

# point-wise least square
get.pw.beta <- function(data, deriv1.data, deriv2.data) {
  time.points <- ncol(data)
  beta.t <- matrix(0, time.points, 2)
  for (i in 1:time.points) {
    y_Zb <- data.frame(y=deriv2.data[,i], x=-data[,i], D1.x=-deriv1.data[,i])
    coefs.i <- coef(lm(y ~ x + D1.x - 1, data=y_Zb))
    beta.t[i,1] <- coefs.i[1]
    beta.t[i,2] <- coefs.i[2]
  }
  return(data.frame(beta0=beta.t[,1], beta1=beta.t[,2]))
}

library(fda.usc)
set.seed(100)
train <- sample(1:n, round(n*0.8))
g1.train <- which(gene$Peak[train] == "G1")
ng1.train <- which(gene$Peak[train] == "non_G1")
g1.fdata <- fdata(gene[g1.train, 2:19])
ng1.fdata <- fdata(gene[ng1.train, 2:19])
gene.fdata <- fdata(gene[, 2:19])

test.fdata <- fdata(gene[-train, 2:19])
test.fdata.d1 <- fdata.deriv(test.fdata, nderiv=1)
test.fdata.d2 <- fdata.deriv(test.fdata, nderiv=2)

# classification with RSQ
pw.g1.beta.t <- get.pw.beta(g1.fdata$data, fdata.deriv(g1.fdata, nderiv=1)$data, fdata.deriv(g1.fdata, nderiv=2)$data)
pw.g1.beta.t <- get.pw.const.beta(g1.fdata, m=2)
pw.ng1.beta.t <- get.pw.beta(ng1.fdata$data, fdata.deriv(ng1.fdata, nderiv=1)$data, fdata.deriv(ng1.fdata, nderiv=2)$data)
pw.ng1.beta.t <- get.pw.const.beta(ng1.fdata, m=2)

g1.beta0 <- data.frame(g1=pw.g1.beta.t[,1])
ggplot(data=g1.beta0, aes(x=seq(0,120,len=18), y=g1)) + geom_line() + ggtitle("pw.g1.beta0(t)") + xlab("time(ms)") +
  geom_line(aes(x=seq(0,120,len=18), y=pw.ng1.beta.t[,1], color="non_G1"))

g1.beta1 <- data.frame(g1=pw.g1.beta.t[,2])
ggplot(data=g1.beta1, aes(x=seq(0,120,len=18), y=g1)) + geom_line() + ggtitle("pw.g1.beta1(t)") + xlab("time(ms)") +
  geom_line(aes(x=seq(0,120,len=18), y=pw.ng1.beta.t[,1], color="non_G1"))

get.RSQ <- function(fdata., coefs) { # for non-constant coefs
  d <- ncol(coefs)
  X <- list(x=fdata.$data)
  for (i in 1:d) {
    if (i >= 4) {n.order <- i + 1} else n.order <- 4
    X[[paste("D", as.character(i), sep="")]] <- fdata.deriv(fdata., nderiv=i, norder=n.order)$data
    X[[i]] <- t( t(X[[i]]) * coefs[,i] )
  }
  PSSE.0 <- apply(X[[d+1]]^2, 1, sum)
  Lx <- Reduce("+", X)
  PSSE.L <- apply(Lx^2, 1, sum)
  RSQ <- (PSSE.0 - PSSE.L) / PSSE.0
  return(RSQ)
}
get.RSQ <- function(fdata., coefs) { # for constant coefs
  d <- length(coefs)
  X <- list(x=fdata.$data)
  for (i in 1:d) {
    if (i >= 4) {n.order <- i + 1} else n.order <- 4
    X[[paste("D", as.character(i), sep="")]] <- fdata.deriv(fdata., nderiv=i, norder=n.order)$data
    X[[i]] <- X[[i]] * coefs[i]
  }
  PSSE.0 <- apply(X[[d+1]]^2, 1, sum)
  Lx <- Reduce("+", X)
  PSSE.L <- apply(Lx^2, 1, sum)
  RSQ <- (PSSE.0 - PSSE.L) / PSSE.0
  return(RSQ)
}
test.RSQ.g1 <- get.RSQ(test.fdata, pw.g1.beta.t)

test.RSQ.ng1 <- get.RSQ(test.fdata, pw.ng1.beta.t)

test.RSQ <- data.frame(g1=test.RSQ.g1, ng1=test.RSQ.ng1)
pred.g <- as.factor(ifelse(apply(test.RSQ, 1, which.min) == 1, "G1", "non_G1"))
table(pred.g, gene$Peak[-train])
mean(pred.g != gene$Peak[-train]) # 0.3852459(non-constant beta) / 0.352459(constant beta)


# classification with pda scores
get.pw.const.beta <- function(fdata., m) {
  library(fda.usc)
  n <- nrow(fdata.$data)
  time.points <- ncol(fdata.$data)
  n.order <- ifelse(m>4, m+1, 4)
  
  y_Zb <- matrix(0, n*time.points, m+1)
  y_Zb[,1] <- as.vector(fdata.deriv(fdata.$data, nderiv=m, norder=n.order)$data) # y
  for (i in 0:(m-1)) { # X
    y_Zb[,i+2] <- -as.vector(fdata.deriv(fdata.$data, nderiv=i, norder=n.order)$data)
  }
  y_Zb <- as.data.frame(y_Zb)
  names(y_Zb)[1] <- "y"
  
  pw.ls.fit <- lm(y ~ .-1, data=y_Zb)
  return(coef(pw.ls.fit))
}

get.psi <- function(coefs, x) { # compute the solution functions(psi) of given coefficient
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

pda.cv.error <- function(cv.fdata, label, max.order, k, seed=100, x.range=c(-1,1)) {
  library(fda.usc)
  set.seed(seed)
  n <- nrow(cv.fdata$data)
  folds <- sample((1:n)%%k + 1, n)
  
  p <- ncol(cv.fdata$data)
  times <- seq(x.range[1], x.range[2], len=p)
  cv.fdata$argvals <- times
  cv.fdata$rangeval <- x.range
  
  cv.error <- matrix(0, max.order, k)
  for (order in 1:max.order) {
    cv.error[order,] <- foreach (i=1:k, .combine="c", .packages="fda.usc",
                                 .export=c("fdata2fd", "fdata.deriv", "get.pw.const.beta", "get.psi")) %dopar% {
                                   if (order >= 4) {n.order <- order + 1} else n.order <- 4
                                   train.i <- cv.fdata
                                   train.i$data <- cv.fdata$data[folds != i,]
                                   coefs <- get.pw.const.beta(train.i, m=order)
                                   psi.list <- get.psi(coefs, times)
                                   mean.curve <- apply(cv.fdata$data[folds != i,], 2, mean)
                                   diff <- as.matrix(cv.fdata$data) - matrix(mean.curve, nrow=n, ncol=p, byrow=T)
                                   scores <- list()
                                   scores[["y"]] <- label
                                   for (j in 1:length(psi.list)) {
                                     scores[[as.character(j)]] <- diff %*% matrix(psi.list[[j]], ncol=1)
                                   }
                                   scores <- as.data.frame(scores)
                                   
                                   logistic.fit <- glm(y~., data=scores[folds != i,], family=binomial)
                                   pi.hat <- predict(logistic.fit, newdata=scores[folds == i,], type="response")
                                   pred <- as.factor(ifelse(pi.hat > 0.5, levels(label)[2], levels(label)[1]))
                                   return(mean(pred != scores$y[folds == i]))
                                 }
    print(paste(as.character(order), "order complete", sep="-"))
  }
  cv.error <- as.data.frame(cv.error)
  cv.error$order <- 1:max.order
  colnames(cv.error)[1:k] <- as.character(1:k)
  return(cv.error)
}

library(doParallel)
nc <- detectCores()
registerDoParallel(nc-1)

train.fdata <- gene.fdata
train.fdata$data <- gene.fdata$data[train,]

gene.pda <- pda.cv.error(train.fdata, gene$Peak[train], k=5, max.order=3, x.range=c(0, 0.12))

gene.pda$mean <- apply(gene.pda[,1:5], 1, mean)

gene.pda$S.E <- apply(gene.pda[,1:5], 1, sd)

gene.pda


# pda after smoothing
gene <- read.csv("yeast_with_label.csv")
gene <- na.omit(gene[, c(2, 7:24, 80)])
n <- nrow(gene)
set.seed(100)
train <- sample(1:n, round(n*0.8))
gene$Peak <- as.factor(ifelse(gene$Peak == "G1" | gene$Peak == "M/G1", "G1", "non_G1"))
gene$Peak <- as.factor(ifelse(gene$Peak == "G1", "G1", "non_G1"))

train.fdata <- fdata(gene[train, 2:19], argvals=seq(0, 119, len=18))
gene.fdata <- fdata(gene[, 2:19], argvals=seq(0, 119, len=18))

lambdas <- 10^(-7:3)

gcv.search <- function(lam, fdata., basis){
  gcv <- smooth.basisPar(argvals=fdata.$argvals, y=t(fdata.$data), fdobj=basis, lambda=lam)$gcv
  return(data.frame(lambda=lam, n.basis=basis$nbasis, gcv=mean(gcv)))
}

bs.params <- data.frame(lambda=0, n.basis=0, gcv=0)
for (i in seq(10, 100, 10)) {
  bs.basis <- create.bspline.basis(c(0,119), i)
  gcv.results <- data.frame(t(sapply(lambdas, gcv.search, fdata.=gene.fdata, basis=bs.basis)))
  bs.params <- rbind(bs.params, gcv.results[which.min(gcv.results$gcv),])
}
bs.params <- bs.params[-1,]
bs.params

lambdas2 <- seq(50, 500, 10)
bs.params <- data.frame(lambda=0, n.basis=0, gcv=0)
for (i in seq(5, 17, 2)) {
  bs.basis <- create.bspline.basis(c(0,119), i)
  gcv.results <- data.frame(t(sapply(lambdas2, gcv.search, fdata.=gene.fdata, basis=bs.basis)))
  bs.params <- rbind(bs.params, gcv.results[which.min(gcv.results$gcv),])
}
bs.params <- bs.params[-1,]
bs.params

bs.basis <- create.bspline.basis(c(0,119), 11)
gene.bs.smooth <- smooth.fdata(gene.fdata, basis=bs.basis, lambda=150)

f.params <- data.frame(lambda=0, n.basis=0, gcv=0)
for (i in seq(5, 17, 2)) {
  f.basis <- create.fourier.basis(c(0,119), i)
  gcv.results <- data.frame(t(sapply(lambdas, gcv.search, fdata.=gene.fdata, basis=f.basis)))
  f.params <- rbind(f.params, gcv.results[which.min(gcv.results$gcv),])
}
f.params <- f.params[-1,]
f.params

lambdas2 <- seq(500, 10000, 500)
f.params <- data.frame(lambda=0, n.basis=0, gcv=0)
for (i in c(3, 5, 7)) {
  f.basis <- create.fourier.basis(c(0,119), i)
  gcv.results <- data.frame(t(sapply(lambdas2, gcv.search, fdata.=gene.fdata, basis=f.basis)))
  f.params <- rbind(f.params, gcv.results[which.min(gcv.results$gcv),])
}
f.params <- f.params[-1,]
f.params

f.basis <- create.fourier.basis(c(0,119), 15)
gene.f.smooth <- smooth.fdata(gene.fdata, basis=f.basis, lambda=10000)

plot.fdata(gene.fdata, col=ifelse(gene$Peak=="G1", 1, 8), main=NULL, xlab="min", ylab="gene expression")
plot.fdata(gene.bs.smooth, col=ifelse(gene$Peak=="G1", 1, 8), main=NULL, xlab="min", ylab="gene expression")
plot.fdata(gene.f.smooth, col=ifelse(gene$Peak=="G1", 1, 8), main=NULL, xlab="min", ylab="gene expression")


library(xtable)
nfold <- 10
repeat.n <- 1
label.sm <- gene$Peak
data.set.smoothed <- gene.fdata
data.set.smoothed <- gene.bs.smooth
data.set.smoothed <- gene.f.smooth

{set.seed(100)
set.sm.RSQ.error <- RSQ.cv.error(data.set.smoothed, label.sm, K=10, max.order=11, deriv.method="bspline", nbasis=12)
set.sm.pda.error <- pda.score.cv.error(data.set.smoothed, label.sm, K=10, max.order=11, deriv.method="bspline", nbasis=12)
set.sm.fpc.error <- fpc.score.cv.error(data.set.smoothed, label.sm, K=10, max.fpc=11)
set.bsm.result0.rt <- cbind(rbind(rbind(set.sm.RSQ.error[,-3], set.sm.pda.error[,-3]), set.sm.fpc.error),
                            max.fpc=c(set.sm.RSQ.error[,3], set.sm.pda.error[,3], 1:11), stringsAsFactors=F)}

{set.seed(100)
set.sm.RSQ.error <- RSQ.cv.error(data.set.smoothed, label.sm, K=10, max.order=11, deriv.method="fourier", nbasis=12)
set.sm.pda.error <- pda.score.cv.error(data.set.smoothed, label.sm, K=10, max.order=11, deriv.method="fourier", nbasis=12)
set.fsm.result0.rt <- cbind(rbind(rbind(set.sm.RSQ.error[,-3], set.sm.pda.error[,-3]), set.sm.fpc.error),
                            max.fpc=c(set.sm.RSQ.error[,3], set.sm.pda.error[,3], 1:11), stringsAsFactors=F)
set.sm.result0.rt <- cbind(paste(round(set.bsm.result0.rt[,1], 2), " (", round(set.bsm.result0.rt[,2], 2), ")", sep=""), set.bsm.result0.rt[,3],
                           paste(round(set.fsm.result0.rt[,1], 2), " (", round(set.fsm.result0.rt[,2], 2), ")", sep=""), set.fsm.result0.rt[,3])
row.names(set.sm.result0.rt) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA", 11), as.character(1:11), sep=""))
xtable(set.sm.result0.rt)}

# repeated classification error of smoothed data set with deriv.method="bspline".
{set.seed(100)
set.sm.RSQ.error <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=10, deriv.method.="bspline", nbasis.=50)
set.sm.pda.error <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=10, deriv.method.="bspline", nbasis.=50)
set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
set.bsm.result0.rt <- cbind(rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error), max.fpc=c("-", "-", 1:5), stringsAsFactors=F)}

# repeated classification error of smoothed data set with deriv.method="fourier".
{set.seed(100)
set.sm.RSQ.error <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=10, deriv.method.="fourier", nbasis.=15)
set.sm.pda.error <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=10, deriv.method.="fourier", nbasis.=15)
set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
set.fsm.result0.rt <- cbind(rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error), max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
set.sm.result0.rt <- cbind(paste(round(set.bsm.result0.rt[,1], 2), " (", round(set.bsm.result0.rt[,2], 2), ")", sep=""), set.bsm.result0.rt[,3],
                           paste(round(set.fsm.result0.rt[,1], 2), " (", round(set.fsm.result0.rt[,2], 2), ")", sep=""), set.fsm.result0.rt[,3])
row.names(set.sm.result0.rt) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA", 5), as.character(1:5), sep=""))
xtable(set.sm.result0.rt)}



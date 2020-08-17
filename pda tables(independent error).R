lambdas <- 10^(-5:3)
library(doParallel)
library(xtable)
nc <- detectCores()
registerDoParallel(nc-1)
nfold <- 8
repeat.n <- 100
basis.length <- round(0.5*length(x)) # number of B-spline basis.
bs.basis <- create.bspline.basis(c(-5,5), basis.length)
f.basis <- create.fourier.basis(c(-5,5), 11)
gcv.mat <- matrix(0, 1, 2)
gcv.search <- function(lam, fdata., basis){
  gcv <- smooth.basisPar(argvals=fdata.$argvals, y=t(fdata.$data), fdobj=basis, lambda=lam)$gcv
  return(c(lam, basis$nbasis, mean(gcv)))
}

err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1

mu <- 1

err <- 1
mu <- 0

mu <- 1#
## data set 4
{set.seed(700)
  data1 <- simul.data(list(cos1, sin1), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(cos1.05, sin1.05), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set4 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label4 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}

# repeat this for 4 scenarios (Fourier basis derivative)
{
gcv.result1 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=bs.basis)))
colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
gcv.result2 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=bs.basis)))
colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
}
{
  data1.smooth <- smooth.fdata(data1.fd, basis=bs.basis, lambda=gcv.lambda1)
  data2.smooth <- smooth.fdata(data2.fd, basis=bs.basis, lambda=gcv.lambda2)
  set.seed(100)
  index <- sample(200, 200)
  data.set.smoothed <- fdata(rbind(data1.smooth$data, data2.smooth$data)[index,], argvals=x)
  label.sm <- as.factor(c(rep("1", 100), rep("2", 100))[index])
  set.seed(100)
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="fourier", nbasis.=11)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="fourier", nbasis.=11)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.bsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.bsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
}
{
  gcv.result1 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=f.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=f.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
}
{
  data1.smooth <- smooth.fdata(data1.fd, basis=f.basis, lambda=gcv.lambda1)
  data2.smooth <- smooth.fdata(data2.fd, basis=f.basis, lambda=gcv.lambda2)
  set.seed(100)
  index <- sample(200, 200)
  data.set.smoothed <- fdata(rbind(data1.smooth$data, data2.smooth$data)[index,], argvals=x)
  label.sm <- as.factor(c(rep("1", 100), rep("2", 100))[index])
  set.seed(100)
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="fourier", nbasis.=11)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="fourier", nbasis.=11)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.fsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.fsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  
  set.sm.result0.rt <- cbind(paste(round(set.bsm.result0.rt[,1], 2), " (", round(set.bsm.result0.rt[,2], 2), ")", sep=""), #set.bsm.result0.rt[,3],
                             paste(round(set.fsm.result0.rt[,1], 2), " (", round(set.fsm.result0.rt[,2], 2), ")", sep=""))#, set.fsm.result0.rt[,3])
  set.sm.result0.bayes <- cbind(paste(round(set.bsm.result0.bayes[,1], 2), " (", round(set.bsm.result0.bayes[,2], 2), ")", sep=""),
                             paste(round(set.fsm.result0.bayes[,1], 2), " (", round(set.fsm.result0.bayes[,2], 2), ")", sep=""))
  
  row.names(set.sm.result0.rt) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  row.names(set.sm.result0.bayes) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  xtable(set.sm.result0.rt)
  xtable(set.sm.result0.bayes)
}





err <- 0.1 # used 0.1 or 1
mu <- 0    # used 0 or 1

mu <- 1

err <- 1
mu <- 0

mu <- 1#

## data set 5
{set.seed(100)
  data1 <- simul.data(list(x.poly3), 100, mu, error=err)
  data1.fd <- fdata(data1, argvals=x)
  data2 <- simul.data(list(3*cos3, 3*sin3), 100, mu, error=err)
  data2.fd <- fdata(data2, argvals=x)
  index <- sample(200, 200)
  data.set5 <- fdata(rbind(data1, data2)[index,], argvals=x)
  label5 <- as.factor(c(rep("1", 100), rep("2", 100))[index])}

# repeat this for 4 scenarios (B-spline basis derivative)
{
  gcv.result1 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=bs.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=bs.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
}
{
  data1.smooth <- smooth.fdata(data1.fd, basis=bs.basis, lambda=gcv.lambda1)
  data2.smooth <- smooth.fdata(data2.fd, basis=bs.basis, lambda=gcv.lambda2)
  set.seed(100)
  index <- sample(200, 200)
  data.set.smoothed <- fdata(rbind(data1.smooth$data, data2.smooth$data)[index,], argvals=x)
  label.sm <- as.factor(c(rep("1", 100), rep("2", 100))[index])
  set.seed(100)
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="bspline", nbasis.=100)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="bspline", nbasis.=100)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.bsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.bsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
}
{
  gcv.result1 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=f.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=f.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
}
{
  data1.smooth <- smooth.fdata(data1.fd, basis=f.basis, lambda=gcv.lambda1)
  data2.smooth <- smooth.fdata(data2.fd, basis=f.basis, lambda=gcv.lambda2)
  set.seed(100)
  index <- sample(200, 200)
  data.set.smoothed <- fdata(rbind(data1.smooth$data, data2.smooth$data)[index,], argvals=x)
  label.sm <- as.factor(c(rep("1", 100), rep("2", 100))[index])
  set.seed(100)
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="bspline", nbasis.=100)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="bspline", nbasis.=100)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.fsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.fsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  
  set.sm.result0.rt <- cbind(paste(round(set.bsm.result0.rt[,1], 2), " (", round(set.bsm.result0.rt[,2], 2), ")", sep=""), #set.bsm.result0.rt[,3],
                             paste(round(set.fsm.result0.rt[,1], 2), " (", round(set.fsm.result0.rt[,2], 2), ")", sep=""))#, set.fsm.result0.rt[,3])
  set.sm.result0.bayes <- cbind(paste(round(set.bsm.result0.bayes[,1], 2), " (", round(set.bsm.result0.bayes[,2], 2), ")", sep=""),
                                paste(round(set.fsm.result0.bayes[,1], 2), " (", round(set.fsm.result0.bayes[,2], 2), ")", sep=""))
  
  row.names(set.sm.result0.rt) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  row.names(set.sm.result0.bayes) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  xtable(set.sm.result0.rt)
  xtable(set.sm.result0.bayes)
}
# write.csv(gcv.mat, "lambda.csv")

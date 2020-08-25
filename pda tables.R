library(doParallel)
library(xtable)
library(stringr)
working.dir <- "C:/Users/gkfrj/Documents/R"
err <- c(0.1, 1) # used 0.1 or 1
mu <- c(0, 1) # used 0 or 1
rho <- c(0.1, 0.4, 0.7) # used +-0.1 or +-

a <- merge(err, mu)
b <- merge(err, rho)
colnames(a) <- c("err", "mu")
colnames(b) <- c("err", "rho")
cov.params <- merge(a, b, by=err)
simul.params <- rbind(cbind(a, rho=rep(0, nrow(a))), cov.params, cov.params)


nc <- detectCores()
registerDoParallel(nc-1)
lambdas <- 10^(-5:3)
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


## 
# data set 3
data.set.list <- data.set3.list
data.set.num <- 3
results.list <- list()
path <- paste0(working.dir, "/data set", data.set.num, " tables.txt")
scenario.num <- 1
start.time <- Sys.time()
for (data.set.i in data.set.list) {
  data1.fd <- data.set.i$data1.fd
  data2.fd <- data.set.i$data2.fd
  params.i <- simul.params[scenario.num,]
  colnames(params.i) <- NULL
  err <- params.i[1]; mu <- params.i[2]; rho <- params.i[3]
  # repeat 28 scenarios (Fourier basis derivative)
  
  gcv.result1 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=bs.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=bs.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
  
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
  
  BS.smoothing.results <- list(RSQ=set.sm.RSQ.result, PDA.score=set.sm.pda.result)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.bsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.bsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  
  gcv.result1 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=f.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=f.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
  
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
  
  Fourier.smoothing.results <- list(RSQ=set.sm.RSQ.result, PDA.score=set.sm.pda.result)
  list.name <- names(data.set.list)[scenario.num]
  results.list[[list.name]] <- list(Bs.smoothing=BS.smoothing.results, F.smoothing=Fourier.smoothing.results)
    
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.fsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.fsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  
  set.bsm.result0 <- cbind(paste(round(set.bsm.result0.rt[,1], 2), " (", round(set.bsm.result0.rt[,2], 2), ")", sep=""),
                           paste(round(set.bsm.result0.bayes[,1], 2), " (", round(set.bsm.result0.bayes[,2], 2), ")", sep=""))
  set.fsm.result0 <- cbind(paste(round(set.fsm.result0.rt[,1], 2), " (", round(set.fsm.result0.rt[,2], 2), ")", sep=""),
                           paste(round(set.fsm.result0.bayes[,1], 2), " (", round(set.fsm.result0.bayes[,2], 2), ")", sep=""))
  
  row.names(set.bsm.result0) <- c("& PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  # row.names(set.sm.result0.bayes) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  # test.error <- xtable(set.sm.result0.rt,
  #                      caption=str_interp("data set ${data.set.num} ${scenario.num}th test error (mu=${mu}, error=${err}, rho=${rho})"))
  # bayes.error <- xtable(set.sm.result0.bayes,
  #                       caption=str_interp("data set ${data.set.num} ${scenario.num}th bayes error (mu=${mu}, error=${err}, rho=${rho})"))
  data.set.i.results <- xtable(cbind(set.bsm.result0, set.fsm.result0),
                               caption=str_interp("data set ${data.set.num} ${scenario.num}th test error (mu=${mu}, error=${err}, rho=${rho})"))
  
  is.append <- scenario.num > 1
  print(data.set.i.results, file = path, compress = FALSE, append=is.append)
  # print(bayes.error, file = path, compress = FALSE, append=T)
  print(paste0(scenario.num, "th scenario over"))
  scenario.num <- scenario.num + 1
  print(Sys.time())
}
saveRDS(results.list, file=paste0(working.dir, "/data set", data.set.num, " results.RData", sep=""))
paste("start :", start.time, "end :", Sys.time())


## data set 5
data.set.list <- data.set5.list
data.set.num <- 5
path <- paste0(working.dir, "/data set", data.set.num, " tables.txt")
scenario.num <- 1
start.time <- Sys.time()
for (data.set.i in data.set.list) {
  data1.fd <- data.set.i$data1.fd
  data2.fd <- data.set.i$data2.fd
  params.i <- simul.params[scenario.num,]
  colnames(params.i) <- NULL
  err <- params.i[1]; mu <- params.i[2]; rho <- params.i[3]
  # repeat 28 scenarios (Fourier basis derivative)
  
  gcv.result1 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=bs.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=bs.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
  
  data1.smooth <- smooth.fdata(data1.fd, basis=bs.basis, lambda=gcv.lambda1)
  data2.smooth <- smooth.fdata(data2.fd, basis=bs.basis, lambda=gcv.lambda2)
  set.seed(100)
  index <- sample(200, 200)
  data.set.smoothed <- fdata(rbind(data1.smooth$data, data2.smooth$data)[index,], argvals=x)
  label.sm <- as.factor(c(rep("1", 100), rep("2", 100))[index])
  set.seed(100)
  get.pw.const.beta
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="bspline", nbasis.=100)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=4, deriv.method.="bspline", nbasis.=100)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  BS.smoothing.results <- list(RSQ=set.sm.RSQ.result, PDA.score=set.sm.pda.result)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.bsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.bsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  
  gcv.result1 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data1.fd, basis=f.basis)))
  colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
  gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]
  gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=data2.fd, basis=f.basis)))
  colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
  gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]
  gcv.mat <- rbind(gcv.mat, c(gcv.lambda1, gcv.lambda2))
  
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
  
  Fourier.smoothing.results <- list(RSQ=set.sm.RSQ.result, PDA.score=set.sm.pda.result)
  list.name <- names(data.set.list)[scenario.num]
  results.list[[list.name]] <- list(Bs.smoothing=BS.smoothing.results, F.smoothing=Fourier.smoothing.results)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  set.fsm.result0.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  set.fsm.result0.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  
  set.bsm.result0 <- cbind(paste(round(set.bsm.result0.rt[,1], 2), " (", round(set.bsm.result0.rt[,2], 2), ")", sep=""),
                           paste(round(set.bsm.result0.bayes[,1], 2), " (", round(set.bsm.result0.bayes[,2], 2), ")", sep=""))
  set.fsm.result0 <- cbind(paste(round(set.fsm.result0.rt[,1], 2), " (", round(set.fsm.result0.rt[,2], 2), ")", sep=""),
                           paste(round(set.fsm.result0.bayes[,1], 2), " (", round(set.fsm.result0.bayes[,2], 2), ")", sep=""))
  
  row.names(set.bsm.result0) <- c("& PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  # row.names(set.sm.result0.bayes) <- c("PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  # test.error <- xtable(set.sm.result0.rt,
  #                      caption=str_interp("data set ${data.set.num} ${scenario.num}th test error (mu=${mu}, error=${err}, rho=${rho})"))
  # bayes.error <- xtable(set.sm.result0.bayes,
  #                       caption=str_interp("data set ${data.set.num} ${scenario.num}th bayes error (mu=${mu}, error=${err}, rho=${rho})"))
  data.set.i.results <- xtable(cbind(set.bsm.result0, set.fsm.result0),
                               caption=str_interp("data set ${data.set.num} ${scenario.num}th test error (mu=${mu}, error=${err}, rho=${rho})"))
  
  is.append <- scenario.num > 1
  print(data.set.i.results, file = path, compress = FALSE, append=is.append)
  # print(bayes.error, file = path, compress = FALSE, append=T)
  print(paste0(scenario.num, "th scenario over"))
  scenario.num <- scenario.num + 1
}
saveRDS(results.list, file=paste0(working.dir, "/data set", data.set.num, " results.RData", sep=""))
paste("start :", start.time, "end :", Sys.time())

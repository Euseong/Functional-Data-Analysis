library(fda.usc)
library(xtable)
library(doParallel)

# aneurisk data
{
  patient <- read.csv("AneuRisk65/Patients.txt", sep=" ")
  patient$type[which(patient$type == "N")] <- "L"
  patient$type <- as.factor(as.character(patient$type))
  FKS1 <- read.csv("AneuRisk65/Rawdata_FKS_1.txt", sep=" ")
  s <- which(FKS1$Curv_Abscissa < -0.74  & FKS1$Curv_Abscissa >= -32.9)
  FKS <- as.data.frame(matrix(FKS1$MISR[s[1]:(s[1]+340)], 1))
}

x.length <- numeric(65)
x.length[1] <- length(s)
for (i in 2:65) {
  file <- paste("AneuRisk65/Rawdata_FKS_", as.character(i), sep="")
  file <- paste(file, ".txt", sep="")
  FKS.i <- read.csv(file, sep=" ")
  s <- which(FKS.i$Curv_Abscissa < -0.74  & FKS.i$Curv_Abscissa >= -32.9)
  FKS <- rbind(FKS, matrix(FKS.i$MISR[s[1]:(s[1]+340)], 1))
  # x.length[i] <- length(s)
}
min(x.length) # 341


FKS.fd <- fdata(FKS)
FKS.fd$argvals <- seq(-32.9, -0.74, len=341)
FKS.fd$rangeval <- c(-32.9, -0.74)
label <- patient$type
pdf(file = paste0('pda simulation results/image/aneurisk', '.pdf'), width=6, height=5)
plot.fdata(FKS.fd, col=ifelse(label=="L", 1, 8), main=NULL, xlab="abscissa", ylab="ICA radius")
legend("topright", c("G1", "non_g1"), col=c("black", "gray"), lty=c(1,1), cex=0.7)

basis.length <- round(0.5*length(FKS.fd$argvals)) # number of B-spline basis.
bs.basis <- create.bspline.basis(FKS.fd$rangeval, basis.length)
f.basis <- create.fourier.basis(FKS.fd$rangeval, 41)
gcv.mat <- matrix(0, 1, 2)
gcv.search <- function(lam, fdata., basis){
  gcv <- smooth.basisPar(argvals=fdata.$argvals, y=t(fdata.$data), fdobj=basis, lambda=lam)$gcv
  return(c(lam, basis$nbasis, mean(gcv)))
}

lambdas <- 10^(-5:3)
gcv.result1 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=FKS.fd, basis=bs.basis)))
colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]

gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=FKS.fd, basis=f.basis)))
colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]

FKS.bs.smooth <- smooth.fdata(FKS.fd, basis=bs.basis, lambda=gcv.lambda1)
FKS.f.smooth <- smooth.fdata(FKS.fd, basis=f.basis, lambda=gcv.lambda2)
plot.fdata(FKS.bs.smooth, col=1)
plot.fdata(FKS.f.smooth, col=1)

nc <- detectCores()
registerDoParallel(nc-1)
nfold <- 10
repeat.n <- 100
max.order.in <- 10
label.sm <- patient$type


sm.basis <- c("without smooth", "B-spline", "Fourier")
results.list <- list()
i <- 1
for (data.set.smoothed in list(FKS.fd, FKS.bs.smooth, FKS.f.smooth)) {
  path <- paste0("AneuRisk2 ", sm.basis[i], " tables.txt")
  set.seed(100)
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=max.order.in, deriv.method.="bspline", nbasis.=41)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=max.order.in, deriv.method.="bspline", nbasis.=41)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  results.list[[sm.basis[i]]] <- list(RSQ=set.sm.RSQ.result, PDA.score=set.sm.pda.result)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  sm.i.result.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  sm.i.result.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  sm.i.results <- cbind(paste(round(sm.i.result.rt[,1], 2), " (", round(sm.i.result.rt[,2], 2), ")", sep=""),
                        paste(round(sm.i.result.bayes[,1], 2), " (", round(sm.i.result.bayes[,2], 2), ")", sep=""))
  
  row.names(sm.i.results) <- c("& PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  results.xtable <- xtable(sm.i.results, caption=str_interp("${sm.basis[i]} test error"))
  
  is.append <- i > 1
  print(results.xtable, file = path, compress = FALSE, append=is.append)
  print(paste0(i, "th smoothing over"))
  i <- i + 1
  print(Sys.time())
}
saveRDS(results.list, file=paste0("AneuRisk2", " classification results.RData", sep=""))

# growth data
str(growth)
growth.fdata <- fdata(t(cbind(growth$hgtm, growth$hgtf)), argvals=growth$age)
set.seed(100)
new.index <- sample(93, 93)
growth.fdata$data <- growth.fdata$data[new.index,]
gender <- as.factor(c(rep("boy", 39), rep("girl", 54))[new.index])

pdf(file = paste0('pda simulation results/image/growth', '.pdf'), width=6, height=5)
plot.fdata(growth.fdata, col=ifelse(gender=="boy", 1, 3), main=NULL, xlab="age", ylab="heights")
legend("topleft", c("boy", "girl"), col=c("black", "green"), lty=c(1,1), cex=1)

basis.length <- round(0.5*length(growth.fdata$argvals)) # number of B-spline basis.
bs.basis <- create.bspline.basis(growth.fdata$rangeval, basis.length)
f.basis <- create.fourier.basis(growth.fdata$rangeval, 15)
gcv.mat <- matrix(0, 1, 2)
gcv.search <- function(lam, fdata., basis){
  gcv <- smooth.basisPar(argvals=fdata.$argvals, y=t(fdata.$data), fdobj=basis, lambda=lam)$gcv
  return(c(lam, basis$nbasis, mean(gcv)))
}

lambdas <- 10^(-5:3)
gcv.result1 <- as.data.frame(t(sapply(lambdas, gcv.search, fdata.=growth.fdata, basis=bs.basis)))
colnames(gcv.result1) <- c("lambda", "n.basis", "gcv")
gcv.lambda1 <- gcv.result1$lambda[which.min(gcv.result1$gcv)]

gcv.result2 <- data.frame(t(sapply(lambdas, gcv.search, fdata.=growth.fdata, basis=f.basis)))
colnames(gcv.result2) <- c("lambda", "n.basis", "gcv")
gcv.lambda2 <- gcv.result2$lambda[which.min(gcv.result2$gcv)]

growth.bs.smooth <- smooth.fdata(growth.fdata, basis=bs.basis, lambda=gcv.lambda1)
growth.f.smooth <- smooth.fdata(growth.fdata, basis=f.basis, lambda=gcv.lambda2)
plot.fdata(growth.bs.smooth, col=1)
plot.fdata(growth.f.smooth, col=1)


nc <- detectCores()
registerDoParallel(nc-1)
nfold <- 10
repeat.n <- 100
max.order.in <- 10
label.sm <- gender


sm.basis <- c("without smooth", "B-spline", "Fourier")
results.list <- list()
i <- 1
for (data.set.smoothed in list(growth.fdata, growth.bs.smooth, growth.f.smooth)) {
  path <- paste0("growth ", sm.basis[i], " tables.txt")
  set.seed(100)
  set.sm.RSQ.result <- RSQ.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=max.order.in, deriv.method.="fourier", nbasis.=15)
  set.sm.pda.result <- pda.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, K.=nfold, max.order.=max.order.in, deriv.method.="fourier", nbasis.=15)
  set.sm.fpc.error <- fpc.score.error(data.set.smoothed, label.sm, n.repeat=repeat.n, max.fpc=5)
  
  results.list[[sm.basis[i]]] <- list(RSQ=set.sm.RSQ.result, PDA.score=set.sm.pda.result)
  
  set.sm.RSQ.error <- data.frame(mean=mean(set.sm.RSQ.result$test.error), standard.error=sd(set.sm.RSQ.result$test.error))
  set.sm.RSQ.bayes.error <- data.frame(mean=mean(set.sm.RSQ.result$bayes.error), standard.error=sd(set.sm.RSQ.result$bayes.error))
  set.sm.pda.error <- data.frame(mean=mean(set.sm.pda.result$test.error), standard.error=sd(set.sm.pda.result$test.error))
  set.sm.pda.bayes.error <- data.frame(mean=mean(set.sm.pda.result$bayes.error), standard.error=sd(set.sm.pda.result$bayes.error))
  
  sm.i.result.rt <- rbind(rbind(set.sm.RSQ.error, set.sm.pda.error), set.sm.fpc.error$test.errors)#, max.fpc=c("-", "-", 1:5), stringsAsFactors=F)
  sm.i.result.bayes <- rbind(rbind(set.sm.RSQ.bayes.error, set.sm.pda.bayes.error), set.sm.fpc.error$bayes.error)
  sm.i.results <- cbind(paste(round(sm.i.result.rt[,1], 2), " (", round(sm.i.result.rt[,2], 2), ")", sep=""),
                        paste(round(sm.i.result.bayes[,1], 2), " (", round(sm.i.result.bayes[,2], 2), ")", sep=""))
  
  row.names(sm.i.results) <- c("& PDA(RSQ)", "& & PDA(scores)", paste(rep("& & FPCA (num of PCs=", 5), as.character(1:5), ")", sep=""))
  results.xtable <- xtable(sm.i.results, caption=str_interp("${sm.basis[i]} test error"))
  
  is.append <- i > 1
  print(results.xtable, file = path, compress = FALSE, append=is.append)
  print(paste0(i, "th smoothing over"))
  i <- i + 1
  print(Sys.time())
}
saveRDS(results.list, file=paste0("growth", " classification results.RData", sep=""))

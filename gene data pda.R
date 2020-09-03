# pda after smoothing
gene <- read.csv("yeast_with_label.csv")
gene <- na.omit(gene[, c(2, 7:24, 80)])

gene$Peak <- as.factor(ifelse(gene$Peak == "G1" | gene$Peak == "M/G1", "G1", "non_G1"))
gene$Peak <- as.factor(ifelse(gene$Peak == "G1", "G1", "non_G1"))

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
repeat.n <- 100
max.order.in <- 10
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

# raw.RSQ.order1 <- RSQ.cv.order(train.fdata, label=gene$Peak[train], K=nfold, max.order=10, deriv.method="bspline", s.basis=bs.basis, nbasis=65)
# bs.RSQ.order1 <- RSQ.cv.order(train.fdata, label=gene$Peak[train], K=nfold, max.order=10, deriv.method="bspline", s.basis=bs.basis, nbasis=65)
# f.RSQ.order1 <- RSQ.cv.order(train.fdata, label=gene$Peak[train], K=nfold, max.order=10, deriv.method="bspline", s.basis=bs.basis, nbasis=65)
# 
# gene.betas <- get.pw.const.beta(fdata., m=, method="bspline", n.basis=65)
# 
# raw.RSQ.order2 <- RSQ.cv.order(train.fdata, label=gene$Peak[train], K=nfold, max.order=10, deriv.method="fourier", s.basis=f.basis, nbasis=15)
# bs.RSQ.order2 <- RSQ.cv.order(train.fdata, label=gene$Peak[train], K=nfold, max.order=10, deriv.method="fourier", s.basis=f.basis, nbasis=15)
# f.RSQ.order2 <- RSQ.cv.order(train.fdata, label=gene$Peak[train], K=nfold, max.order=10, deriv.method="fourier", s.basis=f.basis, nbasis=15)
# 
# gene.betas2 <- get.pw.const.beta(fdata., m=, method="fourier", n.basis=65)

# repeated classification error of smoothed data set with deriv.method="bspline".
library(doParallel)
nc <- detectCores()
registerDoParallel(nc-1)



sm.basis <- c("without smooth", "B-spline", "Fourier")
results.list <- list()
i <- 1
for (data.set.smoothed in list(gene.fdata, gene.bs.smooth, gene.f.smooth)) {
  path <- paste0("gene ", sm.basis[i], " tables.txt")
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
saveRDS(results.list, file=paste0("gene", " classification results.RData", sep=""))


gene <- read.csv("yeast_with_label.csv")
gene <- na.omit(gene[, c(2, 7:24, 80)])
n <- nrow(gene)
table(gene$Peak)
index <- as.numeric(rownames(gene))
set.seed(100)
g1.test <- sample(which(gene$Peak == "G1"), 10)
g2m.test <- sample(which(gene$Peak == "G2/M"), 10)
mg1.test <- sample(which(gene$Peak == "M/G1"), 10)
s.test <- sample(which(gene$Peak == "S"), 10)
sg2.test <- sample(which(gene$Peak == "S/G2"), 10)
test <- c(g1.test, g2m.test, mg1.test, s.test, sg2.test)

g1.train <- which(gene$Peak == "G1")[!(which(gene$Peak == "G1") %in% g1.test)]
g2m.train <- which(gene$Peak == "G2/M")[!(which(gene$Peak == "G2/M") %in% g2m.test)]
mg1.train <- which(gene$Peak == "M/G1")[!(which(gene$Peak == "M/G1") %in% mg1.test)]
s.train <- which(gene$Peak == "S")[!(which(gene$Peak == "S") %in% s.test)]
sg2.train <- which(gene$Peak == "S/G2")[!(which(gene$Peak == "S/G2") %in% sg2.test)]

library(fda.usc)
train.g1.fdata <- fdata(gene[c(g1.train, mg1.train), 2:19])
train.g2m.fdata <- fdata(gene[g2m.train, 2:19])
train.mg1.fdata <- fdata(gene[mg1.train, 2:19])
train.s.fdata <- fdata(gene[s.train, 2:19])
train.sg2.fdata <- fdata(gene[sg2.train, 2:19])
train.Xg1.fdata <- fdata(gene[-c(g1.train, mg1.train), 2:19])

par(mfrow=c(3,2))
plot.fdata(train.g1.fdata, col=1)
plot.fdata(train.g2m.fdata, col=1)
plot.fdata(train.mg1.fdata, col=1)
plot.fdata(train.s.fdata, col=1)
plot.fdata(train.sg2.fdata, col=1)
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot.fdata(train.g1.fdata, col=1)
plot.fdata(train.Xg1.fdata, col=1)
par(mfrow=c(1,1))

train.g1.pca <- create.pc.basis(train.g1.fdata, 1:18)
train.g2m.pca <- create.pc.basis(train.g2m.fdata, 1:18)
train.mg1.pca <- create.pc.basis(train.mg1.fdata, 1:18)
train.s.pca <- create.pc.basis(train.s.fdata, 1:18)
train.sg2.pca <- create.pc.basis(train.sg2.fdata, 1:18)

pc.curve <- function(m) {
  g1.pc1 <- data.frame(g1=train.g1.pca$basis$data[m,])
  time <- seq(0,119,len=18)
  ggplot(data=g1.pc1, aes(x=time, y=g1)) + geom_line() + ylim(-1,1) + ggtitle(as.character(m)) +
    geom_line(aes(x=time, y=train.g2m.pca$basis$data[m,], color="G2/M")) +
    geom_line(aes(x=time, y=train.mg1.pca$basis$data[m,], color="M/G1")) +
    geom_line(aes(x=time, y=train.s.pca$basis$data[m,], color="S")) +
    geom_line(aes(x=time, y=train.sg2.pca$basis$data[m,], color="S/G2"))
}
pc.curve(1)#
pc.curve(2)#
pc.curve(3)#
pc.curve(4)
pc.curve(5)
pc.curve(6)
pc.curve(7)
pc.curve(8)
pc.curve(9)
pc.curve(10)
pc.curve(11)
pc.curve(12)
pc.curve(13)
pc.curve(14)
pc.curve(15)
pc.curve(16)
pc.curve(17)
pc.curve(18)#

set.seed(100)
train <- sample(n, round(n*0.8))
train.label <- gene[train, c(1,20)]
test.label <- gene[-train, c(1,20)]
train.fdata <- fdata(gene[train, 2:19])
test.fdata <- fdata(gene[-train, 2:19])

get.epsilon <- function(fdata, fdata.comp) {
  # mean.curve <- as.vector(fdata.comp$mean$data)
  epsilons <- list()
  n <- nrow(fdata$data)
  ncomp <- length(fdata.comp$l)
  for (m in 1:ncomp) { # 각 FPC 마다 epsilon.im 계산
    epsilon.im <- c()
    for (i in 1:n) { # 1개의 train set은 총 100개의 데이터
      # diff <- fdata$data[i,] - mean.curve
      fpc.curve <- fdata.comp$basis$data[m,]
      epsilon.im[i] <- sum(fdata$data[i,]*fpc.curve)
    }
    colname <- paste("e", as.character(m), sep="")
    epsilons[[colname]] <- epsilon.im
  }
  result <- as.data.frame(epsilons)
  return(result)
}

curve2epsilon <- function(train, test, ncomp=5) {
  library(fda.usc)
  train.pca <- create.pc.basis(train, l=1:ncomp)
  
  train.e <- get.epsilon(train, train.pca)
  test.e <- get.epsilon(test, train.pca)
  
  result <- list(train=train.e, test=test.e)
  return(result)
}

train.g1.fpc <- get.epsilon(train.fdata, fdata.comp=train.g1.pca)[,c(1:3,18)]
train.g2m.fpc <- get.epsilon(train.fdata, fdata.comp=train.g2m.pca)[,c(1:3,18)]
train.mg1.fpc <- get.epsilon(train.fdata, fdata.comp=train.mg1.pca)[,c(1:3,18)]
train.s.fpc <- get.epsilon(train.fdata, fdata.comp=train.s.pca)[,c(1:3,18)]
train.sg2.fpc <- get.epsilon(train.fdata, fdata.comp=train.sg2.pca)[,c(1:3,18)]

fpc.each.type <- function(m) {
  return(data.frame(g1=train.g1.fpc[,m], g2m=train.g2m.fpc[,m], mg1=train.mg1.fpc[,m],
                    s=train.s.fpc[,m], sg2=train.sg2.fpc[,m]))
}
train.fpc1 <- fpc.each.type(1)
fpc1.pred <- levels(gene$Peak)[apply(abs(train.fpc1), 1, which.max)]
train.fpc2 <- fpc.each.type(2)
fpc2.pred <- levels(gene$Peak)[apply(abs(train.fpc2), 1, which.max)]
train.fpc3 <- fpc.each.type(3)
fpc3.pred <- levels(gene$Peak)[apply(abs(train.fpc3), 1, which.max)]
train.fpc18 <- fpc.each.type(4)
fpc18.pred <- levels(gene$Peak)[apply(abs(train.fpc18), 1, which.max)]
table(fpc1.pred, train.label$Peak)
table(fpc2.pred, train.label$Peak)
table(fpc3.pred, train.label$Peak)
table(fpc18.pred, train.label$Peak)

library(nnet)
train.fpcs <- data.frame(m1=fpc1.pred, m2=fpc2.pred, m3=fpc3.pred, m18=fpc18.pred, Peak=train.label$Peak)
train.multi <- multinom(Peak~., data=train.fpcs)

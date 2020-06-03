library(fda)
library(fda.usc)
library(dplyr)
library(ggplot2)
# library(data.table)
set.seed(100)
cell <- read.delim("mbc_9_12_3273__CDCDATA.txt")
label <- read.csv("CellCycle98.csv")
cell$gene <- as.character(cell$gene)

str(label)
head(label, 20)
label$ORF <- as.character(label$ORF)
names(label)[2] <- "gene"
nrow(label)
length(unique(label$ORF))
label$Peak[label$ORF == as.character(cell$gene[1])]
a <- 0
for (i in 1:10) {
  a <- a + sum(as.character(cell$gene) == label$ORF[i])
}
labeled.data <- cell[which(cell$gene %in% label$gene),]
nrow(labeled.data)
labeled.data <- merge(labeled.data, label)
summary(labeled.data$Peak)
write.csv(labeled.data, "yeast_with_label.csv")

gene <- cell[, c(1, 6:23)] # 0~120min, 18 points
gene <- na.omit(gene)
# n <- nrow(gene)
# index <- sample(1:n, size=2000)

# name <- as.character(gene[,1])
# name.part <- paste(substr(name, 1, 3), substr(name, nchar(name), nchar(name)))
# paste(substr(name, 1, 3), substr(name, nchar(name), nchar(name))) %>% unique() -> name.unique
# 
# train.index <- c()
# for (gene.name in name.unique) {
#   name.index <- which(name.part==gene.name)
#   name.sample <- sample(name.index, size=2)
#   train.index <- c(train.index, name.sample)
# }

# used fda package
# train <- gene[, -1]
# train <- t(train)
# basis <- create.bspline.basis(c(0,120), nbasis=6)
# train.smooth <- smooth.basis(argvals=seq(0, 119, len=18), y=train[,train.index], fdParobj=basis)
# # train.fd <- Data2fd(argvals=seq(0, 1, len=18), y=train[,-1])
# train.fdPar <- fdPar(basis, lambda=0)
# train.pca <- pca.fd(train.smooth, nharm=5, harmfdPar=train.fdPar)
# train.pca$varprop
# sum(train.pca$varprop)

# used fda.usc package
train <- gene[, -1] # remove the column of gene's names
train.fdata <- fdata(train)
# train.fd <- fdata2fd(train.fdata, nbasis=5)
# train.pca <- create.pc.basis(train.fdata, l=1:5, lambda=3)
train.pca <- create.pc.basis(train.fdata, l=1:5, lambda=1)
train.pca$basis$data <- -train.pca$basis$data

# b.bases <- create.fdata.basis(train.fdata, l=1:10)


plot(train.pca$basis, ylim=c(-0.6, 0.6))
par(mfrow=c(2,3))
for (i in 1:5) {
  plot(train.pca$basis[i,], main=paste("FPC", i), ylim=c(-0.6, 0.6))
}
par(mfrow=c(1,1))
plot(train.pca, ylim=c(-0.4, 0.2))



# epsilon 1~5, group label을 성분으로 하는 100개의 datasets(trian 100, test 100) 생성
means <- seq(0.6, 0.2, -0.1)
variances <- c(2.6950, 0.8850, 0.1957, 0.1266, 0.1079)
random.sets <- list()
for (i in 1:100) {
  e1 <- append( rnorm(100, means[1], sqrt(variances[1])), rnorm(100, -means[1], sqrt(variances[1])) )
  e2 <- append( rnorm(100, means[2], sqrt(variances[2])), rnorm(100, -means[2], sqrt(variances[2])) )
  e3 <- append( rnorm(100, means[3], sqrt(variances[3])), rnorm(100, -means[3], sqrt(variances[3])) )
  e4 <- append( rnorm(100, means[4], sqrt(variances[4])), rnorm(100, -means[4], sqrt(variances[4])) )
  e5 <- append( rnorm(100, means[5], sqrt(variances[5])), rnorm(100, -means[5], sqrt(variances[5])) )
  group <- as.factor(c( rep("g1", 100), rep("g0", 100) ))
  random.set <- list(e1=e1, e2=e2, e3=e3, e4=e4, e5=e5, group=group)
  set.i <- paste("set", as.character(i), sep="")
  random.sets[[set.i]] <- random.set
}
random.sets$set1$e1[1:10]
append( rnorm(10, means[1], sqrt(variances[1])), rnorm(10, -means[1], sqrt(variances[1])) )
e1.new <- as.data.frame((random.sets$set1$e1[c(1:5)] %o% train.pca$basis$data[1,]))
e2.new <- as.data.frame(t(random.sets$set1$e2[c(1:5)] %o% train.pca$basis$data[2,]))

# mean curve에 mth eigenfunction_m과 epsilon_m을 곱한 값을 더한 X samples 생성
X.curves <- list()
mean.curve <- as.vector(train.pca$mean$data)
pca.coefs <- train.pca$basis$data
i <- 1
for (set.i in random.sets) {
  sum.e <- mean.curve # 1 datase합t당 PC X e들의 합
  for (j in 1:5) {
    ej <- paste("e", as.character(j), sep="")
    ej.new <- as.data.frame( (set.i[[ej]] %o% pca.coefs[j,]) )
    sum.e <- ej.new + sum.e
  }
  X.curves[[names(random.sets)[i]]] <- sum.e
  i <- i + 1
}

# train.smooth$coefs에 각 세트별 X.coefs를 대입하여 해당 세트의 X 샘플 생성(train, test 나눠서)
train.index <- c(1:50, 101:150)
group <- random.sets$set1$group
train.group <- data.frame(group=group[train.index])
test.group <- data.frame(group=group[-train.index])
train1 <- fdata(X.curves$set1)[train.index]
test1 <- fdata(X.curves$set1)[-train.index]

train.test.split <- function(dataset) { # dataset = X.curves$set.i
  train.index <- c(1:50, 101:150)
  library(fda.usc)
  dataset.fd <- fdata(dataset)
  df <- list(train=dataset.fd[train.index], test=dataset.fd[-train.index])
  return(df)
}

get.epsilon <- function(data, fdata.comp) { # fdata.comp : create.pc.basis 결과
  mean.curve <- as.vector(fdata.comp$mean$data)
  epsilons <- list()
  n <- nrow(data)
  ncomp <- length(fdata.comp$l)
  for (m in 1:ncomp) { # 각 FPC 마다 epsilon.im 계산
    epsilon.im <- c()
    for (i in 1:n) { # 1개의 train set은 총 100개의 데이터
      diff <- data[i,] - mean.curve
      fpc.curve <- fdata.comp$basis$data[m,]
      epsilon.im[i] <- sum(diff*fpc.curve)
    }
    colname <- paste("e", as.character(m), sep="")
    epsilons[[colname]] <- epsilon.im
  }
  result <- as.data.frame(epsilons)
  return(result)
}

curve2epsilon <- function(train.test.df, ncomp=5) {
  train <- train.test.df$train
  test <- train.test.df$test$data
  
  library(fda.usc)
  train.pca <- create.pc.basis(train, l=1:ncomp, lambda=1)
  
  train.e <- get.epsilon(train$data, train.pca)
  test.e <- get.epsilon(test, train.pca)
  
  result <- list(train=train.e, test=test.e)
  return(result)
}

error.fpca <- data.frame(mean.g1=0, mean.g2=0, mean.overall=0, std.g1=0, std.g2=0, std.overall=0)
M <- 5 # 사용할 FPC의 개수
for (m in 1:M) {
  error.rate <- list(overall=c(), g1=c(), g2=c())
  i <- 1
  for (curve in X.curves) { # X.curves set1 ~ set100까지의 error.rate를 계산하여 평균냄
    set.i <- train.test.split(curve)
    set.e <- curve2epsilon(set.i)
    group <- as.factor(c( rep("g1", 50), rep("g2", 50) ))
    set.e$train$group <- group
    set.e$test$group <- group
    model <- glm(group ~ ., data=set.e$train[,c(1:m, M+1)], family = binomial)
    pi.hat <- predict(model, set.e$test[,c(1:m, M+1)], type="response")
    pred <- ifelse(pi.hat > 0.5, 1, 0)
    c.tab <- table(set.e$test$group, pred)
    error.rate$overall[i] <- 1 - (sum(diag(c.tab)) / sum(c.tab))
    error.rate$g1[i] <- 1 - ((c.tab[1,1]) / sum(c.tab[1,]))
    error.rate$g2[i] <- 1 - ((c.tab[2,2]) / sum(c.tab[2,]))
    i <- i + 1
  }
  fpca.mean <- lapply(error.rate, mean) # 0.0748
  fpca.std <- lapply(error.rate, sd) # 0.0311
  fpca.means <- c(fpca.mean$g1, fpca.mean$g2, fpca.mean$overall)
  fpca.stds <- c(fpca.std$g1, fpca.std$g2, fpca.std$overall)
  error.fpca <- rbind(error.fpca, c(fpca.means, fpca.stds))
}
error.fpca

ncomp <- 5
train1.pca <- create.pc.basis(train1, l=1:ncomp, lambda=1)
train1.mean <- as.vector(train1.pca$mean$data)
train1.pca$basis$data[1,] + train1.mean
epsilons <- list()
for (m in 1:ncomp) { # 각 FPC 마다 epsilon.im 계산
  epsilon.im <- c()
  for (i in 1:100) { # 1개의 train set은 총 100개의 데이터
    diff <- train1.pca$fdataobj$data[i,] - train1.mean
    fpc.curve <- train1.pca$basis$data[m,]
    epsilon.im[i] <- sum(diff*fpc.curve)
  }
  colname <- paste("e", as.character(m), sep="")
  epsilons[[colname]] <- epsilon.im
}
epsilons <- as.data.frame(epsilons)

set1 <- train.test.split(X.curves$set1)
set1.e <- curve2epsilon(set1)
group <- as.factor(c( rep("g1", 50), rep("g2", 50) ))
set1.e$train$group <- group
set1.e$test$group <- group
model <- glm(group ~ ., data=set1.e$train, family = binomial)
pi.hat <- predict(model, set1.e$test[,1:5], type="response")
pred <- ifelse(pi.hat > 0.5, 1, 0)
c.tab <- table(set1.e$test$group, pred)
error.rate <- 1 - sum(diag(c.tab)) / length(pi.hat)


# B-spline method
# X.curves의 set 당 감마 data.frame 계산
get.gamma <- function(data, mean.curve, n.basis=n.basis) { # fdata.comp : create.pc.basis 결과
  library(fda.usc)
  data$data <- t(apply(data$data, 1, function(x) x - mean.curve)) # X.curve에 mean.curve를 빼줌
  data.fd <- fdata2fd(data, nbasis=n.basis)
  result <- as.data.frame(t(data.fd$coefs))
  return(result)
}

curve2gamma <- function(train.test.df, n.basis=5) {
  train <- train.test.df$train
  test <- train.test.df$test
  
  library(fda.usc)
  train.pca <- create.pc.basis(train, lambda=1)
  mean.curve <- as.vector(train.pca$mean$data)
  
  train.gamma <- get.gamma(train, mean.curve, n.basis=n.basis)
  test.gamma <- get.gamma(test, mean.curve, n.basis=n.basis)
  
  result <- list(train=train.gamma, test=test.gamma)
  return(result)
}

train1.fd <- fdata2fd(train1, nbasis=5)
fdata2fd(test1, type.basis = train1.fd$basis)
train1.pca <- create.pc.basis(train1, lambda=1)
mean.curve <- as.vector(train1.pca$mean$data)
diff <- t(apply(train1$data[1:5,], 1, function(x) x - mean.curve))
train1$data[1:5,] - diff



Q <- 5 # 사용할 b-spline의 개수
error.b <- data.frame(mean.g1=0, mean.g2=0, mean.overall=0, std.g1=0, std.g2=0, std.overall=0)
for (q in 1:Q) {
  i <- 1
  error.rate.b <- list(overall=c(), g1=c(), g2=c())
  for (curve in X.curves) { # X.curves set1 ~ set100까지의 error.rate를 계산하여 평균냄
    set.i <- train.test.split(curve)
    set.g <- curve2gamma.b(set.i)
    group <- as.factor(c( rep("g1", 50), rep("g2", 50) ))
    set.g$train$group <- group
    set.g$test$group <- group
    model <- glm(group ~ ., data=set.g$train[,c(1:q, Q+1)], family = binomial)
    pi.hat <- predict(model, set.g$test[,c(1:q, Q+1)], type="response")
    pred <- ifelse(pi.hat > 0.5, 1, 0)
    c.tab <- table(set.g$test$group, pred)
    error.rate.b$overall[i] <- 1 - (sum(diag(c.tab)) / sum(c.tab))
    error.rate.b$g1[i] <- 1 - ((c.tab[1,1]) / sum(c.tab[1,]))
    error.rate.b$g2[i] <- 1 - ((c.tab[2,2]) / sum(c.tab[2,]))
    i <- i + 1
  }
  b.mean <- lapply(error.rate.b, mean) # 0.0781
  b.std <- lapply(error.rate.b, sd) # 0.0308
  b.means <- c(b.mean$g1, b.mean$g2, b.mean$overall)
  b.stds <- c(b.std$g1, b.std$g2, b.std$overall)
  error.b <- rbind(error.b, c(b.means, b.stds))
}
100*error.b
100*error.fpca


train.fd <- fdata2fd(train.fdata, nbasis=10)
mean.coefs <- c()
for (i in 1:10) {
  if (i %% 2 == 0) {
    mean.coefs[i] <- mean(train.fd$coefs[i,])
  } else {mean.coefs[i] <- 0}
}
eval.fd(mean.coefs, train.fd)

curve2gamma.b <- function(train.test.df, n.basis=5) {
  train <- train.test.df$train
  test <- train.test.df$test
  
  library(fda.usc)
  train.fd <- fdata2fd(train, nbasis=10)
  mean.curve <- c()
  for (i in seq(2, n.basis, 2)) {
    mean.curve[i/2] <- mean(train.fd$coefs[i,]))
  }
  mean.curve <- as.vector(mean.curve)
  
  train.gamma <- get.gamma(train, mean.curve, n.basis=n.basis)
  test.gamma <- get.gamma(test, mean.curve, n.basis=n.basis)
  
  result <- list(train=train.gamma, test=test.gamma)
  return(result)
}

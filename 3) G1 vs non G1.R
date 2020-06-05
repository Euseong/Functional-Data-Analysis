gene <- read.csv("yeast_with_label.csv")
gene <- na.omit(gene[, c(2, 7:24, 80)])
gene$Peak <- as.factor(ifelse(gene$Peak == "G1" | gene$Peak == "M/G1", "G1", "non_G1"))
table(gene$Peak)

n.gene <- nrow(gene)
set.seed(100)
train <- sample(n.gene, round(n.gene*0.8))

library(fda.usc)
str(gene)
train.label <- gene[train, c(1,20)]
test.label <- gene[-train, c(1,20)]
train.fdata <- fdata(gene[train, 2:19])
test.fdata <- fdata(gene[-train, 2:19])

get.epsilon <- function(fdata, fdata.comp) {
  mean.curve <- as.vector(fdata.comp$mean$data)
  epsilons <- list()
  n <- nrow(fdata$data)
  ncomp <- length(fdata.comp$l)
  for (m in 1:ncomp) { # 각 FPC 마다 epsilon.im 계산
    epsilon.im <- c()
    for (i in 1:n) { # 1개의 train set은 총 100개의 데이터
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

curve2epsilon <- function(train, test, ncomp=5, l=0) {
  library(fda.usc)
  train.pca <- create.pc.basis(train, l=1:ncomp, lambda=l)
  
  train.e <- get.epsilon(train, train.pca)
  test.e <- get.epsilon(test, train.pca)
  
  result <- list(train=train.e, test=test.e)
  return(result)
}

lambda <- 0
m <- 18
gene.epsilon <- curve2epsilon(train.fdata, test.fdata, ncomp=m, l=lambda)
gene.epsilon$train <- cbind(gene.epsilon$train, train.label)
gene.epsilon$test <- cbind(gene.epsilon$test, test.label)

train.deriv <- fdata.deriv(train.fdata, nderiv=1)
test.deriv <- fdata.deriv(test.fdata, nderiv=1)

gene.deriv.epsilon <- curve2epsilon(train.deriv, test.deriv, ncomp=m, l=lambda)
gene.epsilon$train <- cbind(gene.epsilon$train, gene.deriv.epsilon$train)
gene.epsilon$test <- cbind(gene.epsilon$test, gene.deriv.epsilon$test)

train.deriv2 <- fdata.deriv(train.fdata, nderiv=2)
test.deriv2 <- fdata.deriv(test.fdata, nderiv=2)

gene.deriv2.epsilon <- curve2epsilon(train.deriv2, test.deriv2, ncomp=m)
gene.epsilon$train <- cbind(gene.epsilon$train, gene.deriv2.epsilon$train)
gene.epsilon$test <- cbind(gene.epsilon$test, gene.deriv2.epsilon$test)

plot(train.fdata$data[2,], col=1, type="l", xlab="time points", ylab="gene expression", ylim=c(-2, 2), main="Original Data")
plot(train.deriv$data[2,], col=1, type="l", xlab="time points", ylab="", ylim=c(-1, 1), main="Differentiated Data")

train.fpc0 <- create.pc.basis(train.fdata, l=1:5)
train.fpc1 <- create.pc.basis(train.fdata, l=1:5, lambda=1)
train.pc <- prcomp(train.fdata$data)

plot(1:18, train.fpc0$basis$data[1,], col=1, type="p", xlab="time points", ylab="gene expression", main="FPC lambda=0")
plot(1:18, train.fpc1$basis$data[1,], col=1, type="p", xlab="time points", ylab="gene expression", main="FPC1 lambda=1", ylim=c(-0.5, 0.3))
plot(1:18, train.pc$rotation[,1], col=1, xlab="time points", ylab="gene expression", main="PC1", ylim=c(-0.5, 0.3))

K <- 10
index <- (1:nrow(gene.epsilon$train))%%K + 1
set.seed(100)
folds <- sample(index, nrow(gene.epsilon$train))
cv.error.for.m <- numeric(m)
# curve only
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- glm(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20)], family="binomial")
    pi.hat <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20)], type="response")
    pred <- as.factor(ifelse(pi.hat > 0.5, "non_G1", "G1"))
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m # 0.09417517 for m=5
logistic.fit1 <- glm(Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], family="binomial")
logistic.pred1 <- as.factor(ifelse(predict(logistic.fit1, gene.epsilon$test[, c(1:5, 20)], type="response") > 0.5, "non_G1", "G1"))
mean(logistic.pred1 != gene.epsilon$test[,20]) # 0.09836066
table(logistic.pred1, gene.epsilon$test[,20])


# with first derivative
set.seed(100)
cv.error.for.m <- matrix(0,m,m)
cv.error.sd <- matrix(0,m,m)
for (i in 1:m) {
  for (k in 1:m) {
    pred.error <- numeric(K)
    for (j in 1:K) {
      model <- glm(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20, 21:(20+k))], family="binomial")
      pi.hat <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20, 21:(20+k))], type="response")
      pred <- as.factor(ifelse(pi.hat > 0.5, "non_G1", "G1"))
      pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
    } 
    cv.error.for.m[k,i] <- mean(pred.error)
    cv.error.sd[k,i] <- sd(pred.error)
  }
}
cv.error.for.m
min(cv.error.for.m) # 0.08596939
m1.m0 <- as.data.frame(which(cv.error.for.m == min(cv.error.for.m), arr.ind=T))
m1.m0

test.error <- numeric(nrow(m1.m0))
for (i in 1:nrow(m1.m0)) {
  logistic.fit2 <- glm(Peak ~ ., data=gene.epsilon$train[, c(1:(0+m1.m0$col[i]), 20, 21:(20+m1.m0$row[i]))], family="binomial")
  pi.hat <- predict(logistic.fit2, gene.epsilon$test[, c(1:(0+m1.m0$col[i]), 20, 21:(20+m1.m0$row[i]))], type="response")
  logistic.pred2 <- as.factor(ifelse(pi.hat > 0.5, "non_G1", "G1"))
  test.error[i] <-  mean(logistic.pred2 != gene.epsilon$test[,20])
  print(table(logistic.pred2, gene.epsilon$test$Peak))
}
test.error # 0.09836066
cv.error.sd[m1.m0$row, m1.m0$col] # 0.03337084

best.m <- 7
logistic.fit2 <- glm(Peak ~ ., data=gene.epsilon$train[, c(1:(0+best.m), 20, 21:(20+best.m))], family="binomial")
logistic.pred2 <- as.factor(ifelse(predict(logistic.fit2, gene.epsilon$test[, c(1:(0+best.m), 20, 21:(20+best.m))], type="response") > 0.5,
                                   "non_G1", "G1"))
mean(logistic.pred2 != gene.epsilon$test[,20])
table(logistic.pred2, gene.epsilon$test[,20])

# with first & second derivative
library(doParallel)
nc <- detectCores()
registerDoParallel(nc-1)
cv.error.for.m <- array(0, c(m,m,m))
result <- foreach (j = 1:K) %dopar% {
  cv.error.for.m <- array(0, c(m,m,m))
  for (i in 1:m) {
    for( k in 1:m) {
      # pred.error <- numeric(K)
      for (h in 1:m) {
        model <- glm(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20, 21:(20+k), 39:(38+h))], family="binomial")
        pi.hat <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20, 21:(20+k), 39:(38+h))], type="response")
        pred <- as.factor(ifelse(pi.hat > 0.5, "non_G1", "G1"))
        # pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
        cv.error.for.m[h,k,i] <- mean(pred != gene.epsilon$train[folds == j, 20])# mean(pred.error)
      } 
    }
  }
  return(cv.error.for.m)
}
str(result)
# saveRDS(result, file="deriv2_cv_error.RData")
result <- readRDS("deriv2_cv_error.RData")
cv.error.for.m <- array(0, c(m,m,m))
for(i in 1:K) {
  cv.error.for.m <- cv.error.for.m + result[[i]]
}
cv.error.for.m <- cv.error.for.m / 10
m2.m1.m0 <- as.data.frame(which(cv.error.for.m == min(cv.error.for.m), arr.ind=T))
m2.m1.m0
min(cv.error.for.m) # 0.0817602

test.error <- numeric(nrow(m2.m1.m0))
for (i in 1:nrow(m2.m1.m0)) {
  logistic.fit2 <- glm(Peak ~ .,
                       data=gene.epsilon$train[,  c(1:( 0+m2.m1.m0$dim3[i]), 20,
                                                   21:(20+m2.m1.m0$dim2[i]),
                                                   39:(38+m2.m1.m0$dim1[i]))], family="binomial")
  pi.hat <- predict(logistic.fit2, gene.epsilon$test[, c(1:( 0+m2.m1.m0$dim3[i]), 20,
                                                        21:(20+m2.m1.m0$dim2[i]),
                                                        39:(38+m2.m1.m0$dim1[i]))], type="response")
  logistic.pred2 <- as.factor(ifelse(pi.hat > 0.5, "non_G1", "G1"))
  test.error[i] <-  mean(logistic.pred2 != gene.epsilon$test[,20])
  print(table(logistic.pred2, gene.epsilon$test$Peak))
}
test.error # 0.10655738 0.10655738 0.10655738 0.10655738 0.09836066 0.09836066 0.09836066 0.09836066 0.10655738

best.m <- 7
logistic.fit2 <- glm(Peak ~ ., data=gene.epsilon$train[, c(1:(0+best.m), 20, 21:(20+best.m))], family="binomial")
logistic.pred2 <- as.factor(ifelse(predict(logistic.fit2, gene.epsilon$test[, c(1:(0+best.m), 20, 21:(20+best.m))], type="response") > 0.5,
                                   "non_G1", "G1"))
mean(logistic.pred2 != gene.epsilon$test[,20])
table(logistic.pred2, gene.epsilon$test[,20])


# svm
library(e1071)
set.seed(100)
radial.svm <-tune(svm, Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], kernel="radial",
                  ranges=list(cost=c(0.001, 0.01, 0.1, 1, 5, 10, 20, 50), gamma=c(0.001, 0.01, 0.1, 1, 10, 50)) )
radial.svm # 0.09812925

set.seed(100)
radial.svm <-tune(svm, Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], kernel="radial",
                  ranges=list(cost=c(30, 40, 50, 60, 70), gamma=c(1e-5, 1e-4, 0.001, 0.01)) )
radial.svm # 0.08796769

set.seed(100)
radial.svm <-tune(svm, Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], kernel="radial",
                  ranges=list(cost=c(55, 57, 60, 62, 65), gamma=c(0.00005, 0.0001, 0.0005)) )
radial.svm # 0.08792517
svm.pred <- predict(radial.svm$best.model, gene.epsilon$test[, c(1:5, 20)])
mean(best.pred !=  gene.epsilon$test[,20])
table(best.pred, gene.epsilon$test[,20])

linear.svm <-tune(svm, Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], kernel="linear",
                  ranges=list(cost=c(0.001, 0.01, 0.1, 1, 5, 10, 20, 50)))
linear.svm # 0.09005102
best.pred <- predict(linear.svm$best.model, gene.epsilon$test[, c(1:5, 20)])
mean(best.pred !=  gene.epsilon$test[,20]) # 0.08196721
table(best.pred, gene.epsilon$test[,20])

polynomial.svm <-tune(svm, Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], kernel="polynomial",
                      ranges=list(cost=c(0.001, 0.01, 0.1, 1,10, 15, 20, 40), d=c(2, 3, 4, 5)) )
polynomial.svm # 0.1085034

sigmoid.svm <-tune(svm, Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], kernel="sigmoid",
                      ranges=list(cost=c(0.001, 0.01, 0.1, 1,10, 15, 20, 40), gamma=c(0.001, 0.01, 0.1, 1, 10, 50)) )
sigmoid.svm # 0.09808673

# qda
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- qda(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20)])
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20)], type="response")$class
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m

# lda
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- lda(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20)])
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20)], type="response")$class
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m

# random forest
library(ranger)
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    rf.fit <- ranger(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20)], mtry=i, num.trees=200, max.depth=3)
    pred <- predict(rf.fit, data=gene.epsilon$train[folds == j, c(1:i, 20)], type='response')$predictions
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m # 0.1718963 0.1575255 0.1494473 0.1494048 0.1514031 0.1431973 0.1473214 0.1514031 0.1596088

# misclassification analysis
logistic.false <- which(logistic.pred1 != gene.epsilon$test$Peak)
true.peak <- gene.epsilon$test$Peak[logistic.false]
false.G1 <- logistic.false[which(logistic.pred1[logistic.false] == "G1")]
false.non_G1 <- logistic.false[which(logistic.pred1[logistic.false] == "non_G1")]
test.false.G1 <- test.fdata
test.false.G1$data <- test.fdata$data[false.G1,]
test.false.non_G1 <- test.fdata
test.false.non_G1$data <- test.fdata$data[false.non_G1,]

# green line:TRUE G1, red line:TRUE non_G1
test.G1 <- test.fdata
test.G1$data <- test.fdata$data[test.label$Peak == "G1",]
test.non_G1 <- test.fdata
test.non_G1$data <- test.fdata$data[test.label$Peak == "non_G1",]

plot.fdata(test.G1, col=1, xlab="time points", ylab="gene expression",
           ylim=c(-3, 3), main="TRUE G1(black) and Misclassified G1(red)")
for (i in 1:length(false.G1)) {
  lines(1:18, test.false.G1$data[i,], type="l", col=2)
}
plot.fdata(test.non_G1, col=1, xlab="time points", ylab="gene expression",
           ylim=c(-3, 3), main="TRUE non-G1(black) and Misclassified G1(red)")
for (i in 1:length(false.G1)) {
  lines(1:18, test.false.G1$data[i,], type="l", col=2)
}


plot.fdata(test.G1, col=1, xlab="time points", ylab="gene expression",
           ylim=c(-3, 3), main="TRUE G1(black) and Misclassified non-G1(red)")
for (i in 1:length(false.non_G1)) {
  lines(1:18, test.false.non_G1$data[i,], type="l", col=2)
}
plot.fdata(test.non_G1, col=1, xlab="time points", ylab="gene expression",
           ylim=c(-3, 3), main="TRUE non-G1(black) and Misclassified non-G1(red)")
for (i in 1:length(false.non_G1)) {
  lines(1:18, test.false.non_G1$data[i,], type="l", col=2)
}

svm.false <- which(svm.pred != gene.epsilon$test$Peak)
true.peak <- gene.epsilon$test$Peak[svm.false]
false.G1 <- svm.false[which(svm.pred[svm.false] == "G1")]
false.non_G1 <- svm.false[which(svm.pred[svm.false] == "non_G1")]
test.false.G1 <- test.fdata
test.false.G1$data <- test.fdata$data[false.G1,]
test.false.non_G1 <- test.fdata
test.false.non_G1$data <- test.fdata$data[false.non_G1,]


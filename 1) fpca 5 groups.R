gene <- read.csv("yeast_with_label.csv")
str(gene)
head(gene)
gene <- na.omit(gene[, c(2, 7:24, 80)])
n <- nrow(gene)
set.seed(100)
train <- sample(n, round(n*0.8))

library(fda.usc)
str(gene)
train.label <- gene[train, c(1,20)]
test.label <- gene[-train, c(1,20)]
train.fdata <- fdata(gene[train, 2:19])
test.fdata <- fdata(gene[-train, 2:19])
train.deriv <- fdata.deriv(train.fdata, nderiv=1)
test.deriv <- fdata.deriv(test.fdata, nderiv=1)
train.deriv2 <- fdata.deriv(train.fdata, nderiv=2)
test.deriv2 <- fdata.deriv(test.fdata, nderiv=2)

train.deriv1 <- train.deriv
train.deriv1$data <- train.deriv$data[1:5,]
plot.fdata(train.fdata, col=1)
plot.fdata(train.deriv1)

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

curve2epsilon <- function(train, test, ncomp=5) {
  library(fda.usc)
  train.pca <- create.pc.basis(train, l=1:ncomp)
  
  train.e <- get.epsilon(train, train.pca)
  test.e <- get.epsilon(test, train.pca)
  
  result <- list(train=train.e, test=test.e)
  return(result)
}

m <- 18
gene.epsilon <- curve2epsilon(train.fdata, test.fdata, ncomp=m)
fpc.train <- fdata2pc(train.fdata, ncomp=18, lambda=0)
gene.epsilon$train <- fpc.train$x
gene.deriv.epsilon <- curve2epsilon(train.deriv, test.deriv, ncomp=m)
gene.deriv2.epsilon <- curve2epsilon(train.deriv2, test.deriv2, ncomp=m)

gene.epsilon$train <- cbind(gene.epsilon$train, train.label)
gene.epsilon$train <- cbind(gene.epsilon$train, gene.deriv.epsilon$train)
gene.epsilon$train <- cbind(gene.epsilon$train, gene.deriv2.epsilon$train)
gene.epsilon$test <- cbind(gene.epsilon$test, test.label)
gene.epsilon$test <- cbind(gene.epsilon$test, gene.deriv.epsilon$test)
gene.epsilon$test <- cbind(gene.epsilon$test, gene.deriv2.epsilon$test)


library(nnet)
K <- 10
index <- (1:nrow(gene.epsilon$train))%%K + 1
set.seed(100)
folds <- sample(index, nrow(gene.epsilon$train))
cv.error.for.m <- numeric(m)
# curve only
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- multinom(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20)], maxit=1000)
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20)], type="class")
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m
multi.fit <- multinom(Peak ~ ., data=gene.epsilon$train[, c(1:5, 20)], maxit=1000)
multi.pred <- predict(multi.fit, gene.epsilon$test[, c(1:5, 20)], type="class")
mean(multi.pred != gene.epsilon$test[,20])
table(multi.pred, gene.epsilon$test[,20])

incorrect.index <- which(multi.pred != gene.epsilon$test[,20])
test.incorrect <- test.fdata
test.incorrect$data <- test.fdata$data[incorrect.index,]
train.pca <- create.pc.basis(train.fdata)
test.label[incorrect.index, ]
test.label[-incorrect.index, ][test.label[-incorrect.index, 2] == "M/G1",]
plot.train.data <- as.data.frame(t(train.fdata$data))
plot.test.data <- as.data.frame(t(test.fdata$data))

sub.data <- train.fdata
sub.data$data <- train.fdata$data[1:5,]
plot.fdata(sub.data, col=1:5, xlab="time points", ylab="gene expression", ylim=c(-2, 2), main=NULL)

library(ggplot2)
colnames(plot.train.data)
ggplot(data=plot.train.data, aes(x=seq(0,119,len=18), y=plot.train.data[,2])) + geom_line() + 
  ylim(-2, 2) + xlab("time (min)") + ylab("gene expression")
  geom_line(aes(y=`106`)) +
  geom_line(aes(y=`127`)) +
  geom_line(aes(y=`386`)) +
# geom_line(data=as.data.frame(t(test.incorrect$data)), aes(x=1:18, y=`30`, color='test 30'))
# ggplot(data=plot.test.data, aes(x=1:18, y=`609`)) + geom_line() + ylim(-4, 4) +
  geom_line(aes(y=`650`)) +
  geom_line(aes(y=`708`)) +
  geom_line(aes(y=`709`)) +
  geom_line(data=as.data.frame(t(test.incorrect$data)), aes(x=1:18, y=`195`, color='test 195'))

# svm
library(e1071)
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- svm(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20)], kernel="sigmoid")
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20)], type="response")
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  }
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m
# linear, radial은 최소 에러가 0.2986395로 multinom보다 조금 높음
# polynomial은 0.4315901, sigmoid는 0.3354592로 안좋음 m=5일 때 에러 최소
# radial kernel tuning
radial.tune <- tune(svm, Peak ~ ., data=gene.epsilon$train[folds != j, c(1:5, 20)], kernel="radial",
                    ranges=list(cost=c(4, 5, 6, 7, 8, 9), gamma=c(0.01, 0.02, 0.03, 0.04, 0.05)) )
# - best parameters:
#   cost gamma
#      6  0.05
# 
# - best performance: 0.2636364 
radial.tune <- tune(svm, Peak ~ ., data=gene.epsilon$train[folds != j, c(1:5, 20)], kernel="radial",
                    ranges=list(cost=c(4, 5, 6, 7, 8, 9), gamma=c(0.05, 0.1, 0.5, 1, 10)) )

best.pred <- predict(radial.tune$best.model, gene.epsilon$test[, c(1:5, 20)])
mean(best.pred !=  gene.epsilon$test[,20])
table(best.pred, gene.epsilon$test[,20])


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
cv.error.for.m # multinom보다 안좋음

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
cv.error.for.m

## classification with derivatives (no derevative data랑 별 차이 없음)
# multinomial with first derivative
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- multinom(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20, 21:(20+i))])
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20, 21:(20+i))], type="class")
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m
multi.fit <- multinom(Peak ~ ., data=gene.epsilon$train[, c(1:5, 20, 21:25)])
multi.pred <- predict(multi.fit, gene.epsilon$test[, c(1:5, 20, 21:25)], type="class")
mean(multi.pred != gene.epsilon$test[,20])
table(multi.pred, gene.epsilon$test[,20])

# multinomial with first, second derivative
colnames(gene.epsilon$train)[39] # e1 -> FPC1 of 2nd derivative curves
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- multinom(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20, 21, 39)])
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20, 21, 39)], type="class")
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m
multi.fit <- multinom(Peak ~ ., data=gene.epsilon$train[, c(1:5, 20, 21, 39)])
multi.pred <- predict(multi.fit, gene.epsilon$test[, c(1:5, 20, 21, 39)], type="class")
mean(multi.pred != gene.epsilon$test[,20])
table(multi.pred, gene.epsilon$test[,20])

# svm with first derivative
cv.error.for.m <- numeric(m)
for (i in 1:m) {
  pred.error <- numeric(K)
  for (j in 1:K) {
    model <- svm(Peak ~ ., data=gene.epsilon$train[folds != j, c(1:i, 20, 21:(20+i))], kernel="radial", cost=6, gamma=0.05)
    pred <- predict(model, gene.epsilon$train[folds == j, c(1:i, 20, 21:(20+i))], type="response")
    pred.error[j] <- mean(pred != gene.epsilon$train[folds == j, 20])
  } 
  cv.error.for.m[i] <- mean(pred.error)
}
cv.error.for.m

radial.tune <- tune(svm, Peak ~ ., data=gene.epsilon$train[folds != j, c(1:5, 20, 21:25, 39:44)], kernel="polynomial",
                    ranges=list(cost=c(10, 15, 20, 40), d=c(2, 3, 4, 5)) )
# - best parameters: for radial
#   cost gamma
#      7  0.03
# - best performance: 0.2659091 

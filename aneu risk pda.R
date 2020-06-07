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

get.pw.const.beta <- function(fdata., m) {
  library(fda.usc)
  n <- nrow(fdata.$data)
  time.points <- ncol(fdata.$data)
  n.order <- ifelse(m>=4, m+1, 4)
  
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

get.RSQ <- function(fdata., coefs) {
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

library(fda.usc)
FKS.fd <- fdata(FKS)
FKS.fd$argvals <- seq(-32.9, -0.74, len=341)
FKS.fd$rangeval <- c(-32.9, -0.74)
plot.fdata(FKS.fd, col=1)

n <- nrow(FKS.fd$data) # 65
test.errors <- matrix(0, 50, 4)
for (order in 1:4) {
  set.seed(100)
    for (i in 1:50) {
    train <- sample(1:n, 55)
    FKS.train <- FKS.fd
    FKS.train.U <- FKS.fd
    FKS.train.L <- FKS.fd
    FKS.train$data <- FKS.fd$data[train,]
    FKS.train.U$data <- FKS.fd$data[train[which(patient$type[train] == "U")],]
    FKS.train.L$data <- FKS.fd$data[train[which(patient$type[train] == "L")],]
    FKS.test <- FKS.fd
    FKS.test$data <- FKS.fd$data[-train,]
    U.beta <- get.pw.const.beta(FKS.train.U, m=order)
    L.beta <- get.pw.const.beta(FKS.train.L, m=order)
    U.RSQ <- get.RSQ(FKS.test, U.beta)
    L.RSQ <- get.RSQ(FKS.test, L.beta)
    RSQ.df <- data.frame(U=U.RSQ, L=L.RSQ)
    pred.g <- as.factor(append(ifelse(apply(RSQ.df, 1, which.max) == 1, "L", "U"), c("L", "U")))[1:10]
    test.errors[i, order] <- mean(pred.g != patient$type[-train])
    }
  print(paste("order ", " complete", sep=as.character(order)))
}
test.errors
apply(test.errors, 2, mean)
apply(test.errors, 2, sd)

FKS.test <- FKS.fd
FKS.test$data <- FKS.fd$data[-train,]
FKS.test.d1 <- fdata.deriv(FKS.test, nderiv=1)
FKS.test.d2 <- fdata.deriv(FKS.test, nderiv=2)
FKS.test.d3 <- fdata.deriv(FKS.test, nderiv=3)

FKS.train.U.2nd <- get.pw.const.beta(FKS.train.U, m=2)
FKS.train.L.2nd <- get.pw.const.beta(FKS.train.L, m=2)
test.RSQ.U <- (apply((FKS.test.d2$data)^2, 1, sum) - 
                 apply(((FKS.test$data * FKS.train.U.2nd[1]) + 
                          ((FKS.test.d1$data) * FKS.train.U.2nd[2]) + 
                          FKS.test.d2$data)^2, 1, sum)) / apply((FKS.test.d2$data)^2, 1, sum)

test.RSQ.L <- (apply((FKS.test.d2$data)^2, 1, sum) - 
                 apply(((FKS.test$data * FKS.train.L.2nd[1]) + 
                          ((FKS.test.d1$data) * FKS.train.L.2nd[2]) + 
                          FKS.test.d2$data)^2, 1, sum)) / apply((FKS.test.d2$data)^2, 1, sum)
test.RSQ <- data.frame(U=test.RSQ.U, L=test.RSQ.L)
pred.g <- as.factor(ifelse(apply(test.RSQ, 1, which.max) == 1, "U", "L"))
table(pred.g, patient$type[-train])
mean(pred.g != patient$type[-train]) # 0.2

FKS.train.U.3rd <- get.pw.const.beta(FKS.train.U, m=3)
FKS.train.L.3rd <- get.pw.const.beta(FKS.train.L, m=3)
test.RSQ.U <- (apply((FKS.test.d3$data)^2, 1, sum) - 
               apply(((FKS.test$data * FKS.train.U.3rd[1]) + 
                     ((FKS.test.d1$data) * FKS.train.U.3rd[2]) +
                     ((FKS.test.d2$data) * FKS.train.U.3rd[3]) +
                       FKS.test.d3$data)^2, 1, sum)) / apply((FKS.test.d3$data)^2, 1, sum)

test.RSQ.L <- (apply((FKS.test.d3$data)^2, 1, sum) - 
                 apply(((FKS.test$data * FKS.train.L.3rd[1]) + 
                          ((FKS.test.d1$data) * FKS.train.L.3rd[2]) +
                          ((FKS.test.d2$data) * FKS.train.L.3rd[3]) +
                          FKS.test.d3$data)^2, 1, sum)) / apply((FKS.test.d3$data)^2, 1, sum)
test.RSQ <- data.frame(U=test.RSQ.U, L=test.RSQ.L)
pred.g <- as.factor(ifelse(apply(test.RSQ, 1, which.max) == 1, "U", "L")) # SSE는 U로만 예측함
table(pred.g, patient$type[-train])
mean(pred.g != patient$type[-train]) # 0.6


# pda score classification
FKS.train.2nd <- get.pw.const.beta(FKS.train, m=2)
FKS.train.3rd <- get.pw.const.beta(FKS.train, m=3)
sol.2nd <- polyroot(c(-FKS.train.2nd, 1))
sol.3rd <- polyroot(c(-FKS.train.3rd, 1))

mean.curve <- apply(FKS.train$data, 2, mean)
diff <- as.matrix(FKS.fd$data) - matrix(mean.curve, nrow=65, ncol=341, byrow=T)
times <- FKS.train$argvals
psi.1.2nd <- exp(Re(sol.2nd[1])*times)*cos(Im(sol.2nd[1])*times)
psi.2.2nd <- exp(Re(sol.2nd[1])*times)*sin(Im(sol.2nd[1])*times)
score1.2nd <- diff %*% matrix(psi.1.2nd, ncol=1)
score2.2nd <- diff %*% matrix(psi.2.2nd, ncol=1)

psi.1.3rd <- exp(Re(sol.3rd[1])*times)
psi.2.3rd <- exp(Re(sol.3rd[2])*times)*cos(Im(sol.3rd[2])*times)
psi.3.3rd <- exp(Re(sol.3rd[2])*times)*sin(Im(sol.3rd[2])*times)
score1.3rd <- diff %*% matrix(psi.1.3rd, ncol=1)
score2.3rd <- diff %*% matrix(psi.2.3rd, ncol=1)
score3.3rd <- diff %*% matrix(psi.3.3rd, ncol=1)

scores.2nd <- data.frame(y=patient$type, score1=score1.2nd, score2=score2.2nd)
scores.3rd <- data.frame(y=patient$type, score1=score1.3rd, score2=score2.3rd, score3=score3.3rd)

logistic.2nd <- glm(y~., data=scores.2nd[train,], family=binomial)
logistic.3rd <- glm(y~., data=scores.3rd[train,], family=binomial)
pi.hat <- predict(logistic.2nd, newdata=scores.2nd[-train,], type="response")
pred.2nd <- as.factor(ifelse(pi.hat > 0.5, "U", "L"))
pi.hat <- predict(logistic.3rd, newdata=scores.3rd[-train,], type="response")
pred.3rd <- as.factor(ifelse(pi.hat > 0.5, "U", "L"))
table(pred.2nd, scores.2nd$y[-train]) # 2차 미방 error rate= 0.2
table(pred.3rd, scores.3rd$y[-train]) # 3차 미방 error rate= 0.1

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

test.errors <- matrix(0, 50, 4)
for (order in 1:4) {
  set.seed(100)
  for (i in 1:50) {
    train <- sample(1:n, 55)
    FKS.train <- FKS.fd
    FKS.train$data <- FKS.fd$data[train,]
    FKS.test <- FKS.fd
    FKS.test$data <- FKS.fd$data[-train,]
    coefs <- get.pw.const.beta(FKS.train, m=order)
    psi.list <- get.psi(coefs, FKS.train$argvals)
    mean.curve <- apply(FKS.train$data, 2, mean)
    diff <- as.matrix(FKS.fd$data) - matrix(mean.curve, nrow=n, ncol=length(FKS.fd$argvals), byrow=T)
    scores <- list()
    scores[["y"]] <- patient$type
    for (j in 1:length(psi.list)) {
      scores[[as.character(j)]] <- diff %*% matrix(psi.list[[j]], ncol=1)
    }
    scores <- as.data.frame(scores)
    
    qda.fit <- qda(y~., data=scores[train,])
    pi.hat <- predict(qda.fit, newdata=scores[-train,])$posterior[,1]
    pred <- as.factor(c(ifelse(pi.hat < 0.5, "U", "L"), "U", "L"))[1:10]
    # logistic.fit <- glm(y~., data=scores[train,], family=binomial)
    # pi.hat <- predict(logistic.fit, newdata=scores[-train,], type="response")
    # pred <- as.factor(ifelse(pi.hat > 0.5, "U", "L"))
    test.errors[i, order] <- mean(pred != scores$y[-train])
  }
  print(paste("order ", " complete", sep=as.character(order)))
}
test.errors
apply(test.errors, 2, mean)
apply(test.errors, 2, sd)
qda.test.errors <- apply(test.errors, 2, mean)
logistic.test.errors
# FPCA classification
get.fpc.score <- function(fdata, fdata.comp) {
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

fpc.test.errors <- matrix(0, 50, 5)
for (fpcs in 1:5) {
  set.seed(100)
  for (i in 1:50) {
    train <- sample(1:n, 55)
    FKS.train <- FKS.fd
    FKS.train$data <- FKS.fd$data[train,]
    FKS.train.fpca <- create.pc.basis(FKS.train, l=1:fpcs, lambda=1)
    scores <- get.fpc.score(FKS.fd, FKS.train.fpca)
    scores$y <- patient$type
    
    logistic.fit <- glm(y~., data=scores[train,], family=binomial)
    pi.hat <- predict(logistic.fit, newdata=scores[-train,], type="response")
    pred <- as.factor(ifelse(pi.hat > 0.5, "U", "L"))
    fpc.test.errors[i, fpcs] <- mean(pred != scores$y[-train])
  }
  print(paste("fpc ", " complete", sep=as.character(fpcs)))
}
fpc.test.errors
fpc.test.errors. <- apply(fpc.test.errors, 2, mean)
apply(fpc.test.errors, 2, sd)

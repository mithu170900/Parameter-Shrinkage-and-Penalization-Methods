#-----------------------------------------------------------------------------------------------
# 1. Spacing delta: Make the first two values of mu zero, then the last one depends on delta
# 2. Penalty function: Bring mu values closer together
#-----------------------------------------------------------------------------------------------

# Example of cross-validation
# Two datasets from normal distributions with common variance

# Normal distribution log-likelihood function
norm_llh <- function(data,par) {
  x <- data
  mu <- par
  n <- length(x)
  llh <- -n/2*log(2*pi) - sum((x-mu)^2)
  return(llh)
}

# Penalized distance function --> distance to force mu-values to be close
pen_dist <- function(data,par) {
  X1 <- data$X1
  X2 <- data$X2
  X3 <- data$X3
  n <- length(X1)
  lambda <- data$lambda
  mu1 <- par[1]
  mu2 <- par[2]
  mu3 <- par[3]
  #D <- -norm_llh(X1,mu1) -norm_llh(X2,mu2) -norm_llh(X3,mu3) + lambda*n*((mu1-mu2)^2+(mu1-mu3)^2+(mu2-mu3)^2)
  D <- -norm_llh(X1,mu1) -norm_llh(X2,mu2) -norm_llh(X3,mu3) + lambda*n*((mu1)^2+(mu2)^2+(mu3)^2)
  return(D)
}

# Creates random permutation for test and validation data
perm_score <- function(X1,X2,X3,lambda.grid) {
  n <- length(X1)
  m <- floor(0.75*n)
  X1.perm <- sample(X1)
  X2.perm <- sample(X2)
  X3.perm <- sample(X3)
  X1.test <- X1.perm[1:m]
  X1.valid <- X1.perm[(m+1):n]
  X2.test <- X2.perm[1:m]
  X2.valid <- X2.perm[(m+1):n]
  X3.test <- X3.perm[1:m]
  X3.valid <- X3.perm[(m+1):n]
  n.grid <- length(lambda.grid)
  V <- rep(0,n.grid)
  for (j in 1:n.grid) {
    data.in <- list(X1 = X1.test, X2 = X2.test, X3 = X3.test, lambda = lambda.grid[j])
    D <- optim(c(mean(X1.test),mean(X2.test),mean(X3.test)), fn = pen_dist, data = data.in, method = "BFGS")
    V[j] <- -norm_llh(X1.valid,D$par[1]) -norm_llh(X2.valid,D$par[2]) -norm_llh(X3.valid,D$par[3])
  }
  return(V)
}

#lambda.grid <- seq(0,10,length.out = 100)
#V <- perm_score(X1,X2,lambda.grid)
#plot(lambda.grid,V,type = "l")

# Simulate data, then use above code to generate 
# M random partitions and for each store the distance
# function evaluated on the random partition validation data
X1 <- rnorm(n = 50, mean = 0, sd = 1)
X2 <- rnorm(n = 50, mean = 0, sd = 1)
X3 <- rnorm(n = 50, mean = 1, sd = 1)
mu.mle <- c(mean(X1),mean(X2),mean(X3))
M <- 25
nL <- 100
lambda.grid <- seq(0,0.5,length.out = nL)
V.store = matrix(rep(0,M*nL),nrow = M)
for (m in 1:M) {
  V <- perm_score(X1,X2,X3,lambda.grid)
  V.store[m,] <- V
}
# Now average the "validation" scores and find the minimum
V.mean <- colMeans(V.store)
plot(lambda.grid,V.mean,type = "l")
lambda.out <- lambda.grid[V.mean==min(V.mean)]
# Finally, evaluate estimators using full data at the optimal lambda
data.opt <- list(X1 = X1, X2 = X2, X3 = X3, lambda = lambda.out)
mu.out <- optim(c(mean(X1),mean(X2),mean(X3)),fn = pen_dist, data = data.opt)$par
mu.mle
mu.out

# Adding a simulation study component:
# (1) Simulate a dataset
# (2) Find both MLE and penalized cross-validation est's
# (3) Store the estimators for each simulated dataset
# (4) Repeat this a large # of times so we can calculate
#     MSE (or RMSE) for the two methods

nDelta <- 5
delta.grid <- seq(0,3,length.out = nDelta)
MSEpen <- rep(0,nDelta)
MSEmle <- rep(0,nDelta)
for (j in 1:nDelta) {
  K <- 5 # Number of datasets to simulate
  mle.est <- matrix(rep(0,3*K),nrow = K)
  pen.est <- matrix(rep(0,3*K),nrow = K)
  mu.true <- c(0, 0, delta.grid[j])
  for (k in 1:K) {
    X1 <- rnorm(n = 50, mean = mu.true[1], sd = 1)
    X2 <- rnorm(n = 50, mean = mu.true[2], sd = 1)
    X3 <- rnorm(n = 50, mean = mu.true[3], sd = 1)
    mu.mle <- c(mean(X1),mean(X2),mean(X3))
    M <- 50 # Number of partitions
    nL <- 100 # Number of lambda values
    lambda.grid <- seq(0,0.5,length.out = nL)
    V.store = matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(X1,X2,X3,lambda.grid)
      V.store[m,] <- V
    }
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    data.opt <- list(X1 = X1, X2 = X2, X3 = X3, lambda = lambda.out)
    mu.out <- optim(c(mean(X1),mean(X2), mean(X3)),fn = pen_dist, data = data.opt)$par
    mle.est[k,] <- mu.mle
    pen.est[k,] <- mu.out
  }
mle.bias <- colMeans(mle.est) - mu.true
pen.bias <- colMeans(pen.est) - mu.true
mle.var <- apply(mle.est,2,var) # calculates var for each column
pen.var <- apply(pen.est,2,var) # calculates var for each column
  #mle.bias
  #pen.bias
  #mle.var
  #pen.var
mle.MSE <- mle.bias^2 + mle.var
pen.MSE <- pen.bias^2 + pen.var
  #mle.MSE
  #pen.MSE
MSEmle[j] <- sum(mle.MSE)
MSEpen[j] <- sum(pen.MSE)
print(j)
}
plot(delta.grid,MSEpen,type = "l")
lines(delta.grid,MSEmle,col = "red")
nllh <- function(x, mu) {
  n <- length(x)
  negllh <- - n*log(1/sqrt(2*pi)) + sum((x-mu)^2)
  return(negllh)
}

pen_dist <- function(data,lambda,pen_type,par) {
  # Adapt this function so that if pen_type = "ToZero" then it shrinks all the mu values closer to zero and if 
  # pen_type = "Closer" then it shrinks all the mu values closer to one another. Use an "if" function
  
  X <- data$X
  index<- data$index
  n <- sum(index==1) # Assuming equal group sizes
  mu <- par
  J <- max(index)
  D.vals <- rep(0,J)
  for (j in 1:J) {
    D.vals[j] <- nllh(X[index==j],mu[j])
  }

  if (pen_type == "ToZero") {
    D <- sum(D.vals) + lambda*n*sum(mu^2) # Penalty on squared mu-values
  } else if (pen_type == "Closer") {
    mu.matrix <- matrix(rep(mu,J),nrow = J,byrow = TRUE)
    D <- sum(D.vals) + lambda*n*sum((mu.matrix-t(mu.matrix))^2)/2
  }
  return(D)
}
#-----------------------------------------------------------------------------------------------
#mu.mle <- rep(0,J)
#for (j in 1:J) {
#  mu.mle[j] = mean(data$X[data$index==j])
#}
#
#optim(mu.mle, fn = pen_dist, data = data, lambda = 10)
#-----------------------------------------------------------------------------------------------
# Creates random permutation for test and validation data

perm_score <- function(dataset, lambda.grid, pen_type) {
  index <- dataset$index
  X <- dataset$X
  J <- max(index)
  n <- length(X[index==1]) # Number of values
  m <- floor(0.75*n)
  # Iterate through each index
  # For each index, we have 50 values and we will split accordingly 
  # Don't iterate through the values

dataset.test <- data.frame(index= c(), X = c())
dataset.valid <- data.frame(index= c(), X = c())

  for (a in 1:J) {
    X_ran <- sample(X[index == a])
    X_test <- data.frame(index= rep(a, m), X = X_ran[1:m])
    dataset.test <- rbind(dataset.test, X_test)
    X_valid <- data.frame(index= rep(a, n-m), X = X_ran[(m+1):n])
    dataset.valid <- rbind(dataset.valid, X_valid)
  }
  nL <- length(lambda.grid)
  mu.in <- c()
  for (a in 1:J) {
    mu.in <- c(mu.in,mean(dataset.test$X[dataset.test$index== a])) #Calculate mean of each subset
  }
  V <- rep(0,nL)
  for (b in 1:nL) {
    # data.in <- list(X = , lambda = lambda.grid[b])
    D <- optim(par = mu.in, fn = pen_dist, data = dataset.test,lambda = lambda.grid[b],pen_type = pen_type, method = "BFGS")
    # Validation score calculated on the validation data for each lambda value
    V[b] <- pen_dist(data = dataset.valid,lambda = 0, pen_type = pen_type,par = D$par )
  }
  return(V)
}

# Adding a simulation study component:
# (1) Simulate a dataset
# (2) Findexboth MLE and penalized cross-validation est's
# (3) Store the estimators for each simulated dataset
# (4) Repeat this a large # of times so we can calculate
#     MSE (or RMSE) for the two methods
J <- 10 # Number of variables
delta <- 3
mu.vals <- seq(0,delta,length.out = J) # Equally spaced
N <- 50 # Number of values in each variable

# Generate data -- insert code here (specify sample size)
index <- rep(1,N)
X <- rnorm(N,mu.vals[1],1)
for (j in 2:10) {
  index <- c(index,rep(j,N))
  X <- c(X,rnorm(n, mean = mu.vals[j], sd = 1))
}
dataset <- data.frame(index = index, X = X)

lambdaMax = 0.5
nL = 100
M <- 15 # Number of partitions
lambda.grid <- seq(0,lambdaMax,length.out = nL)

V.store <- matrix(rep(0,M*nL),nrow = M)
for (m in 1:M) {
  V <- perm_score(dataset,lambda.grid, "Closer")
  V.store[m,] <- V
}
V.mean <- colMeans(V.store)
plot(lambda.grid, V.mean, type = "l")

nDelta <- 5 # Number of delta values
#delta.grid <- seq(0,3,length.out = nDelta)

MSEpen <- rep(0,nDelta)
MSEmle <- rep(0,nDelta)

for (j in 1:nDelta) {
  start_time <- Sys.time()
  K <- 2 # Number of datasets to simulate
  delta.grid <- seq(0,1,length.out = nDelta)
  mle.est <- matrix(rep(0,J*K),nrow = K)
  pen.est <- matrix(rep(0,J*K),nrow = K)
  mu.true <- seq(0,max(delta.grid), length.out = J)
  mu.mle <- c()
  for (a in 1:J) {
    mu.mle <- c(mu.mle, mean(dataset.test$X[dataset.test$index== a])) #Calculate mean of each subset
  }
  for (k in 1:K) {
    lambda.grid <- seq(0,lambdaMax,length.out = nL)
    V.store <- matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(dataset,lambda.grid, "Closer")
      V.store[m,] <- V
    }
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    data.opt <- list(dataset, lambda = lambda.out)
    mu.out <- optim(par = mu.mle,fn = pen_dist, data = data.test, lambda = lambda.grid[k], pen_type = pen_type, method = "BFGS")$par
    mle.est[k,] <- mu.mle  ######
    pen.est[k,] <- mu.out
  }
  mle.bias <- colMeans(mle.est) - mu.true
  pen.bias <- colMeans(pen.est) - mu.true
  mle.var <- apply(mle.est,2,var) # calculates var for each column ########
  pen.var <- apply(pen.est,2,var) # calculates var for each column
  mle.MSE <- mle.bias^2 + mle.var 
  pen.MSE <- pen.bias^2 + pen.var
  MSEmle[j] <- sum(mle.MSE)######
  MSEpen[j] <- sum(pen.MSE)######
  end_time <- Sys.time()
  timer <- end_time - start_time
  print(j)
  print(timer)
}
plot(delta.grid,MSEpen,type = "l")
lines(delta.grid,MSEmle,col = "red")

results <- list(delta = delta.grid,
                MSEpen = MSEpen,
                MSEmle = MSEmle)

#name <- "delta_eq_spaced_closer_n50"
name <- "delta_one_nonzero_closer_n50"
save(results, file = paste(name, ".RData", sep=""))

# ASSIGNMENT: 
#   - Work on finding the min of each pen_type by changing delta and lambda values
#   - Automate code for the simulation part


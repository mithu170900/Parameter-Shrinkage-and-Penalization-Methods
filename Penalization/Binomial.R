binomial_nllh <- function(data, par, N) {
  y <- data
  p <- par
  nllh <- - sum(dbinom(y, N, p, log = TRUE))
  return(nllh)
}

pen_dist <- function(data,lambda,pen_type,par,N) {
  # Adapt this function so that if pen_type = "ToZero" then it shrinks all the p values closer to zero and if 
  # pen_type = "Closer" then it shrinks all the p values closer to one another. Use an "if" function
  
  X <- data$X
  index<- data$index
  n <- sum(index==1) # Assuming equal group sizes
  p <- par
  J <- max(index)
  D.vals <- rep(0,J)
  for (j in 1:J) {
    D.vals[j] <- binomial_nllh(X[index==j],p[j],N)
  }
  
  if (pen_type == "ToZero") {
    D <- sum(D.vals) + lambda*n*sum(p^2) # Penalty on squared p-values
  } else if (pen_type == "Closer") {
    p.matrix <- matrix(rep(p,J),nrow = J,byrow = TRUE)
    D <- sum(D.vals) + lambda*n*sum((p.matrix-t(p.matrix))^2)/2
  }
  return(D)
}

# Creates random permutation for test and validation data

perm_score <- function(dataset, lambda.grid, pen_type, N) {
  
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
  p.in <- c()
  for (a in 1:J) {
    p.in <- c(p.in,mean(dataset.test$X[dataset.test$index== a])/N) #Calculate mean of each subset
  }
  V <- rep(0,nL)
  for (b in 1:nL) {
    D <- optim(par = p.in, fn = pen_dist, data = dataset.test,lambda = lambda.grid[b],pen_type = pen_type, 
               N = N, method = "BFGS") #, lower = rep(0,J), upper = rep(1,J))
    # Validation score calculated on the validation data for each lambda value
    V[b] <- pen_dist(data = dataset.valid,lambda = 0, pen_type = pen_type,par = D$par, N = N )
  }
  return(V)
}

# Adding a siplation study component:
# (1) simulate a datasets
# (2) Find both MLE and penalized cross-validation est's
# (3) Store the estimators for each simulated dataset
# (4) Repeat this a large # of times so we can calculate
#     MSE (or RMSE) for the two methods


# Generate data 

sim.data <- function(n, trials, p.vals) {
  n.p <- length(p.vals)
  index <- rep(1,n)
  X <- rbinom(n, trials, p.vals)
  if (n.p > 1) {
    for (j in 2:n.p) {
      index <- c(index,rep(j,n))
      X <- c(X,rbinom(n, trials, p.vals))
    }
  }
  dataset <- data.frame(index = index, X = X)
  return(dataset)
}

lambdaMax = 50
nL = 15
M = 25 # Number of partitions
n.pvals <- 3
lambda.grid <- seq(0,lambdaMax,length.out = nL)
p.bin <- seq(0.1,0.9,length.out = n.pvals) # define lambda values for Binomial

MSEpen <- rep(1,n.pvals)
MSEmle <- rep(1,n.pvals)

n <- 50
n.var <- 10
pen_type <- "Closer"
trials <- 25

for (j in 1:n.pvals) {
  start_time <- Sys.time()
  K <- 10 # Number of datasets to simulate
  mle.est <- matrix(rep(0,n.var*K),nrow = K)
  pen.est <- matrix(rep(0,n.var*K),nrow = K)
  p.true <- seq(0.1,p.bin[j], length.out = n.var)
  
  
  for (k in 1:K) {
    dataset_simulation <- sim.data(n, trials, p.true)
    
    lambda.grid <- seq(0,lambdaMax,length.out = nL)
    V.store <- matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(dataset_simulation, lambda.grid, "Closer", N = trials)
      V.store[m,] <- V
    }
    
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    
    p.mle <- c()
    for (a in 1:n.var) {
      p.mle <- c(p.mle, mean(dataset_simulation$X[dataset_simulation$index == a])/trials) #Calculate mean of each subset
    }
    
    p.out <- optim(par = p.mle,fn = pen_dist, data = dataset_simulation, lambda = lambda.out, pen_type = pen_type, N = trials, method = "BFGS")$par
    
    mle.est[k,] <- p.mle
    pen.est[k,] <- p.out
  }
  mle.bias <- colMeans(mle.est) - p.true
  pen.bias <- colMeans(pen.est) - p.true
  mle.var <- apply(mle.est,2,var) # calculates var for each column #------
  pen.var <- apply(pen.est,2,var) # calculates var for each column
  mle.MSE <- mle.bias^2 + mle.var
  pen.MSE <- pen.bias^2 + pen.var
  MSEmle[j] <- sum(mle.MSE)
  MSEpen[j] <- sum(pen.MSE)
  end_time <- Sys.time()
  timer <- end_time - start_time
  print(j)
  print(timer)
}
plot(p.bin,log(MSEpen),type = "l")
lines(p.bin,log(MSEmle),col = "red") #-------

results <- list(delta = lambda.grid,
                MSEpen = MSEpen,
                MSEmle = MSEmle)

#name <- "delta_eq_spaced_closer_n50"
name <- "delta_one_nonzero_closer_n50"
save(results, file = paste(name, ".RData", sep=""))

# Very similar results: We don't know a good value for lambda max. Now we need to plot V.mean and lambda.grid to find 
# optimal values of lambda. Play around with p values and lambda vals.

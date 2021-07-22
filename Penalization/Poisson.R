# Poisson model: 
# log-likelihood function: l = sum(y)*log(mu) - n * mu

# Declaring log-likelihood function:
poisson_nllh <- function(data, par) {
  y <- data
  mu <- par
  #n <- nrow(y)  # get the sample size of the data
  llh <- sum(dpois(y, mu, log = TRUE))
  return(-llh)
}

pen_dist <- function(data,lambda,pen_type,par) {
  # Adapt this function so that if pen_type = "ToZero" then it shrinks all the p values closer to zero and if 
  # pen_type = "Closer" then it shrinks all the p values closer to one another. Use an "if" function
  
  X <- data$X
  index<- data$index
  n <- sum(index==1) # Assuming equal group sizes
  mu <- par
  J <- max(index)
  D.vals <- rep(0,J)
  for (j in 1:J) {
    D.vals[j] <- poisson_nllh(X[index==j],mu[j])
  }
  
  if (pen_type == "ToZero") {
    D <- sum(D.vals) + lambda*n*sum(mu^2) # Penalty on squared p-values
  } else if (pen_type == "Closer") {
    mu.matrix <- matrix(rep(mu,J),nrow = J,byrow = TRUE)
    D <- sum(D.vals) + lambda*n*sum((mu.matrix-t(mu.matrix))^2)/2
  }
  return(D)
}

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
    D <- optim(par = mu.in, fn = pen_dist, data = dataset.test,lambda = lambda.grid[b],pen_type = pen_type, method = "L-BFGS-B", lower = rep(0,J), upper = rep(Inf,J))
    # Validation score calculated on the validation data for each lambda value
    V[b] <- pen_dist(data = dataset.valid,lambda = 0, pen_type = pen_type,par = D$par )
  }
  return(V)
}

# Adding a siplation study component:
# (1) simulate a dataset
# (2) Find both MLE and penalized cross-validation est's
# (3) Store the estimators for each simulated dataset
# (4) Repeat this a large # of times so we can calculate
#     MSE (or RMSE) for the two methods


# Generate data 

sim.data <- function(N, mu.pois) {
  nmu.pois <- length(mu.pois)
  index <- rep(1,N)
  X <- rpois(N,mu.pois)
  if (nmu.pois > 1) {
    for (j in 2:nmu.pois) {
      index <- c(index,rep(j,N))
      X <- c(X,rpois(N,mu.pois))
    }
  }
  dataset <- data.frame(index = index, X = X)
  return(dataset)
}

lambdaMax = 1
nL = 15
M <- 25 # Number of partitions
n.muvals <- 3
lambda.grid <- seq(0,lambdaMax,length.out = nL)
mu.pois <- seq(0.5,2,length.out = n.muvals) # define lambda values for Poisson 

# dataset_simulation <- sim.data(20,c(1,1.1,1.2))
# 
# V.store <- matrix(rep(0,M*nL),nrow = M)
# for (m in 1:M) {
#   V <- perm_score(dataset_simulation, lambda.grid, "Closer")
#   V.store[m,] <- V
# }
# V.mean <- colMeans(V.store)
# plot(lambda.grid, V.mean, type = "l")


MSEpen <- rep(1,n.muvals)
MSEmle <- rep(1,n.muvals)

N <- 50
n.var <- 10
pen_type <- "Closer"
for (j in 1:n.muvals) {
  start_time <- Sys.time()
  K <- 10 # Number of datasets to simulate
  mle.est <- matrix(rep(0,n.var*K),nrow = K)
  pen.est <- matrix(rep(0,n.var*K),nrow = K)
  mu.true <- seq(0.5,mu.pois[j], length.out = n.var)
  
  
  for (k in 1:K) {
    dataset_simulation <- sim.data(N,mu.true)
    
    lambda.grid <- seq(0,lambdaMax,length.out = nL)
    V.store <- matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(dataset_simulation, lambda.grid, "Closer")
      V.store[m,] <- V
    }
    
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    
    mu.mle <- c()
    for (a in 1:n.var) {
      mu.mle <- c(mu.mle, mean(dataset_simulation$X[dataset_simulation$index == a])) #Calculate mean of each subset
    }

    mu.out <- optim(par = mu.mle,fn = pen_dist, data = dataset_simulation, lambda = lambda.out, pen_type = pen_type, method = "L-BFGS-B", lower = rep(0,n.var), upper = rep(Inf, n.var))$par
    
    mle.est[k,] <- mu.mle
    pen.est[k,] <- mu.out
  }
  mle.bias <- colMeans(mle.est) - mu.true
  pen.bias <- colMeans(pen.est) - mu.true
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
plot(mu.pois,MSEpen,type = "l")
lines(mu.pois,MSEmle,col = "red") #-------

results <- list(delta = lambda.grid,
                MSEpen = MSEpen,
                MSEmle = MSEmle)

#name <- "delta_eq_spaced_closer_n50"
name <- "delta_one_nonzero_closer_n50"
save(results, file = paste(name, ".RData", sep=""))



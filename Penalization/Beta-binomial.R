library(extraDistr)
#n: how many obs you have per variable
#N: passage size, number of words in the passage

betabin_nllh <- function(data, N, par) {
  x <- data
  alpha <- par[1]
  beta <- par[2]
  bb_llh <- sum(dbbinom(x,N,alpha, beta, log = TRUE))
  return(-bb_llh)
}
pen_dist <- function(data,lambda,pen_type,par,N) {
  X <- data$X
  index<- data$index
  n <- sum(index==1) # Assuming equal group sizes
  pen_type <- "Closer"
  J <- max(index)
  alpha <- par[1:J]
  beta <- par[(J+1):(2*J)]
  p <- alpha/(alpha + beta)
  D.vals <- rep(0,J)
  for (j in 1:J) {
    D.vals[j] <- betabin_nllh(X[index==j],N,c(alpha[j],beta[j]))
  }
  p.matrix <- matrix(rep(p,J),nrow = J,byrow = TRUE)
  D <- sum(D.vals) + lambda*n*sum((p.matrix-t(p.matrix))^2)/2
  
  return(D)
}

# Creates random permutation for test and validation data

perm_score <- function(dataset, lambda.grid, pen_type, N) {
  
  index <- dataset$index
  X <- dataset$X
  J <- max(index)
  n <- length(X[index==1]) # Number of values
  m <- floor(0.75*n)
  pen_type <- "Closer"
  
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
  alpha.in <- c()
  beta.in <- c()
  for (a in 1:J) {
    MLE_bb <- optim(c(1,1), betabin_nllh, data = dataset.test$X[dataset.test$index== a], N = N )$par
    alpha.in <- c(alpha.in, MLE_bb[1])
    beta.in <- c(beta.in, MLE_bb[2])
  }
  V <- rep(0,nL)
  for (b in 1:nL) {
    D <- optim(c(alpha.in,beta.in), fn = pen_dist, data = dataset.test,lambda = lambda.grid[b],pen_type = pen_type, 
               N = N, method = "BFGS") 
    # Validation score calculated on the validation data for each lambda value
    V[b] <- pen_dist(data = dataset.valid,lambda = 0, pen_type = pen_type,par = D$par, N = N )
  }
  return(V)
}
# Generate data 

sim.data <- function(n, N, alpha.vals, beta.vals) {
  n.par <- length(alpha.vals)
  index <- rep(1,n)
  X <- rbbinom(n,N,alpha.vals,beta.vals)
  if (n.par > 1) {
    for (j in 2:n.par) {
      index <- c(index,rep(j,n))
      X <- c(X,rbbinom(n, N, alpha.vals, beta.vals))
    }
  }
  dataset <- data.frame(index = index, X = X)
  return(dataset)
}

lambdaMax = 10
nL = 5
M = 10 # Number of partitions
n.pvals <- 2
lambda.grid <- seq(0,lambdaMax,length.out = nL)

alpha.bin <- seq(0.1,0.9,length.out = n.pvals) 
beta.bin <- seq(0.1,0.9,length.out = n.pvals)

MSEpen <- rep(1,n.pvals)
MSEmle <- rep(1,n.pvals)

n <- 20
n.var <- 5
pen_type <- "Closer"
N <- 25

for (j in 1:n.pvals) {
  start_time <- Sys.time()
  K <- 3 # Number of datasets to simulate
  mle.est <- matrix(rep(0,n.var*K*2),nrow = K)
  pen.est <- matrix(rep(0,n.var*K*2),nrow = K)
  alpha.true <- seq(0.1,alpha.bin[j], length.out = n.var)
  beta.true <- seq(0.1, beta.bin[j], length.out = n.var)
  p.true <- c(alpha.true, beta.true)
  
  for (k in 1:K) {
    dataset_simulation <- sim.data(n, N, alpha.true, beta.true)
    
    lambda.grid <- seq(0,lambdaMax,length.out = nL)
    V.store <- matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(dataset_simulation, lambda.grid, pen_type, N)
      V.store[m,] <- V
    }
    
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    
    alpha.mle <- c()
    beta.mle <- c()
    for (a in 1:n.var) {
      MLE_bb <- optim(c(1,1), betabin_nllh, data = dataset_simulation$X[dataset_simulation$index== a], N = N )$par
      alpha.mle <- c(alpha.mle, MLE_bb[1]) 
      beta.mle <- c(beta.mle, MLE_bb[2])
    }
    par.out <- optim(c(alpha.mle, beta.mle),fn = pen_dist, data = dataset_simulation, lambda = lambda.out, pen_type = pen_type, 
                     N = N, method = "BFGS")$par
    
    mle.est[k,] <- c(alpha.mle, beta.mle)
    pen.est[k,] <- par.out
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
p.bin <- alpha.bin/(alpha.bin + beta.bin)
plot(p.bin,log(MSEpen),type = "l")
lines(p.bin,log(MSEmle),col = "red") #-------

results <- list(delta = lambda.grid,
                MSEpen = MSEpen,
                MSEmle = MSEmle)

#name <- "delta_eq_spaced_closer_n50"
name <- "delta_one_nonzero_closer_n50"
save(results, file = paste(name, ".RData", sep=""))
# Report: Beta - binomial, using equations when possible
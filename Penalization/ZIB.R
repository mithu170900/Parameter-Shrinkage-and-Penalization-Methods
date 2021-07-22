library(extraDistr)
#n: how many obs you have per variable
#N: passage size, number of words in the passage

zifbin_nllh <- function(data, par, N) {
  y <- data
  gam <- par[1]
  p <- par[2]
  #llh0 <- log(gam+(1-gam)*dbinom(y,N,p))*(y==0) # We don't have to change this setup,
  #llh1 <- log((1-gam)*dbinom(y,N,p))*(y>0) # but there is a dzib function which can be used for llh
  #zifbin_llh <- sum(llh0)+sum(llh1)
  zifbin_llh <- sum(dzib(y, N, prob = p, pi = gam, log = TRUE))
  return(-zifbin_llh)
}

pen_dist <- function(data,lambda,pen_type,par,N) { # Can remove the pen_type argument here
  X <- data$X
  index<- data$index
  n <- sum(index==1) # Assuming equal group sizes
  pen_type <- "Closer"
  J <- max(index)
  gam <- par[1:J]
  p <- par[(J+1):(2*J)]
  q <- (1-gam)*p
  D.vals <- rep(0,J)
  for (j in 1:J) {
    D.vals[j] <- zifbin_nllh(X[index==j],c(gam[j],p[j]),N)
  }
  q.matrix <- matrix(rep(q,J),nrow = J,byrow = TRUE)
  D <- sum(D.vals) + lambda*n*sum((q.matrix-t(q.matrix))^2)/2
  
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
  gam.in <- c()
  p.in <- c()
  for (a in 1:J) {
    p.try <- mean(dataset.test$X[dataset.test$index== a])/N
    MLE_zeroBinomial <- optim(c(0.01, p.try),zifbin_nllh, data = dataset.test$X[dataset.test$index== a], N = N, method = "L-BFGS-B", 
                              lower = c(0.00000001,0.00000001), upper = c(0.9999999,0.9999999))$par
    gam.in <- c(gam.in, MLE_zeroBinomial[1])
    p.in <- c(p.in, MLE_zeroBinomial[2])
  }
  V <- rep(0,nL)
  for (b in 1:nL) {
    D <- optim(c(gam.in,p.in), fn = pen_dist, data = dataset.test,lambda = lambda.grid[b],pen_type = pen_type, 
               N = N, method = "L-BFGS-B", lower = c(0.00000001,0.00000001), upper = c(0.9999999,0.9999999)) 
    # Validation score calculated on the validation data for each lambda value
    V[b] <- pen_dist(data = dataset.valid,lambda = 0, pen_type = pen_type,par = D$par, N = N )
  }
  return(V)
}

# Generate data 

sim.data <- function(n, N, gam.vals, p.vals) {
  # n is the number of observations per variable
  # N is the binomial sample size
  # gam.vals is an array of zero probabilities (inflation component)
  # p.vals is an array of binomial success probabilities 
  n.par <- length(gam.vals)
  index <- rep(1,n)
  X <- rzib(n,N,p.vals[1],gam.vals[1]) # Flipped order
  if (n.par > 1) {
    for (j in 2:n.par) {
      index <- c(index,rep(j,n))
      X <- c(X,rzib(n, N, p.vals[j], gam.vals[j])) # Added indices
    }
  }
  dataset <- data.frame(index = index, X = X)
  return(dataset)
}

lambdaMax = 2
nL = 50
M = 20 # Number of partitions
n.pvals <- 3
lambda.grid <- seq(0,lambdaMax,length.out = nL)

gam.bin <- seq(0.05,0.25,length.out = n.pvals) 
p.bin <- seq(0.05,0.25,length.out = n.pvals)

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
  gam.true <- seq(0.01,gam.bin[j], length.out = n.var)
  p.true <- seq(0.01, p.bin[j], length.out = n.var)
  q.true <- c(gam.true, p.true)
  
  for (k in 1:K) {
    dataset_simulation <- sim.data(n, N, gam.true, p.true)
    J <- max(dataset_simulation$index)
    
    lambda.grid <- seq(0,lambdaMax,length.out = nL)
    V.store <- matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(dataset_simulation, lambda.grid, pen_type, N)
      V.store[m,] <- V
    }
    
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    
    gam.mle <- c()
    p.mle <- c()
    for (a in 1:n.var) {
      p.try <- mean(dataset_simulation$X[dataset_simulation$index== a])/N
      MLE_zeroBinomial <- optim(c(0.01, p.try), zifbin_nllh, data = dataset_simulation$X[dataset_simulation$index==a], 
                          N = N, method = "L-BFGS-B", lower = c(0.00000001,0.00000001), upper = c(0.9999999,0.9999999))$par
      gam.mle <- c(gam.mle, MLE_zeroBinomial[1]) 
      p.mle <- c(p.mle, MLE_zeroBinomial[2])
    }
    par.out <- optim(c(gam.mle, p.mle),fn = pen_dist, data = dataset_simulation, lambda = lambda.out, pen_type = pen_type, 
                     N = N, method = "L-BFGS-B", lower = c(rep(0.00000001,J),rep(0.00000001,J)), upper = c(rep(0.9999999,J),rep(0.9999999,J)))$par
    mle.est[k,] <- c(gam.mle, p.mle)
    pen.est[k,] <- par.out
  }
  mle.bias <- colMeans(mle.est) - q.true
  pen.bias <- colMeans(pen.est) - q.true
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
q.bin <- (1-gam.bin)*p.bin
plot(q.bin,MSEpen,type = "l")
lines(q.bin,MSEmle,col = "red") #-------

results <- list(delta = lambda.grid,
                MSEpen = MSEpen,
                MSEmle = MSEmle)
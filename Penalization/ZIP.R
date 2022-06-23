library(extraDistr)
zifpois_nllh <- function(data,par) {
  y <- data
  gam <- par[1] # P(Observing 0)
  mu <- par[2] # Mean of Poisson distribution
  #llh0 <- log(gam+(1-gam)*dpois(y, mu))*(y==0) 
  #llh1 <- log((1-gam)*dpois(y, mu))*(y>0)
  #zifpois_llh <- sum(llh0)+sum(llh1)
  zifpois_llh <- sum(dzip(y, mu,gam,log = TRUE))
  return(-zifpois_llh)
}

pen_dist <- function(data,lambda,pen_type,par) {
  X <- data$X
  index<- data$index
  n <- sum(index==1) # Assuming equal group sizes
  pen_type <- "Closer"
  J <- max(index)
  gam <- par[1:J]
  mu <- par[(J+1):(2*J)]
  q <- (1-gam)*mu
  D.vals <- rep(0,J)
  for (j in 1:J) {
    D.vals[j] <- zifpois_nllh(X[index==j],c(gam[j],mu[j]))
  }
  q.matrix <- matrix(rep(q,J),nrow = J,byrow = TRUE)
  D <- sum(D.vals) + lambda*n*sum((q.matrix-t(q.matrix))^2)/2
  
  return(D)
}

# Creates random permutation for test and validation data

perm_score <- function(dataset, lambda.grid, pen_type) {
  
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
  mu.in <- c()
  for (a in 1:J) {
    mu.try <- mean(dataset.test$X[dataset.test$index== a])
    MLE_zeroPoisson <- optim(c(0.5,mu.try), fn = zifpois_nllh, data = dataset.test$X[dataset.test$index== a])$par
    gam.in <- c(gam.in, MLE_zeroPoisson[1])
    mu.in <- c(mu.in, MLE_zeroPoisson[2])
  }
  V <- rep(0,nL)
  for (b in 1:nL) {
    D <- optim(c(gam.in,mu.in), fn = pen_dist, data = dataset.test,lambda = lambda.grid[b],pen_type = pen_type) 
    # Validation score calculated on the validation data for each lambda value
    V[b] <- pen_dist(data = dataset.valid,lambda = 0, pen_type = pen_type,par = D$par)
  }
  return(V)
}
# Generate data 

sim.data <- function(n, gam.vals, mu.vals) {
  n.par <- length(mu.vals)
  index <- rep(1,n)
  X <- rzip(n,mu.vals[1],gam.vals[1])
  if (n.par > 1) {
    for (j in 2:n.par) {
      index <- c(index,rep(j,n))
      X <- c(X,rzip(n, mu.vals[j],gam.vals[j]))
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

gam.bin <- seq(0.8,1,length.out = n.pvals) 
mu.bin <- seq(0.5,2,length.out = n.pvals)

MSEpen <- rep(1,n.pvals)
MSEmle <- rep(1,n.pvals)

n <- 20
n.var <- 5
pen_type <- "Closer"

for (j in 1:n.pvals) {
  start_time <- Sys.time()
  K <- 3 # Number of datasets to simulate
  mle.est <- matrix(rep(0,n.var*K*2),nrow = K)
  pen.est <- matrix(rep(0,n.var*K*2),nrow = K)
  gam.true <- seq(0.8,gam.bin[j], length.out = n.var)
  mu.true <- seq(0.5, mu.bin[j], length.out = n.var)
  q.true <- c(gam.true, mu.true)
  
  for (k in 1:K) {
    dataset_simulation <- sim.data(n, gam.true, mu.true)
    
    lambda.grid <- seq(0,lambdaMax,length.out = nL)
    V.store <- matrix(rep(0,M*nL),nrow = M)
    for (m in 1:M) {
      V <- perm_score(dataset_simulation, lambda.grid, pen_type)
      V.store[m,] <- V
    }
    
    V.mean <- colMeans(V.store)
    lambda.out <- lambda.grid[V.mean==min(V.mean)]
    
    gam.mle <- c()
    mu.mle <- c()
    for (a in 1:n.var) {
      mu.try2 <- mean(dataset_simulation$X[dataset_simulation$index== a])
      MLE_zeroPoisson <- optim(c(0.5,mu.try2), fn = zifpois_nllh, data = mean(dataset_simulation$X[dataset_simulation$index== a]))$par
      gam.mle <- c(gam.mle, MLE_zeroPoisson[1]) 
      mu.mle <- c(mu.mle, MLE_zeroPoisson[2])
    }
    par.out <- optim(c(gam.mle, mu.mle),fn = pen_dist, data = dataset_simulation, lambda = lambda.out, pen_type = pen_type)$par
    mle.est[k,] <- c(gam.mle, mu.mle)
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
q.bin <- (1-gam.bin)*mu.bin
plot(q.bin,MSEpen,type = "l")
lines(q.bin,MSEmle,col = "red") #-------

results <- list(delta = lambda.grid,
                MSEpen = MSEpen,
                MSEmle = MSEmle)
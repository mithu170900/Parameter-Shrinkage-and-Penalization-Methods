# Poisson model: 
# log-likelihood function: l = sum(y)*log(mu) - n * mu
Count.data <- passages$Count1
n.data <- passages$N1
# Declaring log-likelihood function:
poisson_nllh <- function(data, par) {
  y <- data
  mu <- par
  #n <- nrow(y)  # get the sample size of the data
  llh <- sum(dpois(y, mu, log = TRUE))
  return(-llh)
}
# Maximum likelihood function:
MLE_poisson <- optim(1,poisson_nllh, data = Count.data)

# Calculating AIC:
aic_poisson <- 2 * poisson_nllh(Count.data, MLE_poisson$par) + 2 * 2

# Graph pmf:

# Binomial model:
# log-likelihood function: l = x*log p + (n−x)*log(1−p)
# x is the number of successes (x0) and n-x is the number of failures (x1)

#Declaring log-likelihood function:
binomial_nllh <- function(data, par, n) {
  y <- data
  p <- par
  bi_llh <- sum(dbinom(y, n, p, log = TRUE))
  return(-bi_llh)
}
#Maximum likelihood function:
MLE_binomial <- optim(0.5,binomial_nllh, data=Count.data, n = n.data , method = "Brent", lower = 0, upper = 1)

# Calculating AIC:
count_sucesses <- length(which(Count.data == 0))
prob <- count_sucesses/n.data 
aic_binomial <- 2 * binomial_nllh(Count.data, MLE_binomial$par, n.data ) + 2 * 1

# Beta-binomial:
# No closed-form likelihood function. Therefore, we'll be using a built-in function 
#    in R to find MLE, bb

# Estimation of parameters:
# x = number of successes and n = number of trials
library(extraDistr)
betabin_nllh <- function(data, par, n) {
  x <- data
  alpha <- par[1]
  beta <- par[2]
  bb_llh <- sum(dbbinom(x,n,alpha, beta, log = TRUE))
  return(-bb_llh)
}
#Maximum likelihood approach:
MLE_bb <- optim(c(1,1), betabin_nllh, data = Count.data, n = n.data )

#Calculating AIC:
aic_bb <- 2 * betabin_nllh(Count.data, MLE_bb$par, n.data ) + 2 * 2

# Plotting empirical probabilities:
ecdf(Count.data)
vector1 <- c(unique(Count.data))
prob1 <- c((vector1/n.data ) * 100) 
plot(x = unique(Count.data), y = prob1 , xlab="Count 1 Values", ylab = "Percentages", verticals = 'True', col.01line = "gray70", pch = 20)

grid <- seq(0,max(Count.data)+2, by = 1)
n.grid <- length(grid)
prob1 <- rep(0,n.grid)
for (j in 1:length(grid)) {
  prob1[j] <- mean(Count.data==grid[j])
}
plot(x = grid, y = prob1 , type = "h", xlab="Count 1 Values", ylab = "Probabilities", ylim=c(0, 1))
points(grid,prob1)

# Zero-inflated binomial distribution:
zifbin_nllh <- function(data, par, n) {
  y <- data
  gam <- par[1]
  p <- par[2]
  llh0 <- log(gam+(1-gam)*dbinom(y,n,p))*(y==0) 
  llh1 <- log((1-gam)*dbinom(y,n,p))*(y>0)
  zifbin_llh <- sum(llh0)+sum(llh1)
  return(-zifbin_llh)
}
MLE_zeroBinomial <- optim(c(0.5, 0.5),fn = zifbin_nllh, data=Count.data, n = n.data, method = "L-BFGS-B", lower = c(0.00000001,0.00000001), upper = c(0.9999999,0.9999999))
aic_zeroBinomial <- 2 * zifbin_nllh(Count.data, MLE_zeroBinomial$par, n.data ) + 2 * 2

# Zero-inflated Poisson distribution:
zifpois_nllh <- function(par, data) {
  y <- data
  mu <- par[1]
  gam <- par[2]
  llh0 <- log(gam+(1-gam)*dpois(y, mu))*(y==0) 
  llh1 <- log((1-gam)*dpois(y, mu))*(y>0)
  zifpois_llh <- sum(llh0)+sum(llh1)
  return(-zifpois_llh)
}
MLE_zeroPoisson <- optim(c(1,0.5), fn = zifpois_nllh, data = Count.data)
aic_zeroPoisson <- 2 * zifpois_nllh( MLE_zeroPoisson$par, Count.data) + 2 * 2

# Plotting model-based probabilities for Count1:

# Poisson distribution for Count1:
prob_poisson <- dpois(grid, MLE_poisson$par)
plot(x = grid, y = prob_poisson, type = "h", xlab = "Count 1 Values", ylab = "Poisson probabilities", ylim=c(0, 1))
points(grid,prob_poisson)

# Binomial distribution for Count1:
prob_binom <- dbinom(grid, n.data, MLE_binomial$par)
plot(x = grid, y = prob_binom, type = "h", xlab = "Count 1 Values", ylab = "Binomial probabilities", ylim=c(0, 1))
points(grid,prob_binom)

# Beta-binomial distribution for Count1:
prob_bb <- dbbinom(grid, n.data , alpha = MLE_bb$par[1], beta = MLE_bb$par[2])
plot(x = grid, y = prob_bb, type = "h", xlab = "Count 1 Values", 
     ylab = "ZBeta-binomial probabilities", ylim=c(0, 1))
points(grid,prob_bb)

# Zero-inflated Poisson Distribution:
prob.zeroPoisson <- (MLE_zeroPoisson$par[2]) + ((1 - MLE_zerobBinomial$par[2] ) * dpois(grid[1], MLE_zeroPoisson$par[1]))
prob_zeroPoisson <- (1 - MLE_zeroPoisson$par[2] ) * dpois(grid[2:n.grid], MLE_zeroPoisson$par[1])
plot(x = grid, y = c(prob.zeroPoisson,prob_zeroPoisson), type = "h", xlab = "Count 1 Values", ylab = "Zero-inflated Poisson probabilities", ylim=c(0, 1))
points(grid[1:n.grid],c(prob.zeroPoisson,prob_zeroPoisson))

# Zero-inflated Binomial Distribution:
prob.zero <- (MLE_zeroBinomial$par[1]) + ((1 - MLE_zeroBinomial$par[1] ) * dbinom(grid[1], MLE_zeroBinomial$par[2],size=n.data))
prob_zeroBinom <- (1 - MLE_zeroBinomial$par[1] ) * dbinom(grid[2:n.grid], MLE_zeroBinomial$par[2],size=n.data)
plot(x = grid, y = c(prob.zero,prob_zeroBinom), type = "h", xlab = "Count 1 Values", ylab = "Zero-inflated Binomial probabilities", ylim=c(0, 1))
points(grid[1:n.grid],c(prob.zero,prob_zeroBinom))



par(mfrow=c(2,3))
plot(x = grid, y = prob1 , type = "h", xlab="Count 1 Values", ylab = "Empirical Probabilities", ylim=c(0, 1))
points(grid,prob1)
plot(x = grid, y = prob_poisson, type = "h", xlab = "Count 1 Values", ylab = "Poisson probabilities", ylim=c(0, 1))
points(grid,prob_poisson)
plot(x = grid, y = prob_binom, type = "h", xlab = "Count 1 Values", ylab = "Binomial probabilities", ylim=c(0, 1))
points(grid,prob_binom)
plot(x = grid, y = prob_bb, type = "h", xlab = "Count 1 Values", ylab = "Beta-binomial probabilities", ylim=c(0, 1))
points(grid,prob_bb)
plot(x = grid, y = c(prob.zeroPoisson,prob_zeroPoisson), type = "h", xlab = "Count 1 Values", ylab = "Zero-inflated Poisson probabilities", ylim=c(0, 1))
points(grid[1:n.grid],c(prob.zeroPoisson,prob_zeroPoisson))
plot(x = grid, y = c(prob.zero,prob_zeroBinom), type = "h", xlab = "Count 1 Values", ylab = "Zero-inflated Binomial probabilities", ylim=c(0, 1))
points(grid[1:n.grid],c(prob.zero,prob_zeroBinom))

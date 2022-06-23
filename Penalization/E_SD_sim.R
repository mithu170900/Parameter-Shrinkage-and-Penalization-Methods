# Binomial (this code is "confirmatory" since you already did calcs by hand)
K <- 1000000 # Number of simulated data points
alpha <- 10
beta <- 10
a <- 0.31
b <- 0.35
Count <- rep(0,K)
for (j in 1:K) {
  p.sim <- a + (b-a)*rbeta(1,alpha,beta) # Simulate beta rand #, rescales to (a,b)
  Count[j] <- rbinom(1,size=40,prob = p.sim)
}
c(mean(Count),sd(Count))

# Code for zero-inflated binomial
library(extraDistr)
K <- 1000000
alpha <- 10
beta <- 10
a1 <- 0.01
b1 <- 0.05
a2 <- 0.10
b2 <- 0.14
Count <- rep(0,K)
for (j in 1:K) {
  p.sim <- a1 + (b1-a1)*rbeta(1,alpha,beta)
  gamma.sim <- a2 + (b2-a2)*rbeta(1,alpha,beta)
  Count[j] <- rzib(1,size=40,prob = p.sim, pi = gamma.sim)
}
c(mean(Count),sd(Count))

# Code for beta-binomial 
K <- 1000000
alpha <- 10
beta <- 10
# pi values
a1 <- 0.05
b1 <- 0.06
# nu values
a2 <- 2
b2 <- 10
N <- 40
Count <- rep(0,K)
for (j in 1:K) {
  p.sim <- a1 + (b1-a1)*rbeta(1,alpha,beta)
  nu.sim <- a2 + (b2-a2)*rbeta(1,alpha,beta)
  Count[j] <- rbbinom(1, size = 40, alpha = (p.sim*(N - 1))/(nu.sim -1) - p.sim, beta= ((1-p.sim)*(N-1)/(nu.sim-1)) - p.sim)
}
c(mean(Count),sd(Count))

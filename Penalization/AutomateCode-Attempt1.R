#-----------------------------------------------------------------------------------------------
# 1. Spacing delta: Make them equally spaced
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
# Sample size 
n <- 50
# What is the mean for many variables?
# Set up data frame for data
index <- c(rep(1,n))
X <- rnorm(n, 0, 1)
for (j in 2:10) {
  
}
# Set up mu values for parameters
# # Penalized distance function using for-loops for data and sample size
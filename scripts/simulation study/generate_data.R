###---------------------------------------------###
###   Author: Juliana Freitas de Mello e Silva  ###
###---------------------------------------------###

###--------------------------------
###   Removing already allocated objects
###--------------------------------

rm(list = ls())

###--------------------------------
###   Required packages
###--------------------------------

require(MASS)

###--------------------------------
###   Setting seed
###--------------------------------

set.seed(123456789)

###--------------------------------
###   Defining work directory
###--------------------------------

directory <- paste0(getwd(), "/") ; setwd(directory)
dir.create(paste0(directory, "data"))

###--------------------------------###--------------------------------###

###--------------------------------
###   Generating data
###--------------------------------

### Number of data sets, sample size and number of measurements (same)
L <- 50 ; n <- 30 
J_i <- sample(x = 3:8, size = n, replace = TRUE)

### Mean function (supposing that it is known)
mu <- function(t){
  f <- 10 + 10*sin(2*pi*t)
  return(f)
}

### Parameters of a Beta distribution
a <- 3 ; b <- 3

### Variance
sigma <- 3 ; cov <- 6

###--------------------------------###--------------------------------###

### Directory for each configuration
results_directory <- paste0(directory, "data/beta", a, b, "/")
dir.create(results_directory)

list_par<- list(a = a, b = b, mu = mu, sigma = sigma, cov = cov, L = L, n = n, J_i = J_i)
save(list_par, file = paste0(directory, "list_parameters_beta", a, b, ".R"))


for(l in 1:L){
  
  ### Measurement times
  t_ij_wide <- matrix(data = NA, nrow = n, ncol = max(J_i))
  for(i in 1:n){
    t_ij_wide[i,1:J_i[i]] <- c(0, sort(rbeta(n = J_i[i] - 1, shape1 = a, shape2 = b)))
  }
  t_ij <- na.omit(as.numeric(t(t_ij_wide)))
  
  ### Observations
  W <- matrix(data = NA, nrow = sum(J_i), ncol = 1)
  rownames(W) <- paste0(unlist(t(sapply(1:n, function(x)(rep(x = x, times = J_i[x]))))), "_",
                       unlist(t(sapply(1:n, function(x)(1:J_i[x])))))
  W_wide <- matrix(data = NA, nrow = n, ncol = max(J_i)) ; colnames(W_wide)<- 1:max(J_i) ; rownames(W_wide)<- 1:n
  
  for(i in 1:n){
    Sigma <- matrix(data = cov, ncol = J_i[i], nrow = J_i[i]) ; diag(Sigma) <- sigma^2
    W_wide[i,1:J_i[i]] <- mvrnorm(n = 1,  mu = mu(as.numeric(t_ij_wide[i,1:J_i[i]])), Sigma = Sigma)
  }
  W <- na.omit(as.numeric(t(W_wide)))
  
  ### Observatios on the long format
  id <- sort(as.numeric(unlist(t(sapply(1:n, function(x)(rep(x = x, times = J_i[x])))))))
  
  data <- data.frame(id = id, w_obs = W, t_ij = t_ij)
  
  ### Saving datasets
  write.table(x = data, file = paste0(results_directory, "/data", l, ".txt"))
  
  
}







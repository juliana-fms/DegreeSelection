mode <- function (x) {
  if (!is.numeric(x)) {
    return(NA)
    stop('Not numeric x vector')
  } else {
    fit <- density(x)
    fit$x[which.max(fit$y)]
  }
}

###----------------------------------------------

hmean<- function(x){
  return(length(x) / sum(1 / x))
}

### --------------------------------------------------------

CPO <- function(loglik, n)
{ 
  lik <- apply(loglik, 2, exp)
  CPO <- apply(lik, 2, hmean)
  LPML <- sum(log(CPO))
  aLPML <- mean(log(CPO))
  return(c(LPML, aLPML)) 
}

### --------------------------------------------------------
DIC <- function(loglik, n)
{
  D <- apply( -2*loglik, 1, sum);
  pD <- 0.5*var(D)
  DIC <- mean(D) + pD;
  return( matrix(c(DIC, pD), ncol=2) ) 
}

### --------------------------------------------------------
WAIC <- function(loglik, n)
{
  lpd <- sum( log( apply(exp(loglik), 2, mean) ) )
  pD <- sum( apply(loglik,  2, var) )
  WAIC <- lpd - pD
  return( matrix(c(WAIC, pD),ncol=2) ) 
}

#---------------------------------


BP_Basis<- function(t, m){
  return(t(sapply(t, function(x)(sapply(1:m, function(k)(dbeta(x = x, shape1 = k, shape2 = m - k + 1))) / m))))
}

###----------------------------

degree_elevation_r1<- function(xi){
  ###-------------------------------------------------------###
  ### Function to use degree elevation property when r = 1  ###
  ### xi -> matrix of Bernstein coefficients of degree m    ###
  ### Returns a matrix with (m+1) Bernstein coefficients    ###
  ###-------------------------------------------------------###
  m<- ncol(xi)
  xi_elevation<- matrix(data = NA, nrow = nrow(xi), ncol = m + 1)
  xi_elevation[,1]<- xi[,1] ; xi_elevation[,m+1]<- xi[,m]
  xi_elevation[,2:m]<- sapply(2:m, function(k)(((k / (m+1))*xi[,k]) + ((1 - (k / (m+1)))*xi[,k])))
  return(xi_elevation)
}





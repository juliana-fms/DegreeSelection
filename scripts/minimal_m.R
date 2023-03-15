###-----------------------------------------------------------------------###
###   Script to find minimum degree m of Bernstein Polynomials            ###
###   Author: Juliana Freitas de Mello e Silva                            ###
###                                                                       ###
###   M is the random variable representing the minimum degree that is    ###
###   necessary to capture the change point in the interval (u(1), u(2))  ###
###-----------------------------------------------------------------------###

###-------------------------------
###   Removing objects
###-------------------------------

rm(list = ls())

###-------------------------------
###   Setting seed
###-------------------------------

set.seed(123456789)

###-------------------------------
###   Requiring packages
###-------------------------------


###-------------------------------###-------------------------------###


### Probability function of random variable 

p_M <- function(m, theta1, theta2){
  
  prob <- ((pbeta(q = (m-2)/(m-1), shape1 = theta1, shape2 = theta2) -
              pbeta(q = 1/(m-1), shape1 = theta1, shape2 = theta2))^2) -
    (pbeta(q = (m-3)/(m-2), shape1 = theta1, shape2 = theta2) -
       pbeta(q = 1/(m-2), shape1 = theta1, shape2 = theta2))^2
  
  return(prob)
  
}

### Cumulative distribution function of random variable M

F_M <- function(m, theta1, theta2){
  
  F_aux <- (pbeta(q = (m-2)/(m-1), shape1 = theta1, shape2 = theta2) -
          pbeta(q = 1/(m-1), shape1 = theta1, shape2 = theta2))^2
  
  return(F_aux)
  
}


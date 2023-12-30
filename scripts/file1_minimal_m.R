###---------------------------------------------###
###   File #1: minimum degree                   ###
###   Bernstein Polynomials                     ###
###   Author: Juliana Freitas de Mello e Silva  ###
###   Email: julianafms27@gmail.com             ###
###---------------------------------------------###

###--------------------------------
###   Removing objects
###--------------------------------

rm(list = ls())

###--------------------------------###--------------------------------###

###--------------------------------
###   Probability function of M
###--------------------------------

p_M_Beta <- function(m, theta1, theta2){
  
  if(any(m < 4)){
    stop("m must be greater or equal 4")
    } else {
      if (any(theta1 <= 0, theta2 <= 0)){
        stop("theta1 > 0 and theta2 > 0") 
      }
    }
  
  prob <- ((pbeta(q = (m-2)/(m-1), shape1 = theta1, shape2 = theta2) -
              pbeta(q = 1/(m-1), shape1 = theta1, shape2 = theta2))^2) -
    (pbeta(q = (m-3)/(m-2), shape1 = theta1, shape2 = theta2) -
       pbeta(q = 1/(m-2), shape1 = theta1, shape2 = theta2))^2
  return(prob)
}

# p_M_Beta(m = 4:20, theta1 = 0, theta2 = 1)

###--------------------------------
###   Cumulative distribution of M given U1 and U2
###--------------------------------

F_M_Beta <- function(m, theta1, theta2){
  
  if(any(m < 4)){
    stop("m must be greater or equal 4")
  } else {
    if (any(theta1 <= 0, theta2 <= 0)){
      stop("theta1 > 0 and theta2 > 0") 
    }
  }
  
  prob <- (pbeta(q = (m-2)/(m-1), shape1 = theta1, shape2 = theta2) -
             pbeta(q = 1/(m-1), shape1 = theta1, shape2 = theta2))^2
  return(prob)
}

# F_M_Beta(m = 8, theta1 = 11, theta2 = 34)


### Next, you can run models


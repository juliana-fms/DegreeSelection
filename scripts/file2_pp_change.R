###---------------------------------------------------------###
###   File #2: posterior probability of changing concavity  ###
###   Bernstein Polynomials                                 ###
###   Author: Juliana Freitas de Mello e Silva              ###
###   Email: julianafms27@gmail.com                         ###
###---------------------------------------------------------###

###--------------------------------
###   Removing objects already allocated
###--------------------------------

rm(list = ls())

###--------------------------------
###   Loading packages
###--------------------------------


###--------------------------------
###   Setting seed
###--------------------------------

set.seed(123456789)

###--------------------------------###--------------------------------###

# sample <- read.table(file = "data\\posterior_mu.txt")

###--------------------------------
###   Function that calculates the posterior probability of change in the
###   estimated curve
###--------------------------------

pp_change <- function(sample = sample, m = NULL, S = NULL){
  
  ### sample: posterior sample of the vector of coefficients of dimensions Sxm
  
  ### Degree of the BP
  m <- ncol(sample)
  
  ### Size of posterior chain
  S <- nrow(sample)
  
  ###  Calculating differences between xi[k] - xi[k-1]
  difference_coef <- sapply(2:m, function(k)(sample[,k] - sample[,k-1]))
  colnames(difference_coef) <- sapply(2:m, function(k)(paste0(k, "-", k-1)))
  
  ###  Verifying when the signal of the differences changes
  change_coef <- t(sapply(1:S,
                          function(l)(sapply(2:(m-1),
                                             function(k)(ifelse(sign(difference_coef[l,k] *
                                                                       difference_coef[l,k-1]) < 0, 1, 0))))))
  colnames(change_coef) <- sapply(strsplit(x = colnames(difference_coef), split = "-"), "[[", 2)[-1]
  
  ###  Saving results
  table_change <- rbind(apply(change_coef, 2, sum), S - apply(change_coef, 2, sum)) / S
  rownames(table_change) <- c("yes", "no")
  
  return(table_change)
  
}


# pp_change(sample = sample)

### Now you can go to file #4 to perform the proposed degree elevation methods,












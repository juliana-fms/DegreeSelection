###---------------------------------------------###
###   File #3: degree selection methods         ###
###   Bernstein Polynomials                     ###
###   Author: Juliana Freitas de Mello e Silva  ###
###   Email: julianafms27@gmail.com             ###
###---------------------------------------------###

###--------------------------------
###   Removing objects already allocated
###--------------------------------

rm(list = ls())

###--------------------------------
###   Packages required
###--------------------------------

library(BSDA)
# require(coda)

###--------------------------------
###   Loading necessary functions
###--------------------------------

source("utils\\functions.R")

###--------------------------------
###   Setting seed
###--------------------------------

set.seed(123456789)


###--------------------------------
###   Minimum and maximum degrees (min_m should be chosen so that we can
###   compare estimated results based on direct estimation of the BP using
###   degree m x BP with degree m bases on the degree elevation based in direct
###   estimation using the degree m-1.
###   So, for example, if you are interested in studying degrees greater or
###   equal to 10, min_m should be equal to 10 - 1
###--------------------------------

min_m <- 4 ; max_m <- 7
S <- 5000
level <- 0.10
sample_directory <- "data\\"
test <- "both"

optimal_m <- function(min_m, max_m, S, dir, level = 0.10, test = NULL){
  
  
  ###--------------------------------
  ###   Creating objects
  ###--------------------------------
  
  S <- 5000
  p_value_test1_sign <-
    p_value_test1_wilcox <-
    p_value_test2_sign <-
    p_value_test2_wilcox <- rep(x = NA, times = max_m - min_m - 1)
  names(p_value_test1_sign) <-
    names(p_value_test1_wilcox) <-
    names(p_value_test2_sign) <-
    names(p_value_test2_wilcox) <-  (min_m + 1):(max_m-1)
  
  opt_m_test1_sign <- opt_m_test1_wilcox <-
    opt_m_test2_sign <- opt_m_test2_wilcox <- NULL
  
  D_test1 <- D_test2 <- matrix(data = NA,
                               nrow = S,
                               ncol = max_m - min_m)
  colnames(D_test1) <- colnames(D_test2) <- (min_m + 1):(max_m)
  
  
  message("program in running ...")
  
  for(m in (min_m+1):max_m){
    
    
    ###--------------------------------
    ###   Reading results for m-1
    ###--------------------------------
    
    sample_coef_prev <- as.matrix(read.table(file = paste0(sample_directory,
                                                           "/posterior_m",
                                                           m-1,
                                                           ".txt")))
    
    ###--------------------------------
    ###   Reading results for m
    ###--------------------------------
    
    sample_coef <- as.matrix(read.table(file = paste0(sample_directory,
                                                            "/posterior_m",
                                                            m,
                                                            ".txt")))
    
    ###--------------------------------
    ###   Names of parameters
    ###--------------------------------
    
    names_coef_prev <- paste0("coef", 1:(m-1))
    colnames(sample_coef_prev) <- names_coef_prev
    
    names_coef_m <- paste0("coef", 1:m)
    colnames(sample_coef) <- names_coef_m
    
    ###--------------------------------
    ###   Coefficients via degree elevation
    ###--------------------------------
    
    sample_m_elevation <- degree_elevation_r1(xi = sample_coef_prev)
    
    ###--------------------------------###--------------------------------###
    
    ###--------------------------------
    ###   Test 1
    ###--------------------------------
    
    D_test1[,m-min_m] <- as.numeric(apply(abs(sample_coef - sample_m_elevation), 1, mean))
    
    ###--------------------------------
    ### `Test 2
    ###--------------------------------
    
    A <- sapply(1:m, function(l)(sapply(1:m, function(k)((1 / (2*m - 1)) *
                                                           dhyper(x = k-1,
                                                                  m = k+l-2,
                                                                  n = 2*m - k - l,
                                                                  k = m - 1)))))
    D_test2[,m-min_m] <- sapply(1:S,
                                function(a)(as.numeric((t((sample_coef -
                                                             sample_m_elevation)[a,]) %*% A) %*%
                                                         ((sample_coef -
                                                             sample_m_elevation)[a,]))))
    
    rm(sample_coef_prev,
       sample_coef,
       sample_m_elevation,
       A)
    
  }
  
  ###--------------------------------
  ###   P values - Sign test and Wilcoxon test
  ###--------------------------------
  
  for(m in (min_m+1):(max_m - 1)){
    
    p_value_test1_sign[m-min_m] <- SIGN.test(x = D_test1[,m-min_m], y = D_test1[,(m-min_m)+1], alternative = "greater")$p.value
    p_value_test1_wilcox[m-min_m] <- wilcox.test(x = D_test1[,m-min_m], y = D_test1[,(m-min_m)+1], alternative = "greater")$p.value
    p_value_test2_sign[m-min_m] <- SIGN.test(x = D_test2[,m-min_m], y = D_test2[,(m-min_m)+1], alternative = "greater")$p.value
    p_value_test2_wilcox[m-min_m] <- wilcox.test(x = D_test2[,m-min_m], y = D_test2[,(m-min_m)+1], alternative = "greater")$p.value
    
  }
  
  reference_p_value <- level
  
  ###--------------------------------
  ###   Optimal m according to criterion 1 based on both Sign and Wilcoxon tests
  ###--------------------------------
  
  opt_m_test1_sign <-  as.numeric(names(which(p_value_test1_sign > reference_p_value)[1])) + 1
  opt_m_test1_wilcox <-  as.numeric(names(which(p_value_test1_wilcox > reference_p_value)[1])) + 1
  
  
  ###--------------------------------
  ###   Optimal m according to criterion 2 based on both Sign and Wilcoxon tests
  ###--------------------------------
  
  opt_m_test2_sign <- as.numeric(names(which(p_value_test2_sign > reference_p_value)[1])) + 1
  opt_m_test2_wilcox <- as.numeric(names(which(p_value_test2_wilcox > reference_p_value)[1])) + 1
  
  if(test == "sign"){
    
    return(list(opt_m_criterion1 = opt_m_test1_sign,
                opt_m_criterion2 = opt_m_test2_sign))
    
  }else if(test == "wilcox"){
    
    return(list(opt_m_criterion1_wilcox = opt_m_test1_wilcox,
                opt_m_criterion2_wilcox = opt_m_test2_wilcox))
    
  }else if(test == "both") {
    return(list(opt_m_criterion1_sign = opt_m_test1_sign,
                opt_m_criterion1_wilcox = opt_m_test1_wilcox,
                opt_m_criterion2_sign = opt_m_test2_sign,
                opt_m_criterion2_wilcox = opt_m_test2_wilcox))
    
  }
  
}

optimal_m(min_m = 4, max_m = 12, S = 5000, dir = "data\\", level = 0.10, test = "both")


### Now you can go to file #5 to evaluate probabilities of change in the
### behavior of increasing/decreasing aspects of the approximation

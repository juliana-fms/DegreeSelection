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

sample_m1 <- as.matrix(read.table(file = "data\\posterior_m8.txt"))
sample_m2 <- as.matrix(read.table(file = "data\\posterior_m9.txt"))
sample_m3 <- as.matrix(read.table(file = "data\\posterior_m10.txt"))
level <- 0.10
test <- "both"
crit <- 1

optimal_m <- function(sample_m1, sample_m2, sample_m3,
                      level = 0.10, test = NULL, crit = NULL){
  
  
  ###--------------------------------
  ###   Creating objects
  ###--------------------------------
  
  S <- nrow(sample_m1)
  m1 <- ncol(sample_m2)
  
  p_value_test1_sign <-
    p_value_test1_wilcox <-
    p_value_test2_sign <-
    p_value_test2_wilcox <- NA
  
  opt_m_test1_sign <- opt_m_test1_wilcox <-
    opt_m_test2_sign <- opt_m_test2_wilcox <-
    opt_m_test1_sign_girls <- opt_m_test1_wilcox_girls <-
    opt_m_test2_sign_girls <- opt_m_test2_wilcox_girls <- NULL
  
  D_test1 <- D_test2 <- matrix(data = NA,
                               nrow = S,
                               ncol = 2)
  colnames(D_test1) <- colnames(D_test2) <- c(m1, m1+1)
  
  
  message("please, wait a few seconds...")
  
  psample_1 <- sample_m1
  psample_2 <- sample_m2
  
  for(m in 1:2){
    
    ###--------------------------------
    ###   Reading results for m-1
    ###--------------------------------
    
    sample_coef_prev <- psample_1
    
    ###--------------------------------
    ###   Reading results for m
    ###--------------------------------
    
    sample_coef <- psample_2
    
    ###--------------------------------
    ###   Names of parameters
    ###--------------------------------
    
    names_coef_prev <- paste0("coef", 1:(m1 + m - 2))
    colnames(sample_coef_prev) <- names_coef_prev
    
    names_coef_m <- paste0("coef", 1:(m1 + m - 1))
    colnames(sample_coef) <- names_coef_m
    
    ###--------------------------------
    ###   Coefficients via degree elevation
    ###--------------------------------
    
    sample_m_elevation <- degree_elevation_r1(xi = sample_coef_prev)
    colnames(sample_m_elevation) <- names_coef_m
    
    ###--------------------------------###--------------------------------###
    
    ###--------------------------------
    ###   Test 1
    ###--------------------------------
    
    D_test1[,m] <- as.numeric(apply(abs(sample_coef - sample_m_elevation), 1, mean))
    
    ###--------------------------------
    ### `Test 2
    ###--------------------------------
    
    m_aux <- m1 + m - 1
    A <- sapply(1:m_aux, function(l)(sapply(1:m_aux, function(k)((1 / (2*m_aux - 1)) *
                                                                   dhyper(x = k-1,
                                                                          m = k+l-2,
                                                                          n = 2*m_aux - k - l,
                                                                          k = m_aux - 1)))))
    D_test2[,m] <- sapply(1:S,
                          function(a)(as.numeric((t((sample_coef -
                                                       sample_m_elevation)[a,]) %*% A) %*%
                                                   ((sample_coef -
                                                       sample_m_elevation)[a,]))))
    
    rm(sample_coef_prev,
       sample_coef,
       sample_m_elevation,
       A)
    
    
    psample_1 <- sample_m2
    psample_2 <- sample_m3
    
  }
  
  ###--------------------------------
  ###   P values - Sign test and Wilcoxon test
  ###--------------------------------
  
  p_value_test1_sign <- SIGN.test(x = D_test1[,1], y = D_test1[,2], alternative = "greater")$p.value
  p_value_test1_wilcox <- wilcox.test(x = D_test1[,1], y = D_test1[,2], alternative = "greater")$p.value
  p_value_test2_sign <- SIGN.test(x = D_test2[,1], y = D_test2[,2], alternative = "greater")$p.value
  p_value_test2_wilcox <- wilcox.test(x = D_test2[,1], y = D_test2[,2], alternative = "greater")$p.value
  
  reference_p_value <- level
  
  ###--------------------------------
  ###   Optimal m according to criterion 1 based on both Sign and Wilcoxon tests
  ###--------------------------------
  
  opt_m_test1_sign <- ifelse(test = p_value_test1_sign < reference_p_value,
                             yes = "it is advised that you increase the degree",
                             no = paste0("degree ", m1 + 1, " is fine"))
  opt_m_test1_wilcox <-  ifelse(test = p_value_test1_wilcox < reference_p_value,
                                yes = "it is advised that you increase the degree",
                                no = paste0("degree ", m1 + 1, " is fine"))
  
  ###--------------------------------
  ###   Optimal m according to criterion 2 based on both Sign and Wilcoxon tests
  ###--------------------------------
  
  opt_m_test2_sign <- ifelse(test = p_value_test2_sign < reference_p_value,
                             yes = "it is advised that you increase the degree",
                             no = paste0("degree ", m1 + 1, " is fine"))
  opt_m_test2_wilcox <-  ifelse(test = p_value_test2_wilcox < reference_p_value,
                                yes = "it is advised that you increase the degree",
                                no = paste0("degree ", m1 + 1, " is fine"))
  
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

optimal_m(sample_m1, sample_m2, sample_m3, level = 0.10, test = "both")


### Now you can go to file #5 to evaluate probabilities of change in the
### behavior of increasing/decreasing aspects of the approximation

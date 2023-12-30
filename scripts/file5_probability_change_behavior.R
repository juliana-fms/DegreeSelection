###---------------------------------------------------###
###   File #6: application to Berkeley Growth Study   ###
###   Longitudinal data - Bernstein Polynomials       ###
###   Author: Juliana Freitas de Mello e Silva        ###
###---------------------------------------------------###


###--------------------------------
###   Removing objects already allocated
###--------------------------------

rm(list = ls())

###--------------------------------
###   Packages required
###--------------------------------

require(BSDA)

###--------------------------------
###   Defining work directory
###--------------------------------

directory <- paste0(getwd(), "\\")
data_directory <- paste0(directory, "data\\data_application\\")
results_directory <- paste0(directory, "results\\results_application\\")

###--------------------------------
###   Loading necessary functions
###--------------------------------

source("C:\\Users\\julia\\Desktop\\Projetos\\Doutorado\\DegreeSelection\\utils\\functions.R")

###--------------------------------
###   Setting seed
###--------------------------------

set.seed(123456789)

###--------------------------------###--------------------------------###

###--------------------------------
###   Minimum and maximum degrees
###--------------------------------

min_m <- 4 ; max_m <- 30
sample_size <- 5000
t_plot <- seq(from = 0, to = 1, length.out = 150)

###--------------------------------###--------------------------------###

###--------------------------------
###   Loading data
###--------------------------------

data_long  <- read.table(file = paste0(data_directory, "growth_long.txt"))
index_boys_long <- which(data_long$gender == 0)
index_girls_long <- which(data_long$gender == 1)
index_boys <- unique(data_long$id[which(data_long$gender == 0)])
index_girls <-unique(data_long$id[which(data_long$gender == 1)])
T_max <- ceiling(max(data_long$age))

###--------------------------------###--------------------------------###

###--------------------------------
###   Probabilities - boys
###--------------------------------

### Boys
m_boys <- 11 # optimal m
basis_plot_boys <- sapply(1:m_boys,
                          function(l)(dbeta(x = t_plot, shape1 = l,
                                            shape2 = m_boys - l + 1) / m_boys)) 

posterior_mu_xi_m7 <- as.matrix(read.table(file = paste0(results_directory, "m",
                                                         m_boys,
                                                         "/posterior_mu_xi_boys.txt")))
names_mu_xi_boys <-  paste0("mu_xi_boys", 1:m_boys)


### Difference between adjacent coefficients
difference_mu_xi_boys <- sapply(2:m_boys, function(a)(posterior_mu_xi_m7[,a] - posterior_mu_xi_m7[,a-1])) 
colnames(difference_mu_xi_boys) <- paste0(2:m_boys, "-", 1:(m_boys-1))

### Sign of the differences
sign_difference_mu_xi_boys <- sign(difference_mu_xi_boys)

### Change or no change
change_mu_xi_boys <- sapply(2:(m_boys-1),
                            function(a)(ifelse(test = (sign_difference_mu_xi_boys[,a-1] ==
                                                         sign_difference_mu_xi_boys[,a]),
                                               yes = 0, no = 1)))
colnames(change_mu_xi_boys) <- 2:(m_boys-1)

### Probabilities
probability_change_boys <- as.numeric(apply(change_mu_xi_boys, 2, sum) / sample_size) 

### Age points 
names(probability_change_boys)  <- (round(x = (((1:m_boys) - 1) / (m_boys - 1))*T_max, digits = 2))[2:(m_boys - 1)]


###--------------------------------
###   Probabilities - girls
###--------------------------------

### girls
m_girls <- 12 # optimal m
basis_plot_girls <- sapply(1:m_girls,
                           function(l)(dbeta(x = t_plot, shape1 = l,
                                             shape2 = m_girls - l + 1) / m_girls)) 

posterior_mu_xi_m7 <- as.matrix(read.table(file = paste0(results_directory, "m",
                                                         m_girls,
                                                         "/posterior_mu_xi_girls.txt")))
names_mu_xi_girls <-  paste0("mu_xi_girls", 1:m_girls)


### Difference between adjacent coefficients
difference_mu_xi_girls <- sapply(2:m_girls, function(a)(posterior_mu_xi_m7[,a] - posterior_mu_xi_m7[,a-1])) 
colnames(difference_mu_xi_girls) <- paste0(2:m_girls, "-", 1:(m_girls-1))

### Sign of the differences
sign_difference_mu_xi_girls <- sign(difference_mu_xi_girls)

### Change or no change
change_mu_xi_girls <- sapply(2:(m_girls-1),
                             function(a)(ifelse(test = (sign_difference_mu_xi_girls[,a-1] ==
                                                          sign_difference_mu_xi_girls[,a]),
                                                yes = 0, no = 1)))
colnames(change_mu_xi_girls) <- 2:(m_girls-1)

### Probabilities
probability_change_girls <- as.numeric(apply(change_mu_xi_girls, 2, sum) / sample_size) 

### Age points 
names(probability_change_girls)  <- (round(x = (((1:m_girls) - 1) / (m_girls - 1))*T_max, digits = 2))[2:(m_girls - 1)]



###--------------------------------###--------------------------------###

###--------------------------------
###   The probability of change inflection (decreasing / increasing)
###--------------------------------

probability_change_boys

probability_change_girls





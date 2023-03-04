###---------------------------------------------------###
###   File #5: application to Berkeley Growth Study   ###
###   Longitudinal data - Bernstein Polynomials       ###
###   Author: Juliana Freitas de Mello e Silva        ###
###---------------------------------------------------###

### Removing objetcs already alocated
rm(list = ls())

### Packages required
require(BSDA)

### Defining work directory
directory<- paste0(getwd(), "/") ; setwd(directory)
results_directory<- paste0(directory, "results/")

### Loading necessary functions
source("functions.R")

### Setting seed
set.seed(123456789)

###---------------------------

### Loading data
data_long<- read.table(file = paste0("accelerated_growth_long.txt"))
index_boys_long<- which(data_long$gender == 0)
index_girls_long<- which(data_long$gender == 1)
index_boys<- unique(data_long$id[which(data_long$gender == 0)])
index_girls<-unique(data_long$id[which(data_long$gender == 1)])

###
sample_size<- 5000
t_plot<- seq(from = 0, to = 1, length.out = 150)
max_m<- 30 ; T_max<- ceiling(max(data_long$age))

###-----------------------###
### Probabilities - boys  ###
###-----------------------###

### Boys
m_boys<- 7 # optimal m
basis_plot_boys<- sapply(1:m_boys, function(l)(dbeta(x = t_plot, shape1 = l, shape2 = m_boys - l + 1) / m_boys)) # Bernstein basis

posterior_sample_m7<- as.matrix(read.table(file = paste0(results_directory, "m", m_boys, "/posterior_sample.txt"))) # reading posterior sample
names_mu_xi_boys<- paste0("mu_ksi_boys.", 1:m_boys, ".") # names of the parameters
xi_boys<- posterior_sample_m7[,names_mu_xi_boys] # coefficients

### Difference between adjacent coefficients
difference_xi_boys<- sapply(2:m_boys, function(a)(xi_boys[,a] - xi_boys[,a-1])) 
colnames(difference_xi_boys)<- paste0(2:m_boys, "-", 1:(m_boys-1))

### Signal  of the differences
signal_difference_xi_boys<- sign(difference_xi_boys)

### Change or no change
change_xi_boys<- sapply(2:(m_boys-1), function(a)(ifelse(test = (signal_difference_xi_boys[,a-1] == signal_difference_xi_boys[,a]), yes = 0, 1)))
colnames(change_xi_boys)<- 2:(m_boys-1)

### Probabilities
probability_change_boys<- as.numeric(apply(change_xi_boys, 2, sum) / sample_size) 

### Age points 
(((1:m_boys) - 1) / (m_boys - 1))*T_max

###-----------------------###
### Probabilities - girls  ###
###-----------------------###

### girls
m_girls<- 8 # optimal m
basis_plot_girls<- sapply(1:m_girls, function(l)(dbeta(x = t_plot, shape1 = l, shape2 = m_girls - l + 1) / m_girls)) # Bernstein basis

posterior_sample_m8<- as.matrix(read.table(file = paste0(results_directory, "m", m_girls, "/posterior_sample.txt"))) # reading posterior sample
names_mu_xi_girls<- paste0("mu_ksi_girls.", 1:m_girls, ".") # names of the parameters
xi_girls<- posterior_sample_m8[,names_mu_xi_girls] # coefficients

### Difference between adjacent coefficients
difference_xi_girls<- sapply(2:m_girls, function(a)(xi_girls[,a] - xi_girls[,a-1]))
colnames(difference_xi_girls)<- paste0(2:m_girls, "-", 1:(m_girls-1))

### Signal  of the differences
signal_difference_xi_girls<- sign(difference_xi_girls)

### Change or no change
change_xi_girls<- sapply(2:(m_girls-1), function(a)(ifelse(test = (signal_difference_xi_girls[,a-1] == signal_difference_xi_girls[,a]), yes = 0, 1)))
colnames(change_xi_girls)<- 2:(m_girls-1)

### Probabilities
probability_change_girls<- as.numeric(apply(change_xi_girls, 2, sum) / sample_size) 

### Age points
round(x = (((1:m_girls) - 1) / (m_girls - 1))*T_max, digits = 2)



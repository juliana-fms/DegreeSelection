###---------------------------------------------------###
###   File #4b: application to Berkeley Growth Study  ###
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
###   Loading data
###--------------------------------

data_long <- read.table(file = paste0(data_directory, "growth_long.txt"))
index_boys_long <- which(data_long$gender == 0)
index_girls_long <- which(data_long$gender == 1)

###--------------------------------
###   Necessary quantities
###--------------------------------

sample_size <- 5000
t_plot <- seq(from = 0, to = 1, length.out = 150)
min_m <- 24 ; max_m <- 30 ; T_max <- ceiling(max(data_long$age))

###--------------------------------
###   Creating objects
###--------------------------------

dic_marg <- lpml_marg <- waic_marg <- rep(x = NA, times = (max_m-min_m+1))
names(dic_marg) <- names(lpml_marg) <- names(waic_marg) <- min_m:max_m

###--------------------------------
###   Reading results
###--------------------------------

for(m in min_m:max_m){
  
  ###--------------------------------
  ###   Bernstein basis
  ###--------------------------------
  
  basis <- BP_Basis(t = data_long$age / T_max, m = m)
  
  ###--------------------------------
  ###   Variables names
  ###--------------------------------
  
  names_mu_ksi_boys <- paste0("mu_ksi_boys.", 1:m, ".")
  names_mu_ksi_girls <- paste0("mu_ksi_girls.", 1:m, ".")
  names_sigma2_e <- "sigma2_e"
  
  ###--------------------------------
  ###   Reading comparison measures
  ###--------------------------------
  
  comparison_measures_marg <- read.table(file = paste0(results_directory, "m", m, "/comparison_measures_marg.txt"), row.names = NULL)
  dic_marg[m-min_m + 1] <- as.numeric(comparison_measures_marg[1,2])
  lpml_marg[m-min_m + 1] <- as.numeric(comparison_measures_marg[2,2])
  waic_marg[m-min_m + 1] <- as.numeric(comparison_measures_marg[3,2])
  
}

###--------------------------------
###   Sorting comparison measures (varying m)
###--------------------------------

sort(dic_marg, decreasing = FALSE)
sort(lpml_marg, decreasing = FALSE)
sort(waic_marg, decreasing = FALSE)


### Now you can go to file #5 to evaluate probabilities of change in the
### behavior of increasing / decreasing aspects of the approximation

















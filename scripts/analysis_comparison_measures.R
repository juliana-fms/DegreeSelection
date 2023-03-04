###---------------------------------------------------###
###   File #4: application to Berkeley Growth Study   ###
###   Longitudinal data - Bernstein Polynomials       ###
###   Author: Juliana Freitas de Mello e Silva        ###
###---------------------------------------------------###

### Removing objects already allocated
rm(list = ls())

### Packages required
require(BSDA)

### Defining work directory
directory<- paste0(getwd(), "/") ; setwd(directory)
data_directory<- paste0(getwd(), "/")
results_directory<- "results/"

### Loading necessary functions
source("functions.R")

### Setting seed
set.seed(123456789)

###---------------------------

### Loading data
data_long<- read.table(file = paste0(data_directory, "accelerated_growth_long.txt"))
index_boys_long<- which(data_long$gender == 0)
index_girls_long<- which(data_long$gender == 1)

### Necessary quantities
sample_size<- 5000
t_plot<- seq(from = 0, to = 1, length.out = 150)
min_m<- 5 ; max_m<- 30 ; T_max<- ceiling(max(data_long$age))

### Creting objects
dic_marg<- rep(x = NA, times = (max_m-4+1)) ; names(dic_marg)<- 4:max_m ; lpml_marg<- waic_marg<- dic_marg


### Reading results
for(m in min_m:max_m){
  
  ### Bernstein basis
  basis<- BP_Basis(t = data_long$age / T_max, m = m)
  
  ### Variables names
  names_mu_ksi_boys<- paste0("mu_ksi_boys.", 1:m, ".")
  names_mu_ksi_girls<- paste0("mu_ksi_girls.", 1:m, ".")
  names_sigma2_e<- "sigma2_e"
  
  ### Reading comparison measures
  comparison_measures_marg<- read.table(file = paste0(results_directory, "m", m, "/comparison_measures_marg.txt"), row.names = NULL)
  dic_marg[m-3]<- as.numeric(comparison_measures_marg[1,2])
  lpml_marg[m-3]<- as.numeric(comparison_measures_marg[2,2])
  waic_marg[m-3]<- as.numeric(comparison_measures_marg[3,2])
  
}

### Sorting comparison measures (varying m)
sort(dic_marg, decreasing = FALSE)
sort(lpml_marg, decreasing = FALSE)
sort(waic_marg, decreasing = FALSE)



















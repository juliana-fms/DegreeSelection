###---------------------------------------------------###
###   File #4: application to Berkeley Growth Study   ###
###   Longitudinal data - Bernstein Polynomials       ###
###   Author: Juliana Freitas de Mello e Silva        ###
###---------------------------------------------------###

### Removing objects already allocated
rm(list = ls())

### Packages required
require(BSDA)
require(coda)

### Defining work directory
directory<- paste0(getwd(), "/") ; setwd(directory)
data_directory<- getwd()
results_directory<- "results/"

### Loading necessary functions
source("functions.R")

### Setting seed
set.seed(123456789)


###-----------------------------
min_m<- 5 ; max_m<- 30 # minimum and maximum degrees
###-----------------------------

###-----------------------------
### Creating objects
sample_size<- 5000
p_value_test1_sign_boys<- rep(x = NA, times = max_m - min_m) ; opt_m_test1_sign_boys<- NULL
p_value_test1_wilcox_boys<- rep(x = NA, times = max_m - min_m) ; opt_m_test1_wilcox_boys<- NULL
p_value_test1_sign_girls<- rep(x = NA, times = max_m - min_m) ; opt_m_test1_sign_girls<- NULL
p_value_test1_wilcox_girls<- rep(x = NA, times = max_m - min_m) ; opt_m_test1_wilcox_girls<- NULL
p_value_test2_sign_boys<- rep(x = NA, times = max_m - min_m) ; opt_m_test2_sign_boys<- NULL
p_value_test2_wilcox_boys<- rep(x = NA, times = max_m - min_m) ; opt_m_test2_wilcox_boys<- NULL
p_value_test2_sign_girls<- rep(x = NA, times = max_m - min_m) ; opt_m_test2_sign_girls<- NULL
p_value_test2_wilcox_girls<- rep(x = NA, times = max_m - min_m) ; opt_m_test2_wilcox_girls<- NULL

D_test1_boys<- matrix(data = NA, nrow = sample_size, ncol = max_m - min_m + 1) ; colnames(D_test1_boys)<- min_m:max_m ; D_test1_girls<- D_test1_boys
D_test2_boys<- matrix(data = NA, nrow = sample_size, ncol = max_m - min_m + 1) ; colnames(D_test2_boys)<- min_m:max_m ; D_test2_girls<- D_test2_boys

for(m in 5:max_m){
  
  ### Reading results for m-1
  posterior_sample_prev<- as.matrix(read.table(file = paste0(results_directory, "m", m-1, "/posterior_sample.txt")))
  ### Reading results for m
  posterior_sample_m<- as.matrix(read.table(file = paste0(results_directory, "m", m, "/posterior_sample.txt")))
  
  ### Names of parameters
  names_mu_xi_boys_prev<- paste0("mu_ksi_boys.", 1:(m-1), ".") ; names_mu_xi_girls_prev<- paste0("mu_ksi_girls.", 1:(m-1), ".")
  names_mu_xi_boys_m<- paste0("mu_ksi_boys.", 1:m, ".") ; names_mu_xi_girls_m<- paste0("mu_ksi_girls.", 1:m, ".")
  
  ### Mean of coefficients for boys and girls, given that degree is m-1
  xi_prev_boys<- posterior_sample_prev[,names_mu_xi_boys_prev] ; xi_prev_girls<- posterior_sample_prev[,names_mu_xi_girls_prev]
  ### Mean of coefficients for boys and girls, via degree elevation
  xi_elevation_boys<- degree_elevation_r1(xi = xi_prev_boys) ; xi_elevation_girls<- degree_elevation_r1(xi = xi_prev_girls)
  ### Mean of coefficients for boys and girls, given that degree is m
  xi_m_boys<- posterior_sample_m[,names_mu_xi_boys_m] ; xi_m_girls<- posterior_sample_m[,names_mu_xi_girls_m]
  
  # Test 1
  D_test1_boys[,m-4]<- as.numeric(apply(abs(xi_m_boys - xi_elevation_boys), 1, mean))
  D_test1_girls[,m-4]<- as.numeric(apply(abs(xi_m_girls - xi_elevation_girls), 1, mean))
  
  # Test 2
  A<- sapply(1:m, function(l)(sapply(1:m, function(k)((1 / (2*m - 1))*dhyper(x = k-1, m = k+l-2, n = 2*m - k - l, k = m - 1)))))
  D_test2_boys[,m-4]<- sapply(1:sample_size, function(a)(as.numeric((t((xi_m_boys - xi_elevation_boys)[a,])%*%A)%*%((xi_m_boys - xi_elevation_boys)[a,]))))
  D_test2_girls[,m-4]<- sapply(1:sample_size, function(a)(as.numeric((t((xi_m_girls - xi_elevation_girls)[a,])%*%A)%*%((xi_m_girls - xi_elevation_girls)[a,]))))
  
  rm(xi_prev_boys, xi_prev_girls, xi_elevation_boys, xi_elevation_girls, xi_m_boys, xi_m_girls, A)
  #print(m)
}

### P values - sign test and wilcoxon test
for(m in min_m:(max_m-1)){
  p_value_test1_sign_boys[m-4]<- SIGN.test(x = D_test1_boys[,m-4], y = D_test1_boys[,(m-4)+1], alternative = "greater")$p.value
  p_value_test1_wilcox_boys[m-4]<- wilcox.test(x = D_test1_boys[,m-4], y = D_test1_boys[,(m-4)+1], alternative = "greater")$p.value
  p_value_test1_sign_girls[m-4]<- SIGN.test(x = D_test1_girls[,m-4], y = D_test1_girls[,(m-4)+1], alternative = "greater")$p.value
  p_value_test1_wilcox_girls[m-4]<- wilcox.test(x = D_test1_girls[,m-4], y = D_test1_girls[,(m-4)+1], alternative = "greater")$p.value
  
  p_value_test2_sign_boys[m-4]<- SIGN.test(x = D_test2_boys[,m-4], y = D_test2_boys[,(m-4)+1], alternative = "greater")$p.value
  p_value_test2_wilcox_boys[m-4]<- wilcox.test(x = D_test2_boys[,m-4], y = D_test2_boys[,(m-4)+1], alternative = "greater")$p.value
  p_value_test2_sign_girls[m-4]<- SIGN.test(x = D_test2_girls[,m-4], y = D_test2_girls[,(m-4)+1], alternative = "greater")$p.value
  p_value_test2_wilcox_girls[m-4]<- wilcox.test(x = D_test2_girls[,m-4], y = D_test2_girls[,(m-4)+1], alternative = "greater")$p.value
}


### Optimal m according to criterion 1 based on both Sign and Wilcoxon tests, for boys
opt_m_test1_sign_boys<- which(p_value_test1_sign_boys > 0.10)[1] + 5
opt_m_test1_wilcox_boys<- which(p_value_test1_wilcox_boys > 0.10)[1] + 5

### Optimal m according to criterion 1 based on both Sign and Wilcoxon tests, for girls
opt_m_test1_sign_girls<- which(p_value_test1_sign_girls > 0.10)[1] + 5
opt_m_test1_wilcox_girls<- which(p_value_test1_wilcox_girls > 0.10)[1] + 5

### Optimal m according to criterion 2 based on both Sign and Wilcoxon tests, for boys
opt_m_test2_sign_boys<- which(p_value_test2_sign_boys > 0.10)[1] + 5
opt_m_test2_wilcox_boys<- which(p_value_test2_wilcox_boys > 0.10)[1] + 5

### Optimal m according to criterion 2 based on both Sign and Wilcoxon tests, for girls
opt_m_test2_sign_girls<- which(p_value_test2_sign_girls > 0.10)[1] + 5
opt_m_test2_wilcox_girls<- which(p_value_test2_wilcox_girls > 0.10)[1] + 5

### Showing otpimal degrees
c(opt_m_test1_sign_boys, opt_m_test1_wilcox_boys, opt_m_test2_sign_boys, opt_m_test2_wilcox_boys)
c(opt_m_test1_sign_girls, opt_m_test1_wilcox_girls, opt_m_test2_sign_girls, opt_m_test2_wilcox_girls)


###---------------------------------------###
### Results based on the optimal degrees  ###
###---------------------------------------###

### Reading data
data_long<- read.table(file = "accelerated_growth_long.txt")
index_boys<- unique(data_long$id[which(data_long$gender == 0)])
index_girls<- unique(data_long$id[which(data_long$gender == 1)])
T_max<- ceiling(max(data_long$age))

### Times for plot
t_plot<- seq(from = min(data_long$age), to = T_max, length.out = 100)

### Colors
colors_plot<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###-----------------------------
### Figure for girls
m<- 8
posterior_sample<- as.matrix(read.table(file = paste0(results_directory, "m", m, "/posterior_sample.txt")))
posterior_sample_mu_girls_m8<- posterior_sample[, paste0("mu_ksi_girls.", 1:m, ".")] ; rm(posterior_sample) ; rm(m)
basis_plot_girls_m8<- BP_Basis(t = t_plot / T_max, m = 8)
mu_plot_girls_m8<- basis_plot_girls_m8%*%t(posterior_sample_mu_girls_m8)


pdf(file = paste0(results_directory, "estimated_mean_girls_m8_degree_elevation.pdf"), width = 8, height = 6)
par(mar = c(5, 5, 3, 1) + 0.1)
m_girls<- 8
### Plotting observed trajectories
plot(x = data_long$age[which(data_long$id == index_girls[1])], y = data_long$height[which(data_long$id == index_girls[1])],
     type = "l",
     xlim = c(min(data_long$age), max(data_long$age)), ylim = c(0, max(data_long$height)),
     xlab = "Age (in years)", ylab = "Acceleration of growth in girls (in cm)", cex.lab = 1.2,
     col = colors_plot[1], lty = 1, lwd = 1,
     xaxt = "n")
axis(side = 1, at = unique(data_long$age))
unlist(sapply(index_girls[-1], function(i)(lines(x = data_long$age[which(data_long$id == i)], y = data_long$height[which(data_long$id == i)],
                                                 col = colors_plot[1], lty = 1))))
axis(side = 3, at = round(x = (((1:m_girls) - 1)/(m_girls - 1))*T_max, digits = 2), col.ticks = colors_plot[1], col.axis = colors_plot[1],
     cex.axis = 0.84)

### Observed mean - girls
lines(x = sort(unique(data_long$age)),
      y = sapply(1:length(unique(data_long$age)),
                 function(j)(mean(data_long$height[which(data_long$age == sort(unique(data_long$age))[j] & data_long$gender == 1)]))),
      col = 1, lty = 3, lwd = 2)

### Estimated median and HPD
lines(x = t_plot, y = apply(mu_plot_girls_m8, 1, median), type = "l", col = 1, lty = 1, lwd = 2)
lines(x = t_plot, y = HPDinterval(obj = as.mcmc(t(mu_plot_girls_m8)), prob = 0.95)[,"lower"], type = "l", col = 1, lty = 2, lwd = 2)
lines(x = t_plot, y = HPDinterval(obj = as.mcmc(t(mu_plot_girls_m8)), prob = 0.95)[,"upper"], type = "l", col = 1, lty = 2, lwd = 2)

legend(x = "topright", legend = c("obs. trajectory", "obs. mean girls", "est. mean", "est. HPD"), col = c(colors_plot[1], rep(x = 1, times = 4)),
       lty = c(1, 3, 1, 2), lwd = c(1, rep(x = 2, times = 3)), bty = "n")
dev.off()
###-----------------------------


###-----------------------------
### Figure for boys
m<- 7 # opimal degree
posterior_sample<- as.matrix(read.table(file = paste0(results_directory, "m", m, "/posterior_sample.txt"))) # reading results 
posterior_sample_mu_boys_m7<- posterior_sample[, paste0("mu_ksi_boys.", 1:m, ".")] ; rm(posterior_sample) ; rm(m) 
basis_plot_boys_m7<- BP_Basis(t = t_plot / T_max, m = 7)
mu_plot_boys_m7<- basis_plot_boys_m7%*%t(posterior_sample_mu_boys_m7)


pdf(file = paste0(results_directory, "estimated_mean_boys_m7_degree_elevation.pdf"), width = 8, height = 6)
par(mar = c(5, 5, 3, 1) + 0.1)
m_boys<- 7
### Plotting observed trajectories
plot(x = data_long$age[which(data_long$id == index_boys[1])], y = data_long$height[which(data_long$id == index_boys[1])],
     type = "l",
     xlim = c(min(data_long$age), max(data_long$age)), ylim = c(0, max(data_long$height)),
     xlab = "Age (in years)", ylab = "Acceleration of growth in boys (in cm)", cex.lab = 1.2,
     col = colors_plot[1], lty = 1, lwd = 1,
     xaxt = "n")
axis(side = 1, at = unique(data_long$age))
unlist(sapply(index_boys[-1], function(i)(lines(x = data_long$age[which(data_long$id == i)], y = data_long$height[which(data_long$id == i)],
                                                 col = colors_plot[1], lty = 1))))

axis(side = 3, at = round(x = (((1:m_boys) - 1) / (m_boys - 1))*T_max, digits = 2),  col.ticks = colors_plot[1], col.axis = colors_plot[1],
     cex.axis = 0.84)

### Observed mean - boys
lines(x = sort(unique(data_long$age)),
      y = sapply(1:length(unique(data_long$age)),
                 function(j)(mean(data_long$height[which(data_long$age == sort(unique(data_long$age))[j] & data_long$gender == 0)]))),
      col = 1, lty = 3, lwd = 2)

### Estimated median and HPD
lines(x = t_plot, y = apply(mu_plot_boys_m7, 1, median), type = "l", col = 1, lty = 1, lwd = 2)
lines(x = t_plot, y = HPDinterval(obj = as.mcmc(t(mu_plot_boys_m7)), prob = 0.95)[,"lower"], type = "l", col = 1, lty = 2, lwd = 2)
lines(x = t_plot, y = HPDinterval(obj = as.mcmc(t(mu_plot_boys_m7)), prob = 0.95)[,"upper"], type = "l", col = 1, lty = 2, lwd = 2)

legend(x = "topright", legend = c("obs. trajectory", "obs. mean boys", "est. mean", "est. HPD"), col = c(colors_plot[1], rep(x = 1, times = 4)),
       lty = c(1, 3, 1, 2), lwd = c(1, rep(x = 2, times = 3)), bty = "n")
dev.off()
###-----------------------------




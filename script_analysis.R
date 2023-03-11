###---------------------------------------------------###
###   Script to run analysis on Berkeley Growth data  ###
###   Author: Juliana Freitas de Mello e Silva        ###
###---------------------------------------------------###


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

# require(data.table)
library(dplyr)
library(ggplot2)
# require(mcmcplots)
# require(MCMCvis)
# require(mvtnorm)
# require(R2jags)

###-------------------------------
###   Working directories
###-------------------------------

directory <- paste0(getwd(), "\\")
data_directory <- paste0(directory, "data\\")
figures_directory <- paste0(directory, "figures\\")
tables_directory <- paste0(directory, "tables\\")

###-------------------------------###-------------------------------###

### Graphical definitions
colour_plot <- c("#999999", "#F0E442", "#E69F00", "#D55E00", "#56B4E9", "#0072B2", "#009E73", "#CC79A7")
size_plot <- 12
vjust_plot <- 1.8
theme_plot <- theme_bw() +
  theme(axis.text.x = element_text(size = size_plot),
        axis.text.y = element_text(size = size_plot),
        axis.title.x = element_text(size = size_plot, vjust = -vjust_plot), 
        axis.title.y = element_text(size = size_plot, vjust = -vjust_plot), 
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = size_plot),
        legend.title = element_text(size = size_plot),
        legend.position = "bottom",
        plot.margin = unit(x = c(1, 1, 1, 0.5), units = "cm"),
        plot.title = element_text(hjust = 0.5)
  )

###------------------------------###------------------------------###

###-------------------------------
###   Important settings for this analysis
###-------------------------------



###------------------------------###------------------------------###

###-------------------------------
###   Calling functions
###-------------------------------

# source("functions.R")

###------------------------------###------------------------------###

###-------------------------------
###   Reading data
###-------------------------------

data_long <- read.table(file = "data_growth_long.txt")


###------------------------------###------------------------------###

###------------------------------
###   MCMC specifications
###------------------------------

### Parameters
parameters <- c("sigma2_e", "mu_ksi_boys", "mu_ksi_girls", "Sigma_ksi", "log_lik")

### Chain specifications
burn_in <- 80000
lag <- 30
sample_size <- 2500
number_iterations <- burn_in + (sample_size)*lag

### JAGS model
model <- file.path("model_growth_study.R")
# file.show(model)

###------------------------------

### Data for JAGS
data_jags <- list() # create object
data_jags$n <- length(unique(data_long$id))
data_jags$J_i <- as.numeric(table(data_long$id)) # number of measurements
data_jags$W <- data_long$height
data_jags$N <- sum(data_jags$J_i)
data_jags$id <- data_long$id
T_max <- max(data_long$age)
data_jags$index_boys <- unique(data_long$id[which(data_long$gender == 0)])
data_jags$index_girls <- unique(data_long$id[which(data_long$gender == 1)])
index_boys_long <- which(data_long$gender == 0)
index_girls_long <- which(data_long$gender == 1)

for(m in 5:30){
   
   data_jags$m <- m
   
   results_directory <- paste0(directory, "/results/m", data_jags$m, "/")
   dir.create(path = results_directory, showWarnings = TRUE)
   
   ###------------------------------
   ###  Bernstein Basis
   ###------------------------------
   
   data_jags$basis<- BP_Basis(t = data_long$age / T_max, m = data_jags$m)
   
   ###------------------------------
   ###  Prior specifications
   ###------------------------------
   
   data_jags$mu_mu_ksi <- rep(x = 0, times = data_jags$m)
   data_jags$Tau_mu_ksi <- diag(x = 1/100, nrow = data_jags$m)
   data_jags$I <- diag(x = 1, nrow = data_jags$m)
   
   ###------------------------------
   ###  Initial values
   ###------------------------------
   
   initial_values <- function(){
     list(tau_e = rexp(n = 1, rate = 1),
          ksi = rmvnorm(n = data_jags$n,
                        mean = rep(x = 0, times = data_jags$m),
                        sigma = diag(x = 1,nrow = data_jags$m)),
          mu_ksi_boys = rnorm(n = data_jags$m, mean = 0, sd = 1),
          mu_ksi_girls = rnorm(n = data_jags$m, mean = 0, sd = 1),
          Tau_ksi = rWishart(n = 1, df = data_jags$m, Sigma = data_jags$m*diag(nrow = data_jags$m))[,,])
   }
   
   ###------------------------------
   ###  Compiling on JAGS
   ###------------------------------
   
   output <- jags(data = data_jags, 
                 inits = initial_values,
                 parameters = parameters,
                 model.file = model,
                 n.chains = 2,
                 n.iter = number_iterations,
                 n.burnin = burn_in,
                 n.thin = lag,
                 working.directory = directory,
                 DIC = FALSE)
   
   ###-----------------------------
   ###  Checking convergence
   ###------------------------------
   
   # mcmcplot(mcmcout = output, parms = parameters[-length(parameters)], dir = results_directory,
   # extension = "html",
   #  greek = TRUE)
   
   ###------------------------------
   ###  Variable names
   ###------------------------------
   
   names_mu_xi_boys <- paste0("mu_ksi_boys[", 1:data_jags$m, "]")
   names_mu_xi_girls <- paste0("mu_ksi_girls[", 1:data_jags$m, "]")
   names_sigma2_e <- "sigma2_e"
   names_Sigma_xi <- paste0("Sigma_ksi[", expand.grid(1:data_jags$m, 1:data_jags$m)[,1], ",",
                            expand.grid(1:data_jags$m, 1:data_jags$m)[,2], "]")

   ###------------------------------
   ###  Saving posterior sample
   ###------------------------------
   
   posterior_sample <- as.matrix(as.mcmc(x = output, start = burn_in, end = number_iterations,
                                         thin = lag))
   posterior_mu_xi_boys <- posterior_sample[,names_mu_xi_boys]
   posterior_mu_xi_girls <- posterior_sample[,names_mu_xi_girls]
   posterior_Sigma_xi <- posterior_sample[,names_Sigma_xi]
   posterior_sigma2_e <- posterior_sample[,names_sigma2_e]
   write.table(x = as.matrix(posterior_sample), file = paste0(results_directory, "posterior_sample.txt"))
   
   ###-----------------------------
   ###  Calculating differences between xi[k] - xi[k-1]
   ###------------------------------
   
   difference_xi_boys <- sapply(2:data_jags$m,
                                function(k)(posterior_mu_xi_boys[,k] - posterior_mu_xi_boys[,k-1]))
   colnames(difference_xi_boys) <- sapply(2:data_jags$m,
                                          function(k)(paste0(k, "-", k-1)))
   difference_xi_girls <- sapply(2:data_jags$m,
                                 function(k)(posterior_mu_xi_girls[,k] - posterior_mu_xi_girls[,k-1]))
   colnames(difference_xi_girls) <- sapply(2:data_jags$m,
                                           function(k)(paste0(k, "-", k-1)))
   
   
   ###-----------------------------
   ###  Verifying when the signal of the differences changes
   ###------------------------------
   
   change_xi_boys <- t(sapply(1:nrow(posterior_mu_xi_boys),
                              function(l)(sapply(2:(data_jags$m-1),
                                                 function(k)(ifelse(sign(difference_xi_boys[l,k] *
                                                                           difference_xi_boys[l,k-1]) < 0, 1, 0))))))
   colnames(change_xi_boys) <- sapply(strsplit(x = colnames(difference_xi_boys), split = "-"),
                                      "[[", 2)[-1]
   
   change_xi_girls <- t(sapply(1:nrow(posterior_mu_xi_girls),
                               function(l)(sapply(2:(data_jags$m-1),
                                                  function(k)(ifelse(sign(difference_xi_girls[l,k] *
                                                                            difference_xi_girls[l,k-1]) < 0, 1, 0))))))
   colnames(change_xi_girls) <- sapply(strsplit(x = colnames(difference_xi_girls), split = "-"),
                                       "[[", 2)[-1]
   
   ###-----------------------------
   ###  Saving results
   ###-----------------------------
   
   change_xi_boys <- t(sapply(1:nrow(posterior_mu_xi_boys),
                              function(l)(sapply(2:(data_jags$m-1),
                                                 function(k)(ifelse(sign(difference_xi_boys[l,k] *
                                                                           difference_xi_boys[l,k-1]) < 0, 1, 0))))))
   colnames(change_xi_boys) <- sapply(strsplit(x = colnames(difference_xi_boys), split = "-"),
                                      "[[", 2)[-1]
   
   change_xi_girls <- t(sapply(1:nrow(posterior_mu_xi_girls),
                               function(l)(sapply(2:(data_jags$m-1),
                                                  function(k)(ifelse(sign(difference_xi_girls[l,k] *
                                                                            difference_xi_girls[l,k-1]) < 0, 1, 0))))))
   colnames(change_xi_girls) <- sapply(strsplit(x = colnames(difference_xi_girls), split = "-"),
                                       "[[", 2)[-1]
   
   table_change_boys <- rbind(apply(change_xi_boys, 2, sum), nrow(change_xi_boys) - apply(change_xi_boys,
                                                                                          2, sum))
   table_change_girls <- rbind(apply(change_xi_girls, 2, sum), nrow(change_xi_girls) -
                                 apply(change_xi_girls, 2, sum))
   
   write.table(x = table_change_boys, file =  paste0(results_directory,
                                                     "table_change_point_mean_curve_boys.txt"),
               row.names = c("change", "no change"))
   write.table(x = table_change_girls, file =  paste0(results_directory,
                                                      "table_change_point_mean_curve_girls.txt"),
               row.names = c("change", "no change"))

   ###-----------------------------
   ###  Log-likelihood and comparison measures
   ###-----------------------------
   
   log_likelihood_marg <- matrix(data = NA, nrow = data_jags$N, ncol = nrow(posterior_sample))
   aux_marginal_variance <- matrix(data = 0, nrow = data_jags$N, ncol = nrow(posterior_sample))
   
   for(ss in 1:nrow(posterior_sample)){
      for(ii in 1:data_jags$m){
         for(jj in 1:data_jags$m){
            aux_marginal_variance[,ss] <- aux_marginal_variance[,ss] + 
               (data_jags$basis[,ii] * data_jags$basis[,jj] *
                  posterior_sample[ss,paste0("Sigma_ksi[",ii, ",", jj, "]")])  
         }
      }
   }
   
   for(ss in 1:nrow(posterior_sample)){
      for(ll in index_boys_long){
         log_likelihood_marg[ll,ss]<- dnorm(x = data_jags$W[ll],
                                            mean = data_jags$basis[ll,] %*%
                                              posterior_sample[ss,names_mu_xi_boys],
                                            sd = sqrt(posterior_sigma2_e[ss] + aux_marginal_variance[ll,ss] ), log = TRUE)
      }
      rm(ll)
      for(ll in index_girls_long){
         log_likelihood_marg[ll,ss]<- dnorm(x = data_jags$W[ll],
                                            mean = data_jags$basis[ll,]%*%
                                              posterior_sample[ss,names_mu_xi_girls], 
                                            sd = sqrt(posterior_sigma2_e[ss] + aux_marginal_variance[ll,ss] ), log = TRUE)
      }
      rm(ll)
   }
   
   ### Calculating comparison measures
   lpml_marg <- CPO(loglik = t(log_likelihood_marg), n = sample_size)[1]
   waic_marg <- WAIC(loglik = t(log_likelihood_marg), n = sample_size)[1]
   dic_marg <- DIC(loglik = t(log_likelihood_marg), n = sample_size)[1]
   
   comparison_measures_marg <- c(dic_marg, -2*lpml_marg, -2*waic_marg)
   names(comparison_measures_marg)<- c("dic", "lpml", "waic")
   
   write.table(x = comparison_measures_marg, file = paste0(results_directory, "comparison_measures_marg.txt"))

}














###--------------------------------
###   Simulation study
###--------------------------------

###--------------------------------
###   Removing already allocated objects
###--------------------------------

rm(list = ls())

###--------------------------------
###   Required packages
###--------------------------------

library(mcmcplots)
library(MCMCvis)
library(mvtnorm)
library(R2jags)

###--------------------------------
###   Defining work directory
###--------------------------------

directory <- paste0(getwd(), "/")
model_directory<- paste0(directory, "scripts\\models\\")

###--------------------------------
###   Loading necessary functions
###--------------------------------

source(paste0(directory, "functions_DS.R"))

###--------------------------------
###   Setting seed
###--------------------------------

set.seed(123456789)

###---------------------------###--------------------------------###

###--------------------------------
###   Reading data
###--------------------------------

a <- 1 ; b <- 3
data_directory<- paste0("C:\\Users\\julia\\Desktop\\Projetos\\Doutorado\\tese\\codigos\\sanduiche\\3_09October2018\\data\\beta", a, b, "\\")

for(r in 1:50){
  
  data<- read.table(file = paste0(data_directory, "data", r, ".txt"))
  dir.create(path = paste0(directory, "beta", a, b, "/data", r), showWarnings = FALSE)
  
  # Values for m
  # m<- ceiling(seq(from = 4, to = max(as.numeric(table(data$id))), by = 1))
  m<- unique(ceiling(nrow(data)^(seq(from = 0.3, to = 0.75, by = 0.05))))
  
  ###---------------------------
  
  ### Data for JAGS
  data_jags<- list() # create object
  data_jags$n<- length(unique(data$id))
  data_jags$J_i<- as.numeric(table(data$id)) # number of meaurements
  data_jags$W<- data[,"w_obs"]
  data_jags$N<- sum(data_jags$J_i)
  data_jags$id<- data[,"id"]
  
  ### Bernstein Basis
  
  for(index_m in 1:length(m)){
    data_jags$m<- m[index_m] # order of BP
    data_jags$basis<- BP_Basis(t = data[,"t_ij"], m = data_jags$m)
    results_directory<- paste0(directory, "beta", a, b, "/data", r, "/m", data_jags$m, "/") ; dir.create(path = results_directory, showWarnings = FALSE)
    
    ### Prior specifications
    data_jags$I<- diag(x = 1, nrow = data_jags$m)
    data_jags$mu_mu_ksi<- rep(x = 0, times = data_jags$m) ; data_jags$Tau_mu_ksi<- diag(x = 1/100, nrow = data_jags$m)
    
    ###---------------------------
    
    ### Parameters
    parameters<- c("sigma2_e", "mu_ksi", "Sigma_ksi")
    
    ### Initial values
    initial_values<- function(){list(tau_e = rexp(n = 1, rate = 1),
                                     ksi = matrix(data = 1, nrow = data_jags$n, ncol = data_jags$m),
                                     mu_ksi = rep(x = 1, times = data_jags$m),
                                     Tau_ksi = rWishart(n = 1, df = data_jags$m, Sigma = data_jags$m*diag(nrow = data_jags$m))[,,])
    }
    
    
    
    ### Chain specifications
    burn_in<- 50000
    lag<- 10
    sample_size<- 5000 # tamanho da amostra a posteriori final
    number_iterations<- burn_in + (sample_size)*lag
    
    # Carregando o modelo:
    model<- file.path(paste0(model_directory, "model.R"))
    # file.show(model) 
    
    #-------------------
    # Compiling on JAGS
    t1<- Sys.time()
    output<- jags(data = data_jags, inits = initial_values, parameters = parameters, model.file = model, n.chains = 2,
                  n.iter = number_iterations, n.burnin = burn_in, n.thin = lag, working.directory = directory, DIC = FALSE)
    t2<- Sys.time()
    t2 - t1
    #----------------
    
    ### Saving time of posterior computing
    # write.table(x = t(c(t1, t2)), file = paste0(results_directory, "running_time.txt"))
    
    ###-----------------------------------------------
    
    ### Variable names
    names_mu_xi<- paste0("mu_ksi[", 1:data_jags$m, "]")
    names_Sigma_xi<- paste0("Sigma_ksi[", expand.grid(1:data_jags$m, 1:data_jags$m)[,1], ",", expand.grid(1:data_jags$m, 1:data_jags$m)[,2],
                            "]")
    names_sigma2_e<- "sigma2_e"
    
    ### Saving posterior sample
    posterior_sample<- as.mcmc(x = output, start = burn_in, end = number_iterations, thin = lag)
    posterior_sample_mu_xi<- as.matrix(posterior_sample[[1]][,names_mu_xi])
    save(posterior_sample_mu_xi, file = paste0(results_directory, "posterior_sample_mu_xi.Rdata"))
    # save(posterior_sample, file = paste0(results_directory, "posterior_sample.Rdata"))
    
    ### Plots
    # mcmcplot(mcmcout = posterior_sample, parms = parameters[-length(parameters)], dir = results_directory, filename = "mcmc_output",
    #          greek = TRUE)
    
    ### Point estimates
    # posterior_summary<- MCMCsummary(object = posterior_sample, params = parameters, digits = 4, n.eff = TRUE, func = mode,
    #                                 func_name = "mode")
    # summary(posterior_summary[,"n.eff"])
    # which(posterior_summary[,"n.eff"] < sample_size) ; length(which(posterior_summary[,"n.eff"] < sample_size))
    # which(posterior_summary[,"Rhat"] < 1.2) ; length(which(posterior_summary[,"Rhat"] < 1.2))
    
    ### Estimates
    # parameters2<- c(names_sigma2_e, names_mu_xi, names_Sigma_xi)
    # posterior_estimates<- cbind(posterior_summary, HPDinterval(obj = posterior_sample[[1]][,parameters2]))[,c("mean", "50%", "mode", "sd",
    #                                                                                                           "lower", "upper", "Rhat",
    #                                                                                                           "n.eff")]
    # colnames(posterior_estimates)<- c("mean", "median", "mode", "sd", "lower", "upper", "Rhat", "n.eff")
    # write.table(x = posterior_estimates, file = paste0(results_directory, "posterior_estimates.txt"))
    
    ###-----------------------------------------------
    ###   Analysis
    ###-----------------------------------------------
    
    ###-----
    ### 
    # time_plot<- seq(from = 0, to = 1, length.out = 150)
    # basis_plot<- BP_Basis(t = time_plot, m = data_jags$m)
    # mean_curve<- t(unlist(sapply(1:sample_size, function(l)(basis_plot%*%posterior_sample[[1]][,names_mu_xi][l,]))))
    ###-----
    
    
    ###-----
    log_likelihood<- matrix(data = NA, nrow = data_jags$N, ncol = sample_size) 
    
    for(i in 1:data_jags$n){
      position<- which(data_jags$id == i)
      for(o in 1:sample_size){
        sigma2_e<- posterior_sample[[1]][o,names_sigma2_e]
        posterior_mu_W<- as.numeric(data_jags$basis[position,]%*%posterior_sample[[1]][o,names_mu_xi])
        posterior_Sigma_xi<- matrix(data = posterior_sample[[1]][o,names_Sigma_xi], nrow = data_jags$m, ncol = data_jags$m, byrow = TRUE)
        posterior_Sigma_W<- (data_jags$basis[position,]%*%posterior_Sigma_xi)%*%t(data_jags$basis[position,])
        for(q in 1:length(position)){
          log_likelihood[position[q],o]<- dnorm(x = data_jags$W[position[q]], mean = posterior_mu_W[q],
                                                sd = sqrt(posterior_Sigma_W[q,q] + sigma2_e), log = TRUE)
        }
        
      }
    }
    ###-----
    
    ###-----
    ### Calculating comparison measures
    lpml<- CPO(loglik = log_likelihood, n = sample_size)[1]
    waic<- WAIC(loglik = log_likelihood, n = sample_size)[1]
    write.table(x = rbind(lpml, waic), file = paste0(results_directory, "comparison_measures.txt"))
  }
  
  
  
}


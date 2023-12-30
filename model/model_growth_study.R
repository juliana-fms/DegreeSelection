model{
  
  for(l in 1:N){ # long format
    
    Y[l] ~ dnorm(X[l], tau_e) # observed values
    X[l]<- inprod(basis[l,1:m], xi[id[l],1:m]) # mean
    log_lik[l]<- dnorm(Y[l], X[l], tau_e) # log-likelihood
    
  }
  
  ###---------------------------
  
  ### Transformations
  Sigma_xi<- inverse(Tau_xi)
  sigma2_e<- 1 / tau_e
  
  ###---------------------------
  
  ### Prior distributions
  tau_e ~ dgamma(0.01, 0.01) # E(tau_e) = 1 ; Var(tau_e) = 100
  
  for(i in index_boys){
    xi[i, 1:m] ~ dmnorm(mu_xi_boys[1:m], Tau_xi[1:m,1:m])
  }
  for(i in index_girls){
    xi[i, 1:m] ~ dmnorm(mu_xi_girls[1:m], Tau_xi[1:m,1:m])
  }
  
  mu_xi_boys[1:m] ~ dmnorm(mu_mu_xi[1:m], Tau_mu_xi[1:m,1:m])
  mu_xi_girls[1:m] ~ dmnorm(mu_mu_xi[1:m], Tau_mu_xi[1:m,1:m])
  Tau_xi ~ dwish(m*I, m+2)
  
}


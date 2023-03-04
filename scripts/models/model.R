model{
  
  for(l in 1:N){ # long format
    
    W[l] ~ dnorm(X[l], tau_e) # observed values
    X[l]<- inprod(basis[l,1:m], ksi[id[l],1:m]) # mean
    
  }
  
  ###---------------------------
  
  ### Transformations
  Sigma_ksi<- inverse(Tau_ksi)
  sigma2_e<- 1 / tau_e
  
  ###---------------------------
  
  ### Prior distributions
  tau_e ~ dgamma(0.01, 0.01) # E(tau_e) = 1 ; Var(tau_e) = 100
  
  for(i in 1:n){
    ksi[i, 1:m] ~ dmnorm(mu_ksi[1:m], Tau_ksi[1:m,1:m])
  }
  
  mu_ksi[1:m] ~ dmnorm(mu_mu_ksi[1:m], Tau_mu_ksi[1:m,1:m])
  Tau_ksi ~ dwish(m*I, m+2)
  
}


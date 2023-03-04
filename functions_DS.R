mode <- function (x) {
  if (!is.numeric(x)) {
    return(NA)
    stop('Not numeric x vector')
  } else {
    fit <- density(x)
    fit$x[which.max(fit$y)]
  }
}

###----------------------------------------------

hmean<- function(x){
  return(length(x) / sum(1 / x))
}

### --------------------------------------------------------

CPO <- function(loglik, n)
{ 
  lik <- apply(loglik, 2, exp)
  CPO <- apply(lik, 2, hmean)
  LPML <- sum(log(CPO))
  aLPML <- mean(log(CPO))
  return(c(LPML, aLPML)) 
}

### --------------------------------------------------------
DIC <- function(loglik, n)
{
  D <- apply( -2*loglik, 1, sum);
  pD <- 0.5*var(D)
  DIC <- mean(D) + pD;
  return( matrix(c(DIC, pD), ncol=2) ) 
}

### --------------------------------------------------------
WAIC <- function(loglik, n)
{
  lpd <- sum( log( apply(exp(loglik), 2, mean) ) )
  pD <- sum( apply(loglik,  2, var) )
  WAIC <- lpd - pD
  return( matrix(c(WAIC, pD),ncol=2) ) 
}

#---------------------------------

grade<- function(tempos, eventos, numero.intervalos = NULL){
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  #	Funcao que retorna um vetor com os limites dos intervalos														#
  #	Entradas: tempos -> vetor com os tempos de falha/censura 														#
  #		eventos -> vetor indicando a ocorrencia de falha/censura													#
  #		numero.intervalos -> um valor, caso o numero de intervalos seja fixado; NULL caso seja o numero de falhas distintas observadas 	#
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  
  ordem<- order(tempos)
  tempos<- tempos[ordem]
  eventos<- eventos[ordem]
  temposfalha.distintos<- unique(tempos[eventos == 1])
  if(is.null(numero.intervalos)){
    numero.intervalos<- length(temposfalha.distintos)
  }
  # m = numero maximo de intervalos
  m<- length(temposfalha.distintos) 
  if(numero.intervalos > m){
    a<- c(0, unique(tempos[eventos == 1]))
    a[length(a)]<- Inf
  }else{
    b<- min(m, numero.intervalos)
    k1<- trunc(m / b)
    r<- m - b*k1
    k2<- k1 + 1 
    idf1<- seq(k1, (b-r)*k1, k1)
    idf2<- sort(seq(m, max(idf1), -k2))
    idf<- unique(c(idf1, idf2))
    a<- c(0, temposfalha.distintos[idf])
    a[length(a)]<- Inf
  }
  return(a = a) 
}

g<- function(t, tau, m){
  #-----------------------------------------------------------------------------------------------------#
  # Funcao que retorna um vetor representando a funcao de base g do polinomio de Bernstein - funcao de  #
  # densidade da Beta / tau.                                                                            #
  # Entradas: t -> vetor com os tempos de falha/censura                                                 #
  #           tau -> tau - maximo entre os tempos                                                       #
  #           m -> grau do polinomio                                                                    #
  #-----------------------------------------------------------------------------------------------------#
  
  funcao.g<- matrix(data = NA, nrow = length(t), ncol = m)
  for(k in 1:m){
    funcao.g[,k]<- dbeta(x = t/tau, shape1 = k, shape2 = m - k + 1) / tau
  }
  return(funcao.g)
}

#--------------------------

G<- function(t, tau, m){
  #---------------------------------------------------------------------------------------------------#
  # Funcao que retorna um vetor representando a funcao de base G do polinomio de Bernstein - funcao   #
  # de distribuicao acumulada da Beta / tau.                                                          #
  # Entradas: t -> vetor com os tempos de falha/censura                                               #
  #           tau -> tau - maximo entre os tempos                                                     #
  #           m -> grau do polinomio                                                                  #
  #---------------------------------------------------------------------------------------------------#
  
  funcao.G<- matrix(data = NA, nrow = length(t), ncol = m)
  for(k in 1:m){
    funcao.G[,k]<- pbeta(q = t/tau, shape1 = k, shape2 = m - k + 1, lower.tail = TRUE)
  }
  return(funcao.G)
}

#--------------------------

BP_Basis<- function(t, m){
  return(t(sapply(t, function(x)(sapply(1:m, function(k)(dbeta(x = x, shape1 = k, shape2 = m - k + 1))) / m))))
}

###----------------------------

degree_elevation_r1<- function(xi){
  ###-------------------------------------------------------###
  ### Function to use degree elevation property when r = 1  ###
  ### xi -> matrix of Bernstein coefficients of degree m    ###
  ### Returns a matrix with (m+1) Bernstein coefficients    ###
  ###-------------------------------------------------------###
  m<- ncol(xi)
  xi_elevation<- matrix(data = NA, nrow = nrow(xi), ncol = m + 1)
  xi_elevation[,1]<- xi[,1] ; xi_elevation[,m+1]<- xi[,m]
  xi_elevation[,2:m]<- sapply(2:m, function(k)(((k / (m+1))*xi[,k]) + ((1 - (k / (m+1)))*xi[,k])))
  return(xi_elevation)
}





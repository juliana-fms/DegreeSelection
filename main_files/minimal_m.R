###---------------------------------------------------###
###   File #2: minimum degree                         ###
###   Longitudinal data - Bernstein Polynomials       ###
###   Author: Juliana Freitas de Mello e Silva        ###
###---------------------------------------------------###

### Removing objects already allocated
rm(list = ls())

### Probability function of M
p_M_Beta<- function(m, teta1, teta2){
  prob<- ((pbeta(q = (m-2)/(m-1), shape1 = teta1, shape2 = teta2) - pbeta(q = 1/(m-1), shape1 = teta1, shape2 = teta2))^2) -
    (pbeta(q = (m-3)/(m-2), shape1 = teta1, shape2 = teta2) - pbeta(q = 1/(m-2), shape1 = teta1, shape2 = teta2))^2
  return(prob)
}

###-----------------------------
### Girls
###-----------------------------

a_girls<- 5 ; b_girls<- 3
par(mar = c(5, 5, 3, 2) + 0.1)
curve(expr = dbeta(x = x, shape1 = a_girls, shape2 = b_girls),
      from = 0, to = 1,
      xlab = "Age", ylab = paste0("Beta(", round(x = a_girls, digits = 4), ",", round(x = b_girls, digits = 4), ")"), cex.lab = 1.2,
      xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 18, by = 0.5) / 18, labels = seq(from = 0, to = 18, by = 0.5))

### Distribution of minimum M for girls given (U1,U2)
p_M_girls<- p_M_Beta(m = 4:30, teta1 = a_girls, teta2 = b_girls) ; names(p_M_girls)<- 4:30
round(x = p_M_girls, digits = 4)
which.max(p_M_girls) # max
round(x = cumsum(x = p_M_girls), digits = 4) # cumulative 

###-----------------------------
### Boys
###-----------------------------
a_boys<- 20 ; b_boys<- 7
par(mar = c(5, 5, 3, 2) + 0.1)
curve(expr = dbeta(x = x, shape1 = a_boys, shape2 = b_boys),
      from = 0, to = 1,
      xlab = "Age", ylab = paste0("Beta(", round(x = a_boys, digits = 4), ",", round(x = b_boys, digits = 4), ")"), cex.lab = 1.2,
      xaxt = "n")
axis(side = 1, at = seq(from = 0, to = 18, by = 0.5) / 18, labels = seq(from = 0, to = 18, by = 0.5))

### Distribution of minimum M for boys given (U1,U2)
p_M_boys<- p_M_Beta(m = 4:30, teta1 = a_boys, teta2 = b_boys) ; names(p_M_boys)<- 4:30
round(x = p_M_boys, digits = 4)
which.max(p_M_boys) # max
cumsum(x = p_M_boys) # cumulative

### Next, you can go to file #3


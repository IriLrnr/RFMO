library(ggplot2)

EquilPop <- function(r, h){
  H <- sum(h)
  h_prop <- h/H
  n_star <- max(0,1 - H/r)
  Y <- H*n_star
  y <- Y*h_prop
  return(list(y=y, n=n_star))
}

ConservationFunction <- function(n, R) {
  C <- R*(exp(2*n)-1)
  return(C)
}

ProfitFunction <- function(y_i, h_i, c) {
  p_i <- pmax(0, y_i - c*h_i)
  return(p_i)
}

UtilityFunction <- function(P, C, w) {
  B <- C^w * P^(1-w)
  return(B)
}

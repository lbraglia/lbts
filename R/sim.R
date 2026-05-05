#'@export
sim_ma1 <- function(n=100, omega=0, theta=0.5, sigma2=1){
  e <- sqrt(sigma2) * rnorm(n) # epsilon_t sim N(0, sigma^2)
  x <- omega # x[1]=omega: initialization of the process
  for(t in seq(1, n-1)){
    x[t+1] <- omega  + theta * e[t] + e[t+1]
  }
  x
}

#'@export
sim_ar1 <- function (n = 100, omega = 0, psi = 0.5, sigma2 = 1) 
{
  e <- sqrt(sigma2) * rnorm(n) # white noise N(0, sigma2)
  x <- omega/(1-psi) # x[1] initialisation matters when omega is nonzero
  for (t in seq(1, n - 1)) {
    x[t+1] <- omega + psi * x[t] + e[t+1]
  }
  x
}


#'@export
sim_rw <- function (n = 100, x0 = 0, xt = rnorm) 
{
  x <- x0
  xts <- xt(n=n) ## iid
  rw <- cumsum(c(x0, xts))
  rw[-length(rw)]
}


#'@export
sim_arma11 <- function(n=100, omega=0, theta=0.5, psi = 0.5, sigma2=1){
  e <- sqrt(sigma2) * rnorm(n) # epsilon_t sim N(0, sigma^2)
  x <- omega/(1-psi) # initialization
  for (t in seq(1, n-1)){
    x[t+1] = omega + psi * x[t]   + theta * e[t]    + e[t+1]
  }
  x
}

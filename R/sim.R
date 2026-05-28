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
sim_ar1 <- function (n = 100, omega = 0, phi = 0.5, sigma2 = 1) 
{
  e <- sqrt(sigma2) * rnorm(n) # white noise N(0, sigma2)
  x <- omega/(1-phi) # x[1] initialisation matters when omega is nonzero
  for (t in seq(1, n - 1)) {
    x[t+1] <- omega + phi * x[t] + e[t+1]
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
sim_arma11 <- function(n=100, omega=0, theta=0.5, phi = 0.5, sigma2=1){
  e <- sqrt(sigma2) * rnorm(n) # epsilon_t sim N(0, sigma^2)
  x <- omega/(1-phi) # initialization
  for (t in seq(1, n-1)){
    x[t+1] = omega + phi * x[t]   + theta * e[t]    + e[t+1]
  }
  x
}



#'@export
sim_ar1noise <- function(n=100,
                         sigma_epsilon=1, # variance error term noise
                         sigma_eta=1, # variance error term signal
                         phi=0.5, # AR phi parameter
                         mu1=0, # unconditional mean of an AR(1) process with no intercept (its mu[1])
                         y1=0, # starting point of series
                         return_mu=TRUE
                         ){

  y <- y1
  mu <- mu1
  eta <- sigma_eta * rnorm(n) # N(0,sigma_eta^2)
  epsilon <- sigma_epsilon * rnorm(n)

  for(t in 1:(n-1)){
    # we start updating the prediction of next signal with initialization parameters
    mu[t+1] <- phi * mu[t] + eta[t]
    # obtain the observed signal + noise
    y[t+1] <-  mu[t+1] + epsilon[t+1]
  }

  if (return_mu) list("y"=y, "mu"=mu) else y
}

################
##### Estimation
################

#' @export
uDCS_t_model_estimator <- function(y, theta_0){
  
  ## ###Take T
  ## T <- length(dati)
  
  ## ###Parameter Selections Dynamic Location
  ## omega <- param[1]
  ## phi   <- param[2]
  ## k     <- param[3]
  ## varsigma <- param[4]
  ## nu       <- param[5]

  ## ###Create a vector with the parameters
  
  ## Take Bounds
  ## t  <- c(omega, phi, k, varsigma, nu)
  lower <- c(-Inf, -0.999, -2, 1e-05, 2.099)
  upper <- c( Inf,  0.999,  2, Inf, 300)
  
  optimizer <- suppressWarnings(nlminb(start = theta_0, objective = interprete_uDCS_t_model, 
                                       dati  = y, gradient = NULL, 
                                       control = list(trace = 0), hessian = NULL,
                                       lower = lower, upper = upper))
  
  ## #------> Save the optimized parameters Dynamic Location
  ## omega_opt <- optimizer$par[1]  
  ## phi_opt <- optimizer$par[2]
  ## k_opt  <- (optimizer$par[3])
  
  ## varsigma_opt <- optimizer$par[4]
  ## nu_opt       <- optimizer$par[5]
  
  ## ###Create a vector with ALL the optimized parameters
  ## theta_opt <- c(omega_opt, phi_opt, k_opt, varsigma_opt, nu_opt)

  ## ###Create a list with ALL the optimized parameters
  ## theta_list <- list(omega = omega_opt,
  ##                    phi = phi_opt,
  ##                    k  = k_opt,
  ##                    varsigma = varsigma_opt,
  ##                    nu    = nu_opt)
  
  ######################
  ####### OUTPUT #######
  ######################
  
  ## #------> Some detail
  ## Elapsed_Time <- Sys.time() - Start
  ## print(paste("Elapsed Time: ", toString(Elapsed_Time)))
  
  ## ###Make List
  ## out <- list(theta_list = theta_list,
  ##             theta      = theta_opt,
  ##             optimizer  = optimizer)
  
  ## return(out)
  optimizer$par
}


################
##### Interprete
################

#' @export
interprete_uDCS_t_model <- function(dati, param){
    
  ###Take T
  T <- length(dati)
  
  ###Parameter Selections Dynamic Location
  omega <- param[1]
  phi   <- param[2]
  k     <- param[3]
  
  varsigma <- param[4]
  nu       <- param[5]
    
   ###Create a new vector with the parameters
   theta_new <- c(omega, phi, k, varsigma, nu)
    
    #------> Fitness Functions
    fitness <- uDCS_t_model_filter(dati, theta_new)$Log_Likelihood
    
    if(is.na(fitness) | !is.finite(fitness)) fitness <- -1e10
    if(fitness != fitness) fitness <- -1e10
    
  return(-fitness)
} 

###############
##### Filtering
###############

#' @export
uDCS_t_model_filter <- function(y, theta){
  
  ###Take T
  T <- length(y)
  
  ###Define LogLikelihoods
  dloglik <- array(data = NA, dim = c(T))
  loglik  <- numeric()
  
  ###Parameter Selections Dynamic Location
  omega <- theta[1]
  phi   <- theta[2]
  k     <- theta[3]
  
  varsigma <- theta[4]
  nu       <- theta[5]
  
  ###Define Dynamic Location and Innovations
  mu_t <- array(data = NA, dim = c(T+1))
  u_t  <- array(data = NA, dim = c(T))
  
  ###Initialize Dynamic Location
  mu_t[1]   <- (omega)
    
    ###Initialize Likelihood
    dloglik[1] <- uSTDT_uDCS_t(y[1], mu_t[1], varsigma = varsigma, nu = nu, log = TRUE)
    loglik     <- dloglik[1]
    
    for(t in 2:(T+1)) {
      ###Dynamic Location Innovations
      u_t[t-1] <- martingale_diff_u_t(y[t-1], mu_t[t-1], varsigma, nu)
      ###Updating Filter                    
      mu_t[t]   <- omega + phi * (mu_t[t-1] - omega) + k * u_t[t-1]
    
      if(t < (T+1)){
        ###Updating Likelihoods
        dloglik[t] <- uSTDT_uDCS_t(y[t], mu_t = mu_t[t], varsigma = varsigma, nu = nu, log = TRUE)
        loglik     <- loglik + dloglik[t]
      }
    }

  ######################
  ####### OUTPUT #######
  ######################
  mu_t <- ts(mu_t, start = start(y), frequency = frequency(y))
  u_t  <- ts(u_t, start = start(y), frequency = frequency(y))
  
  ###Make List
  out <- list(Dynamic_Location = mu_t,
              Innovation_u_t   = u_t,
              Log_Densities_i  = dloglik,
              Log_Likelihood   = loglik)
  
  return(out)
}

############################################
####### ADDITIONAL FUNCTIONS #######
############################################

#' @export
martingale_diff_u_t <- function(y, mu_t, varsigma, nu){
  
  u_t <- c((1 / (1 + (y - mu_t)^2/(nu*varsigma)) * (y - mu_t)))
  
  return(u_t)
}

#' @export
uSTDT_uDCS_t <- function(y, mu_t, varsigma, nu, log = TRUE){
  
  ulpdf <- (lgamma((nu + 1) / 2) - lgamma(nu / 2) - (1/2) * log(varsigma) -
              (1/2)  * log(pi * nu) - ((nu + 1) / 2) * log(1 + (y - mu_t)^2 / (nu*varsigma) ))
  
  if(log != TRUE){
    ulpdf <- exp(ulpdf)
  } 
  
  return(ulpdf)
}

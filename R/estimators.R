#' @export
KF <- function(y, theta_0 = c("phi"=0.8, "sigma_epsilon"=1, "sigma_eta"=1), plot=TRUE){

  KF_pred <- function(y, phi, sigma_epsilon, sigma_eta){
    n <- length(y)
    # allocate space 
    mu_pred <- 0  # this will be updated as mu_{t|t-1}: initialize it
    P <- 1   # this will be updated as P_{t|t-1}: initialize it
    v <- rep(NA, n) # the innovation error 
    K <- 0   # the Kalman gain 
    F <- 0   # the conditional variance of v_t 
    dllk <- rep(NA, n)    #the log-likelihood value 
    llk <- 0 
    # the recursion 
    for(t in 1:(n-1)){
      v[t] = y[t] - mu_pred[t]; 
      F[t] = P[t] + sigma_epsilon^2;
      K[t] = (phi * P[t])/F[t]
      P[t+1] = phi^2 * P[t] + sigma_eta^2 - K[t]*F[t]*K[t]
      mu_pred[t+1] = phi * mu_pred[t] + K[t]*v[t]
      dllk[t] = - 0.5 * log(F[t] + (v[t]^2/F[t]))
      llk  = llk + dllk[t]
    }
    list("mu_pred" = mu_pred, "llk"=llk)
  }

  KF_loglikelihood <- function(par, y){
    phi <- par[1]
    sigma_epsilon <- par[2]
    sigma_eta <- par[3]
    - KF_pred(y, phi=phi, sigma_epsilon=sigma_epsilon, sigma_eta=sigma_eta)$llk
  }

  ## parameters optimization
  ## hat_theta <- nlminb(start = theta_0, objective = KF_loglikelihood, y = y)$par
  hat_theta <- optim(par = theta_0, fn = KF_loglikelihood, y = y)$par

  ## predicted mu by KF
  kf_mu <- KF_pred(y=y,
                   phi=hat_theta["phi"],
                   sigma_epsilon=hat_theta["sigma_epsilon"],
                   sigma_eta=hat_theta["sigma_eta"])$mu_pred

  if (plot) {
    lbts::mts_plot(list(y, kf_mu), col = c("grey", "blue"))
    legend("topleft", col = c("grey", "blue"), lty=1, legend = c("observed y", "KF mu"))
  }
  
  invisible(list("hat_theta" = hat_theta, "y" = y, "KF_mu" = kf_mu))
       
}

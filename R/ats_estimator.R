#' @export
ats_estimator <- function(y,
                          kf_theta0 = c("phi" = 0.8, "sigma_epsilon" = 1, "sigma_eta" = 1),
                          sc_theta0 = c("omega" = 0, "phi" = 0.8, "k" = 0.8, "varsigma"=1, "nu"=8))
{
  ## Kalman filter
  kf <- KF(y, theta_0=kf_theta0, plot=FALSE)
  ## Score driven model is basically the same function but with suppressed
  ## messaging and other minor things changed
  sd_hat_theta  <- uDCS_t_model_estimator(y, sc_theta0)
  sd_filter <- uDCS_t_model_filter(y, sd_hat_theta)
  sct <- c(list("hat_theta" = sd_hat_theta), sd_filter)
  ## wrapping all together
  rval <- list("y"=y,"kf"=kf, "sct"=sct)
  class(rval) <- "ats_estimates"
  rval
}


#' @export
print.ats_estimates <- function(x)
  print(list("kf_hat_theta" = x$kf$hat_theta,
             "sct_hat_theta" = x$sct$hat_theta))


#' @export
plot.ats_estimates <- function(x, ...){
  # browser()
  ts_list <- list(x$y, x$kf$KF_mu, x$sct$Dynamic_Location)
  cols <- c("grey", "red", "blue")
  ltys <- c("solid", "solid", "dashed")
  legend <- c("observed y", "KF mu", "SDt loc")
  lbts::mts_plot(
    lapply(ts_list, as.numeric), # normalize
    col = cols,
    lty = ltys,
    ...
  )
  legend("topleft", col = cols, lty=ltys, legend = legend)
}

#' @export
residuals.ats_estimates <- function(x, lag.max=50){
  ## browser()
  kf_innovation_error <- as.numeric(x$kf$v)
  sct_innovation_error <- as.numeric(x$y - x$sct$Dynamic_Location)
  plTS <- function(x, ...) ts.plot(x, xlab="Time", ...)
  plACF <- function(x) acf(x,lag.max = lag.max, main="")
  plPACF <- function(x) pacf(x,lag.max = lag.max, main="")
  qq <- function(x) {
    test <- lbmisc::pretty_pval(shapiro.test(x)$p.value, equal=TRUE)
    subtitle <- sprintf("(Shapiro's p%s, Kurtosis=%.2f)", test, kurtosis(x))
    qqnorm(x)
    qqline(x)
    title(sub=subtitle)
  }
  par(mfrow=c(2,3))
  # check WN
  plTS(kf_innovation_error, main = "KF innovation error")
  plACF(kf_innovation_error)
  qq(kf_innovation_error)
  # check T
  plTS(sct_innovation_error, main = "SC-t innovation error")
  plACF(sct_innovation_error)
  qq(sct_innovation_error)
}

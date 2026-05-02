# Plot of a single time series with standard diagnostic
#'@export
ts_plot <- function(x, lag.max=50, ...){
  par(mfrow=c(1,3))
  ts.plot(x, xlab="Time", ...)
  acf(x,lag.max = lag.max, drop.lag.0 = FALSE, main="")
  pacf(x,lag.max = lag.max, main="")
}

# Plotting of multiple time series
#'@export
mts_plot <- function(x, col=lbmisc::col2hex("black", 0.1), ...){
  if (is.list(x)) x <- do.call(cbind, x)
  ts.plot(x, xlab="Time", col=col, ...)
}


# Stick plot for theoretical ACF
#'@export
stick_plot <- function(y, x=seq(0,length(y)-1), xlab="Time", ...){
  plot(x = x, y = y, type = 'h', ...)
  abline(h=0)
}

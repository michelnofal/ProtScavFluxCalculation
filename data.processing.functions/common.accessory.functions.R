### COMMON ACCESSORY FUNCTIONS

# generate sequence of xvalues to plot and predict y values using fit
generate.predictions <- function (fit, dat) { 
  max.x <- max(dat$time.pt)
  xvals <- seq(0, max.x, max.x/200)
  data.frame(time.pt=xvals, predicted=predict(fit, data.frame(time.pt=xvals)))
}

# integrate exponential curve using model parameters
integrate.exponential <- function (intercept, slope, max) {
  intercept <- as.numeric(intercept)
  slope <- as.numeric(slope)
  max <- as.numeric(max)
  as.numeric(integrate(function (t) exp(intercept + slope*t), 0, max)[1])
}


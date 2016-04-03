### I BELIEVE THIS FILE IS UNNEEDED (replaced by PCV.curvefit.plot.and.compare.R)

# load common accessory functions
# includes generate.predictions and integrate.exponential
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/common.accessory.functions.R")

# accessory function: converts growth rate into doubling time
growthRate.to.doublingTime <- function (growthRate) {
  1/(log(exp(growthRate),2))
}

###
### FUNCTIONS TO 1) PROCESS 2) PLOT GROWTH DATA
###

# 1 - calculate growth predictions at measured time points (to plot points) - OUTPUT: $preds
# 2 - generate predictions at short, regular time intervals (to plot curve) - OUTPUT: $curve
# 3 - calculate fitted growth curve parameters and area under growth curve - OUTPUT: $params
# need to import name for growth data - for eventual comparison of multiple growth curves
process.growth.dat <- function(data, name) {
  output.list <- list() # any object to be returned will be saved in this list
  growth.dt <- tbl_df(data)
  
  growth.fits <- growth.dt %>% do(fit=nls(pcv ~ exp(intercept + rate*time.pt), data=., start=list(intercept=0, rate=.1)), data=(.))
  growth.augment <- growth.fits %>% do(augment(.$fit[[1]], .$data[[1]]))
  growth.preds <- growth.augment %>% mutate(predicted = .fitted) %>% select(time.pt, pcv, predicted)
  growth.preds$data <- name
  output.list$preds <- growth.preds
  
  growth.preds.curve <- growth.fits %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))
  growth.preds.curve$data <- name
  output.list$curve <- growth.preds.curve
  
  growth.tidy <- growth.fits %>% do(tidy(.$fit[[1]]))
  curve.params <- data.frame(intercept = growth.tidy[growth.tidy$term == "intercept",]$estimate, 
                             rate = growth.tidy[growth.tidy$term == "rate",]$estimate)
  curve.params$doubling.time <- growthRate.to.doublingTime(curve.params$rate)
  
  growth.max <- growth.dt %>% do(data.frame(max = max(.$time.pt)))
  curve.params$max <- growth.max
  growth.integral <- integrate.exponential(curve.params$intercept, curve.params$rate, curve.params$max)
  curve.params$growth.integral <- growth.integral
  curve.params$data <- name
  output.list$params <- curve.params
  
  return(output.list)
}

# plot processed growth data
growth.plotter <- function(dat, curve) {
  growth.dt <- dat
  growth.preds.curve <- curve
  
  growth.max <- max(dat$time.pt)
  
  growth.plot <- ggplot() + 
    geom_point(data=growth.dt, aes(x=time.pt, y=pcv)) + 
    geom_line(data=growth.preds.curve, aes(x=time.pt, y=predicted)) + 
    scale_x_continuous(breaks = seq(0, growth.max, growth.max/3), limits = c(0, growth.max)) + 
    scale_y_continuous(limits = c(0, ceiling(max(growth.dt$pcv + 0.5)))) +
    xlab("Time (hours)") + ylab("Packed Cell Volume (uL)")
  growth.plot
}


###
### Example
###

# setwd("/Users/Michel/Desktop/Research/code/")
# growth.dat_KRPCA_May2015 <- read.table("ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/krpc_HF_lab_timecourse_pcv.csv", header=TRUE, sep=",")
# proc.dat_KRPCA_May2015 <- process.growth.dat(growth.dat_KRPCA_May2015, "KRPCA_May2015_HFwBSA")
# KRPCA_May2015_plot <- growth.plotter(proc.dat_KRPCA_May2015$preds, proc.dat_KRPCA_May2015$curve)

proc.dat_KRPCA_May2015 <- process.growth.dat(growth.dat_KRPCA_May2015, "KRPCA_May2015_HFwBSA")
proc.dat_KRPCA_May2015_other <- growth.func.new(growth.dat_KRPCA_May2015, "KRPCA_May2015_HFwBSA")

# load common accessory functions
# includes generate.predictions and integrate.exponential
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/common.accessory.functions.R")

# convert growth rate into doubling time
growthRate.to.doublingTime <- function (growthRate) {
  1/(log(exp(growthRate),2))
}

# process pcv data
growth.func.new <- function(data, name) {
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

# plot growth curve 
# input is output of growth.func.new
growth.curve.plotter <- function(processed.growth.data) {
  growth.dt <- processed.growth.data$preds
  growth.preds.curve <- processed.growth.data$curve
  
  growth.max <- max(growth.dt$time.pt)
  
  growth.plot <- ggplot() + 
    geom_point(data=growth.dt, aes(x=time.pt, y=pcv)) + 
    geom_line(data=growth.preds.curve, aes(x=time.pt, y=predicted)) + 
    scale_x_continuous(breaks = seq(0, growth.max, growth.max/3), limits = c(0, growth.max)) + 
    scale_y_continuous(limits = c(0, ceiling(max(growth.dt$pcv + 0.5)))) +
    xlab("Time (hours)") + ylab("Packed Cell Volume (uL)")
  growth.plot
}

# compare growth curves
growth.comparer.new <- function (...) {
  output.list <- list()
  
  data.to.plot <- data.frame()
  curve.to.plot <- data.frame()
  agg.growth.dat <- data.frame()
  for(each in list(...)) {
    data.to.plot <- rbind(data.to.plot, each$preds)
    curve.to.plot <- rbind(curve.to.plot, each$curve)
    agg.growth.dat <- rbind(agg.growth.dat, data.frame(doubling.time = each$params$doubling.time,
                                                       growth.integral = each$params$growth.integral,
                                                       data = each$params$data))
  }
  output.list$meta <- agg.growth.dat
  
  growth.max <- max(data.to.plot$time.pt)
  
  growth.plot <- ggplot() + 
    geom_point(data=data.to.plot, aes(x=time.pt, y=pcv, color=data)) + 
    geom_line(data=curve.to.plot, aes(x=time.pt, y=predicted, color=data)) + 
    scale_x_continuous(breaks = seq(0, growth.max, growth.max/3), limits = c(0, growth.max)) + 
    scale_y_continuous(limits = c(0, ceiling(max(data.to.plot$pcv + 0.5)))) +
    xlab("Time (hours)") + ylab("Packed Cell Volume (uL)")
  output.list$plot <- growth.plot
  
  output.list
}



### FUNCTION FOR COMPARISON OF GROWTH CURVES


###
### Example (process & plot)
###

# growth.dat_KRPCA_May2015 <- read.table("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/krpc_HF_lab_timecourse_pcv.csv", header=TRUE, sep=",")
# proc.dat_KRPCA_May2015 <- process.growth.dat(growth.dat_KRPCA_May2015, "KRPCA_May2015_HFwBSA")
# KRPCA_May2015_plot <- growth.plotter(proc.dat_KRPCA_May2015$preds, proc.dat_KRPCA_May2015$curve)

###
### Example (compare)
###

# lipo0.growth <- read.table("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/Data/KRPCA_cholesterol/KRPC_lipo0_pcv.csv", header=TRUE, sep=",")
# lipo0_2.growth <- read.table("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/Data/KRPCA_cholesterol/KRPC_lipo0_startingpcv2_pcv.csv", header=TRUE, sep=",")
# lipo20.growth <- read.table("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/Data/KRPCA_cholesterol/KRPC_lipo20_pcv.csv", header=TRUE, sep=",")
# lipo50.growth <- read.table("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/Data/KRPCA_cholesterol/KRPC_lipo50_pcv.csv", header=TRUE, sep=",")

# lipo0.growth.dat <- growth.func.new(lipo0.growth, "krpc_lipo0_rerun")
# lipo0_2.growth.dat <- growth.func.new(lipo0_2.growth, "krpc_lipo0_2_rerun")
# lipo20.growth.dat <- growth.func.new(lipo20.growth, "krpc_lipo20_rerun")
# lipo50.growth.dat <- growth.func.new(lipo50.growth, "krpc_lipo50_rerun")

# lipo50.growth.plot <- growth.curve.plotter(lipo50.growth.dat)

# lipo.growth.comparison <- growth.comparer.new(lipo0.growth.dat, lipo0_2.growth.dat, lipo20.growth.dat, lipo50.growth.dat)


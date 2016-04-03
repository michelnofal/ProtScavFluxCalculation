# load common accessory functions
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/common.accessory.functions.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/aa.info.R")
# aa.info.R generates default.aa.table ; need to use generate.aa.table() if default aa concentrations are not used

preproc.unlab.aa <- function(unproc.media.data) {
  unproc.media.data %>% mutate(sum = unlabeled + labeled)
}

# convert extracellular amino acid counts to umol in medium
get.extra.umol.new <- function (media.aa.dt, aa.table=default.aa.table, num.mLs=2) {
  if (!("sum" %in% colnames(media.aa.dt))) {
    media.aa.dt <- preproc.unlab.aa(media.aa.dt)
  }
  
  init.extra.aa <- media.aa.dt %>% filter(time.pt == 0) %>% group_by(compound) %>% summarize(mean=mean(sum))
  aa.conv.dt <- merge(init.extra.aa, aa.table, by="compound") 
  aa.conv.dt <- aa.conv.dt %>% mutate(conv.factor = DMEM.conc*(num.mLs*.001)/mean) %>% select(-mean, -DMEM.conc)
  
  extra.sums.dt <- media.aa.dt %>% select(compound, time.pt, sum)
  extra.umol.dt <- merge(extra.sums.dt, aa.conv.dt)
  extra.umol.dt <- extra.umol.dt %>% mutate(umol = sum*conv.factor) %>% select(-sum, -conv.factor)
  return(extra.umol.dt)
}

calc.uptake.rates.new <- function (extra.umol.dt, growth.list, name) {
  output.list <- list()
  
  # make data frame containing pcv predictions and extracellular 
  pcv.v.extra <- merge(unique(growth.list$preds %>% select(-pcv)), extra.umol.dt)

  # make model to predict extracellular aa abundance using predicted pcv
  pve.fits <- pcv.v.extra %>% group_by(compound) %>% do(fit=lm(umol ~ predicted, data=.), data=(.))
  pve.tidy <- pve.fits %>% group_by(compound) %>% do(tidy(.$fit[[1]]))
  pve.tidy$term[pve.tidy$term == "(Intercept)"] = "intercept"
  pve.tidy$term[pve.tidy$term == "predicted"] = "slope"
  # pve.spread <- pve.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
  pve.spread <- pve.tidy %>% select(-(std.error:p.value)) %>% spread(term, estimate)
  
  growth.k <- growth.list$params$rate
  uptake.rates <- pve.spread %>% mutate(uptake.rate = -growth.k*slope) %>% select(compound, uptake.rate)
  uptake.rates$data <- name
  output.list$uptake.rates <- uptake.rates
  
  pve.aug <- pve.fits %>% group_by(compound) %>% do(augment(.$fit[[1]], .$dat[[1]]))
  pve.aug <- pve.aug %>% mutate(pred.umol = .fitted) %>% select(predicted, compound, umol, pred.umol)
  pve.plot <- ggplot(pve.aug) + geom_point(aes(x=predicted,y=umol,color=compound)) + 
    geom_line(aes(x=predicted, y=pred.umol, color=compound)) + geom_hline(yintercept=0) +
    scale_y_continuous(limits=c(0, max(pve.aug$umol)*1.1)) +
    scale_color_hue(c=60, l=70) + 
    xlab("Packed Cell Volume (uL)") + ylab("umol in medium") + 
    facet_grid( ~ compound, labeller=aa.labeller.complete)
  output.list$plot <- pve.plot
  
  return(output.list)
}

# calculate umol unlab amino acid at final time point
get.final.umol.unlab.new <- function (extra.umol.dt, media.aa.dt) {
  frac.dt <- media.aa.dt %>% filter(loc == "Extracellular") %>% select(compound, time.pt, frac)
  extra.umol.unlab <- merge(extra.umol.dt, frac.dt)
  extra.umol.unlab <- extra.umol.unlab %>% mutate(umol.unlab = umol*frac)
  
  max.time.pt <- max(extra.umol.unlab$time.pt)
  final.umol.unlab.dt <- extra.umol.unlab %>% filter(time.pt %in% c(0,max.time.pt)) %>% 
    group_by(compound, time.pt) %>% summarize(mean.umol.unlab=mean(umol.unlab))
  final.umol.unlab.dt$time.pt[final.umol.unlab.dt$time.pt == 0] = "time0"
  final.umol.unlab.dt$time.pt[final.umol.unlab.dt$time.pt != "time0"] = "endtime"
  final.umol.unlab.dt <- final.umol.unlab.dt %>% spread(time.pt, mean.umol.unlab)
  final.umol.unlab.dt <- final.umol.unlab.dt %>% mutate(net.umol.unlab = endtime - time0) %>% select(-time0, -endtime)
  
  return(final.umol.unlab.dt)
}

###
### Example
###

# growth.list.new <- growth.func.new(growth.dat, "whateva")
# media.data <- aa.lab.list$data %>% filter(loc == "Extracellular")
# 
# media.extra.umol <- get.extra.umol.new(media.data)
# 
# final.umol.unlab <- get.final.umol.unlab.new(media.extra.umol, media.data)
# uptake.list <- calc.uptake.rates.new(media.extra.umol, growth.list.new, "whateva")

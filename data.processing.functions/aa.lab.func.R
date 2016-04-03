# note: 08.13.2015 - leaving this function the same as it was in protein_influx_functions_v3.R 
# potential improvement: adding data_name column to facilitate plotting labeling data from different timecourses

source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/common.accessory.functions.R")

nat.C13.frac <- 0.01109
get.exp.unlab.frac <- function(numCarbons) {
  dbinom(0, size=numCarbons, prob=nat.C13.frac)
}
# get.exp.unlab.frac(5)

num.carb.df <- data.frame(compound = c("lysine","phenylalanine","threonine","tyrosine","valine","histidine"), 
                          num.carb = c(6, 9, 4, 9, 5, 6),
                          impure.frac = c(0.023, 0.071, 0.038, 0.023, NA, 0.034)) # from Labeled_medium_exactive_Dec2014
adjust.unlab.and.lab <- num.carb.df %>% mutate(exp.unlab.frac = get.exp.unlab.frac(num.carb), 
                                       adj.unlab.factor = ifelse(compound == "valine", 1, 1/exp.unlab.frac),
                                       adj.lab.factor = ifelse(compound == "valine", 1, 1/(1-impure.frac)))

# left_join(intra.dat, num.carb.df, by="compound")

# 1 - return organized intra- and extracellular data
# 2 - calculate fitted amino acid labeling curve parameters and predictions
# 3 - generate amino acid labeling plots
aa.lab.func <- function (intra.data, extra.data, name, adjust=FALSE) {
  output.list <- list()
  
  if (adjust) {
    intra.data <- left_join(intra.data, adjust.unlab.and.lab, by="compound")
    intra.data <- intra.data %>% mutate(unlabeled = unlabeled*adj.unlab.factor, labeled = labeled*adj.lab.factor)
    extra.data <- left_join(extra.data, adjust.unlab.and.lab, by="compound")
    extra.data <- extra.data %>% mutate(unlabeled = unlabeled*adj.unlab.factor, labeled = labeled*adj.lab.factor)
  }

  intra.dt <- tbl_df(intra.data)
  intra.dt <- intra.dt %>% mutate(sum=unlabeled+labeled, frac=unlabeled/sum, loc="Intracellular")
  
  extra.dt <- tbl_df(extra.data)
  extra.dt <- extra.dt %>% mutate(sum=unlabeled+labeled, frac=unlabeled/sum, loc="Extracellular")
  
  aa.dt <- rbind(intra.dt, extra.dt)
  aa.dt$data <- name
  output.list$data <- aa.dt
  
  aa.fits <- aa.dt %>% group_by(compound, loc) %>% 
    do(fit=nls(frac ~ exp(intercept + rate*time.pt), data=., start=list(intercept=0, rate=.1)), data=(.))
  
  aa.aug <- aa.fits %>% group_by(compound, loc) %>% do(augment(.$fit[[1]], .$data[[1]]))
  aa.tidy <- aa.fits %>% group_by(compound, loc) %>% do(tidy(.$fit[[1]]))
  aa.tidy$estimate <- as.numeric(aa.tidy$estimate)
  aa.tidy.spread <- aa.tidy %>% dplyr::select(compound:estimate) %>% spread(term, estimate)
  aa.tidy.spread$data <- name
  output.list$curve.params <- aa.tidy.spread
  
  aa.preds <- aa.fits %>% group_by(compound, loc) %>% do(generate.predictions(.$fit[[1]], .$dat[[1]]))
  aa.dt$loc <- factor(aa.dt$loc, levels = c("Intracellular","Extracellular"))
  aa.preds$loc <- factor(aa.preds$loc, levels = c("Intracellular","Extracellular"))
  aa.preds$data <- name
  output.list$preds <- aa.preds
  
  aa.plots <- ggplot () + 
    geom_point(data=aa.dt, aes(x=time.pt, y=frac, color=compound)) + 
    geom_line(data=aa.preds, aes(x=time.pt, y=predicted, color=compound)) + geom_hline(yintercept=0) + 
    scale_x_continuous(limits = c(0, max(aa.dt$time.pt)), breaks = seq(0, max(aa.dt$time.pt), max(aa.dt$time.pt)/3)) +
    scale_color_hue(c=60, l=70) + 
    xlab("Time (hours)") + ylab("Unlabeled Fraction") + 
    facet_grid(compound ~ loc, scales="free_y") 
  output.list$plots <- aa.plots
  
  return(output.list)
}

###
### Example
###

# setwd("/Users/Michel/Desktop/Research/code/")
# intra.dat_KRPCA_May2015 <- read.table("ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/Intracellular_data_HFB_processed.csv", header=TRUE, sep=",")
# extra.dat_KRPCA_May2015 <- read.table("ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/Extracellular_data_HFB_processed.csv", header=TRUE, sep=",")
# extra.dat_KRPCA_May2015 <- extra.dat_KRPCA_May2015 %>% filter(time.pt != "0") %>% mutate(time.pt = ifelse(time.pt == "fresh", 0, time.pt))
# extra.dat_KRPCA_May2015 <- extra.dat_KRPCA_May2015 %>% filter(compound %in% c("histidine","lysine","phenylalanine","threonine","tyrosine","valine"))
# extra.dat_KRPCA_May2015$time.pt = as.numeric(extra.dat_KRPCA_May2015$time.pt)

# aa.lab_KRPCA_May2015 <- aa.lab.func(intra.dat_KRPCA_May2015, extra.dat_KRPCA_May2015, "HF_BLAH")

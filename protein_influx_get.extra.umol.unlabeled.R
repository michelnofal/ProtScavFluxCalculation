get.extra.umol.unlabeled.only <- function(aa_table, extra.dt, num_mLs) {
  extra.sums.dt <- extra.dt %>% mutate(sum=unlabeled+labeled)
  init.extra.aa <- extra.sums.dt %>% filter(time.pt == 0) %>% group_by(compound) %>% summarize(mean=mean(sum))
  aa.conv.dt <- merge(init.extra.aa, aa_table, by="compound") 
  aa.conv.dt <- aa.conv.dt %>% mutate(conv.factor = DMEM.conc*(num_mLs*.001)/mean) %>% select(-mean, -DMEM.conc)
  
  extra.sums.dt <- extra.sums.dt %>% select(compound, time.pt, sum)
  extra.umol.dt <- merge(extra.sums.dt, aa.conv.dt)
  extra.umol.dt <- extra.umol.dt %>% mutate(umol = sum*conv.factor) %>% select(-sum, -conv.factor) 
  
  return(extra.umol.dt)
}

# modifications:
# input should only require extracellular data, not intracellular
calc.uptake.rates.unlabeled.only <- function (aa_table, extra.dt, growth.preds, growth.curve.params, num_mLs) {
  output.list <- list()
  
  # convert amino acid 
  extra.umol.dt <- get.extra.umol.unlabeled.only(aa_table, extra.dt, num_mLs)
  
  # make data frame containing pcv predictions and extracellular 
  pcv.v.extra <- merge(unique(growth.preds %>% select(-pcv)), extra.umol.dt)
  
  # make model to predict extracellular aa abundance using predicted pcv
  pve.fits <- pcv.v.extra %>% group_by(compound) %>% do(fit=lm(umol ~ predicted, data=.), data=(.))
  pve.tidy <- pve.fits %>% group_by(compound) %>% do(tidy(.$fit[[1]]))
  pve.tidy$term[pve.tidy$term == "(Intercept)"] = "intercept"
  pve.tidy$term[pve.tidy$term == "predicted"] = "slope"
  pve.spread <- pve.tidy %>% select(-(stderror:p.value)) %>% spread(term, estimate)
  
  growth.k <- growth.curve.params$rate
  uptake.rates <- pve.spread %>% mutate(uptake.rate = -growth.k*slope) %>% select(compound, uptake.rate)
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

# source("ProtScavFluxCalculation/protein_influx_aa_tables.R")
# uptake_uptakeAA <- calc.uptake.rates.unlabeled.only(aa_table_uptake, extra.uptake, growth.list$preds, growth.list$curve.params, num_mLs)

# global variables (no need to make these arguments): unlab_aa_table, alpha.df
# from protein_influx_aa_tables
compute.flux.forUptakeAA <- function (unlab.AA.uptake, known.protein.fluxes) {
  fluxes_uptakeAA <- unlab_aa_table
  fluxes_uptakeAA$umol.protein.flux <- mean(known.protein.fluxes$umol.protein.flux)
  
  merge1 <- merge(fluxes_uptakeAA, alpha.df, by="compound")
  merge1 <- merge1 %>% mutate(resulting.aa.flux = umol.protein.flux * alpha)
  merge2 <- merge(merge1, uptake_uptakeAA$uptake.rates, by="compound")
  
  merge2
}

combine.fluxes <- function (fluxes.labAA, fluxes.unlabAA) {
  rbind(fluxes.labAA %>% select(compound, umol.protein.flux, resulting.aa.flux, uptake.rate),
        fluxes.unlabAA %>% select(compound, umol.protein.flux, resulting.aa.flux, uptake.rate))
}

# foo <- fluxes.list$fluxes
# bar <- compute.flux.forUptakeAA(uptake_uptakeAA, fluxes.list$fluxes)
# foobar <- combine.fluxes(foo, bar)

# load common accessory functions
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/common.accessory.functions.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/aa.info.R")
# aa.info.R generates default.aa.table ; need to use generate.aa.table() if default aa concentrations are not used
# aa.info.R also contains alpha.df, used as a global variable

# integrate product of 2 exponentials
do.integration.product.exponentials <- function (intercept1, slope1, intercept2, slope2, max) {
  intercept1 <- as.numeric(intercept1)
  slope1 <- as.numeric(slope1)
  intercept2 <- as.numeric(intercept2)
  slope2 <- as.numeric(slope2)
  max <- as.numeric(max)
  as.numeric(integrate(function (t) exp((intercept1 + slope1*t) + (intercept2 + slope2*t)), 0, max)[1])
}

# plug calculated concentrations / rates / integrals into final MFA equation
do.prot.scav.flux.equation <- function (alpha, final.unlab.med, uptake.rate, growth.int, prod.int) {
  alpha <- as.numeric(alpha)
  final.unlab.med <- as.numeric(final.unlab.med)
  uptake.rate <- as.numeric(uptake.rate)
  growth.int <- as.numeric(growth.int)
  prod.int <- as.numeric(prod.int)
  numerator <- final.unlab.med + uptake.rate*prod.int
  denominator <- growth.int - prod.int
  (1/alpha)*numerator/denominator
}

# calculate integral of product of exponentials describing growth rate and amino acid labeling
integrate.prod.exponentials.new <- function (growth.params, aa.params) {
  growth.max <- as.numeric(growth.params$max)
  growth.dt <- growth.params %>% select(intercept, rate) %>% mutate(max = growth.max)
  intra.params <- aa.params %>% filter(loc == "Intracellular")
  intra.params <- intra.params %>% mutate(aa.intercept = intercept, aa.slope = rate) %>% 
    select(compound, aa.intercept, aa.slope)
  integrals <- merge(intra.params, growth.dt)
  integrals$integral.prods <- apply(integrals, 1, function (x) {
    do.integration.product.exponentials(x['aa.intercept'], x['aa.slope'], x['intercept'], x['rate'], x['max']) } )
  integrals <- integrals %>% select(compound, integral.prods)
  
  return(integrals)
}

# 1 - calculate flux using MFA equation
# 2 - plot flux estimates from different aas
# 3 - plot influx of amino acids through monomeric uptake vs through protein catabolism 
compute.flux.new <- function (final.umol.unlab, uptake.rates, growth.params, aa.params) {
  output.list <- list()
  
  integral.prods <- integrate.prod.exponentials.new(growth.params, aa.params)
  
  eqn.vars <- merge(alpha.df, merge(final.umol.unlab, uptake.rates))
  eqn.vars$growth.integral <- growth.params$growth.integral
  eqn.vars <- merge(eqn.vars, integral.prods)
  eqn.vars$umol.protein.flux <- apply(eqn.vars, 1, function (x) {
    do.prot.scav.flux.equation(x['alpha'], x['net.umol.unlab'], x['uptake.rate'], x['growth.integral'], x['integral.prods']) } )
  eqn.vars <- eqn.vars %>% mutate(resulting.aa.flux = alpha*umol.protein.flux)
  output.list$fluxes <- eqn.vars
  
  umol.uptake.plot <- ggplot(eqn.vars, aes(x=compound, y=umol.protein.flux*10^6, fill=compound)) + 
    geom_bar(stat="identity", position="dodge", color="black") + scale_fill_hue(c=60, l=70) +
    geom_abline(slope=0) +
    xlab("") + ylab("Flux (pmol / uL cell / hr)")
  output.list$plot.fluxes <- umol.uptake.plot
  
  compare.uptakes <- eqn.vars %>% select(-net.umol.unlab, -growth.integral, -integral.prods) %>%
    gather(uptake.route, flux, c(uptake.rate, resulting.aa.flux))
  compare.uptake.plot <- ggplot(compare.uptakes, aes(x=compound, y=flux, fill=uptake.route)) + 
    geom_bar(stat="identity", position="dodge", color="black") + 
    scale_fill_manual(values=c("grey50","tomato2")) + 
    geom_abline(slope=0) +
    xlab("") + ylab("Flux (umol / uL cell / hr)")
  output.list$plot.compare <- compare.uptake.plot
  
  return(output.list)
}
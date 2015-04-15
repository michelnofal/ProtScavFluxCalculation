# do i need all of these??

library(data.table)
library(reshape2)
library(ggplot2)
library(dplyr)
library(broom)
library(tidyr)
options(stringsAsFactors = FALSE)

setwd("/Users/Michel/Desktop/Research/code/")
source("mn_common/libraries_and_themes_jan15.R")
theme_set(figure_theme)
source("ProtScavFluxCalculation/protein_influx_functions_v3.R")

file_handle <- "Ras_timecourse1"
data_folder <- "ProtScavFluxCalculation/Data/BMK_White_Lab/"

# import data
growth.dat <- read.table("ProtScavFluxCalculation/Data/BMK_White_Lab/BMK_Ras/Ras_pcv_data.csv", header=TRUE, sep=",")
intra.dat <- read.table("ProtScavFluxCalculation/Data/BMK_White_Lab/BMK_Ras/Ras_IC_data.csv", header=TRUE, sep=",")
extra.dat <- read.table("ProtScavFluxCalculation/Data/BMK_White_Lab/BMK_Ras/Ras_media_data.csv", header=TRUE, sep=",")
extra.dat$time.pt[extra.dat$time.pt == "fresh"] <- 0
extra.dat$time.pt <- as.numeric(extra.dat$time.pt)

num_mLs <- 2
media_perc <- 40

# 1 - calculate growth predictions
# 2 - calculate fitted growth curve parameters
# 3 - calculate area under growth curve
# 4 - generate growth plot
growth.list <- growth.func(growth.dat)

# 1 - return organized intra- and extracellular data
# 2 - calculate fitted amino acid labeling curve parameters and predictions
# 3 - generate amino acid labeling plots
aa.lab.list <- aa.lab.func(intra.dat, extra.dat)

# 1 - calculate amino acid uptake rates
# 2 - generate plot of predicted pcv vs amino acid levels in medium (slope = uptake rate / k_growth)
aa_table <- data.table(compound = c("lysine","phenylalanine","threonine","tyrosine","valine"), 
                       DMEM.conc = c(800,400,800,500,800)/1.06*(media_perc/100))
uptake.list <- calc.uptake.rates(aa_table, aa.lab.list$data, growth.list$preds, growth.list$curve.params, num_mLs)

# calculate umol unlab amino acid at final time point
final.umol.unlab <- get.final.umol.unlab(aa_table, aa.lab.list$data, num_mLs)

# calculate integral of product of exponentials describing growth rate and amino acid labeling
integral.prod.curves <- integrate.prod.exponentials(growth.list$curve.params, growth.list$integrals$max, aa.lab.list$curve.params) 

# 1 - calculate flux using MFA equation
# 2 - plot flux estimates from different aas
# 3 - plot influx of amino acids through monomeric uptake vs through protein catabolism 
alpha.df <- data.frame(compound = c("lysine","phenylalanine","threonine","tyrosine","valine"), alpha = c(59, 27, 33, 20, 36))
fluxes.list <- compute.flux(alpha.df, final.umol.unlab, uptake.list$uptake.rates, growth.list$integrals, integral.prod.curves)

showQuickPlots = FALSE
saveQuickPlots = FALSE
saveRData = FALSE

# showQuickPlots = TRUE
if (showQuickPlots) {
  growth.list$plot
  aa.lab.list$plots
  uptake.list$plot
  fluxes.list$plot.fluxes
  fluxes.list$plot.compare
}

# saveQuickPlots = TRUE
if (saveQuickPlots) {
  save.plots(growth.list$plot, aa.lab.list$plots, uptake.list$plot, fluxes.list$plot.fluxes, fluxes.list$plot.compare, 
             directory="ProtScavFluxCalculation/Figures/Quick_plots/", file_handle=file_handle)
}

# saveRData = TRUE
if (saveRData) {
  Ras1_fluxes <- fluxes.list$fluxes
  Ras1_fluxes$data = "BMK_Ras_1"
  save(Ras1_fluxes, file=paste(data_folder,"R_data/Ras1_fluxes.Rda",sep=""))
}


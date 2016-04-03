options(stringsAsFactors = FALSE)
library(ggplot2)
library(dplyr)
library(broom)
library(tidyr)
library(data.table)

source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/aa.info.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/common.accessory.functions.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/PCV.curvefit.plot.and.compare.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/aa.lab.func.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/calc.extracellularAA.mol.and.uptake.R")
source("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/data.processing.functions/compute.flux.R")

calculate.protein.scavenging.flux <- function (pcv.data, intra.data, media.data, data.name, 
                                               num.mLs=2, aa.table=default.aa.table, adjust=FALSE) {
  output <- list()
  
  # process pcv data
  pcv.list <- growth.func.new(pcv.data, data.name)
  output$pcv <- pcv.list
  
  # process aa labeling data
  aa.lab.list <- aa.lab.func(intra.data, media.data, data.name, adjust=adjust)
  # create data frame containing only extracellular aa data (to use below)
  media.aa.data <- aa.lab.list$data %>% filter(loc == "Extracellular")
  output$aa <- aa.lab.list
  
  # convert media aa data from counts to umol
  extra.umol.data <- get.extra.umol.new(media.aa.data, num.mLs=num.mLs)
  # calculate umols of unlab amino acid at final time point
  final.umol.unlab <- get.final.umol.unlab.new(extra.umol.data, media.aa.data)
  # calculate uptake rates for all aa
  uptake.list <- calc.uptake.rates.new(extra.umol.data, pcv.list, data.name)
  output$uptake <- uptake.list
  
  # compute protein scavenging flux
  fluxes.list <- compute.flux.new(final.umol.unlab, uptake.rates=uptake.list$uptake.rates, 
                                  growth.params=pcv.list$params,  aa.params=aa.lab.list$curve.params)
  output$fluxes <- fluxes.list
  
  return(output)
}

###
### Example
### 

# growth.dat <- read.table("ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/krpc_HF_lab_timecourse_pcv.csv", header=TRUE, sep=",")
# intra.dat <- read.table("ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/Intracellular_data_HFB_processed.csv", header=TRUE, sep=",")
# extra.dat <- read.table("ProtScavFluxCalculation/Data/KRPCA_May2015_timecourses/Extracellular_data_HFB_processed.csv", header=TRUE, sep=",")
# num_mLs <- 2
# media_perc <- 40

# extra.dat <- extra.dat %>% filter(time.pt != "0")
# extra.dat <- extra.dat %>% mutate(time.pt = ifelse(time.pt == "fresh", 0, time.pt))
# extra.dat$time.pt = as.numeric(extra.dat$time.pt)

# extra.dat <- extra.dat %>% filter(compound %in% c("histidine","lysine","phenylalanine","threonine","tyrosine","valine"))

# test <- calculate.protein.scavenging.flux(growth.dat, intra.dat, extra.dat, "HF_BLAHBLAH")
# test.fluxes <- test$fluxes

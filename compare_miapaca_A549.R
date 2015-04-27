theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/MFA/code/")

source("ProtScavFluxCalculation/comparison_functions.R")

load(file="ProtScavFluxCalculation/Data/Misc_Human_KRAS/R_data/A549_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/Misc_Human_KRAS/R_data/MiaPaCa_fluxes.Rda")

A549_Mia_fluxes <- rbind(A549_fluxes, MiaPaCa_fluxes)

AAtoplot <- c("lysine","phenylalanine","threonine")
A549_Mia_labels <- c("A549", "MiaPaCa2")

ProtScav.Comparison.plotter(A549_Mia_fluxes, AAs = c(AAtoplot), expt_labels = A549_Mia_labels, save=TRUE, 
                            directory="ProtScavFluxCalculation/Figures/Misc_Human_KRAS/", filename="A549_MiaPaCa")

Compare_Timecourses_and_Fluxes(A549_Mia_fluxes, AAs = c(AAtoplot), expt_labels=A549_Mia_labels, normTo1 = FALSE,
                               save=TRUE, directory="ProtScavFluxCalculation/Figures/Misc_Human_KRAS/", filename="A549_MiaPaCa")

A_M_KRPC_fluxes <- rbind()
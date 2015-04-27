theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/MFA/code/")

source("ProtScavFluxCalculation/comparison_functions.R")

load(file="ProtScavFluxCalculation/Data/Misc_Human_KRAS/R_data/A549_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/Misc_Human_KRAS/R_data/MiaPaCa_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/KRPCA_dropout_comparison/R_data/KRPCA_DB_fluxes.Rda")

AAtoplot <- c("lysine","phenylalanine","threonine")

A_M_KRPC_fluxes <- rbind(KRPCA_DB_fluxes, A549_fluxes, MiaPaCa_fluxes)

A_M_KRPC_fluxes$data <- factor(A_M_KRPC_fluxes$data, levels=c("KRPCA_DB_1","MiaPaCa","A549"))
A_M_KRPC_labels <- c(expression(paste("KRPC"^"A")), "MiaPaCa2", "A549")

ProtScav.Comparison.plotter(A_M_KRPC_fluxes, AAs = c(AAtoplot), expt_labels = A_M_KRPC_labels, save=TRUE, 
                            directory="ProtScavFluxCalculation/Figures/Misc_Human_KRAS/", filename="KRPCA_Mia_A549")

Compare_Timecourses_and_Fluxes(A_M_KRPC_fluxes, AAs = c(AAtoplot), expt_labels=A_M_KRPC_labels, normTo1 = FALSE,
                               save=TRUE, directory="ProtScavFluxCalculation/Figures/Misc_Human_KRAS/", filename="KRPCA_Mia_A549")

Compare_Timecourses_and_Fluxes(A_M_KRPC_fluxes, AAs= c(AAtoplot, "tyrosine", "valine"), expt_labels=A_M_KRPC_labels, normTo1 = FALSE,
                               save=FALSE, directory="ProtScavFluxCalculation/Figures/Misc_Human_KRAS/", filename="KRPCA_Mia_A549")


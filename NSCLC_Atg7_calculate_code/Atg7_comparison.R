theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/MFA/code/")

source("ProtScavFluxCalculation/comparison_functions.R")

load(file="ProtScavFluxCalculation/Data/NSCLC_Atg_WTvNull/R_data/Atg7_28717_wBSA_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/NSCLC_Atg_WTvNull/R_data/Atg7_28717_noBSA_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/NSCLC_Atg_WTvNull/R_data/Atg7_28298_slow_wBSA_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/NSCLC_Atg_WTvNull/R_data/Atg7_28298_slow_noBSA_fluxes.Rda")

AAtoplot <- c("lysine","phenylalanine","threonine","tyrosine","valine")

Atg7_1_fluxes <- rbind(Atg7_28717_wBSA_fluxes, Atg7_28717_noBSA_fluxes, Atg7_28298_slow_wBSA_fluxes, Atg7_28298_slow_noBSA_fluxes)

Atg7_1_fluxes$data <- factor(Atg7_1_fluxes$data, levels=c("Atg7_28717_noBSA","Atg7_28717_wBSA","28298_slow_noBSA","28298_slow_wBSA"))
Atg7_1_fluxes_labels <- c("Atg7+/+ (no BSA)","Atg7+/+","Atg7-/- (no BSA)", "Atg7-/-")

ProtScav.Comparison.plotter(Atg7_1_fluxes, AAs = c(AAtoplot), expt_labels = Atg7_1_fluxes_labels, save=TRUE, 
                            directory="ProtScavFluxCalculation/Figures/NSCLC_Atg7/", filename="Atg7_comparison")

# Compare_Timecourses_and_Fluxes(, AAs = c(AAtoplot), expt_labels=HOPS_fluxes_labels, normTo1 = FALSE,
#                                save=TRUE, directory="ProtScavFluxCalculation/Figures/HOPS complex/", filename="HOPS_comparison")

# Compare_Timecourses_and_Fluxes(HOPS_fluxes, AAs= c(AAtoplot, "tyrosine", "valine"), expt_labels=HOPS_fluxes_labels, normTo1 = FALSE,
#                                save=FALSE, directory="ProtScavFluxCalculation/Figures/HOPS complex/", filename="HOPS_comparison")

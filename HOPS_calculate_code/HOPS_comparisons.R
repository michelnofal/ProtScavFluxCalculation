theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/code/")

source("ProtScavFluxCalculation/comparison_functions.R")
source("ProtScavFluxCalculation/protein_influx_clean_plots.R")

load(file="ProtScavFluxCalculation/Data/KRPCA_dropout_comparison/R_data/KRPCA_DB_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/HOPS complex/R_data/HOPS1_intergenic_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/HOPS complex/R_data/HOPS1_Vps39_1_1C_fluxes.Rda")
load(file="ProtScavFluxCalculation/Data/HOPS complex/R_data/HOPS1_Vps39_8_5C_fluxes.Rda")

AAtoplot <- c("lysine","phenylalanine","threonine","valine")

HOPS_fluxes <- rbind(KRPCA_DB_fluxes, HOPS1_intergenic_fluxes, HOPS1_Vps39_1_1C_fluxes, HOPS1_Vps39_8_5C_fluxes)

HOPS_fluxes$data <- factor(HOPS_fluxes$data, levels=c("KRPCA_DB_1","HOPS1_intergenic","HOPS1_Vps39_1_1C","HOPS1_Vps39_8_5C"))
HOPS_fluxes_labels <- c("Uninfected","Intergenic control","sgVps39_1_1C","sgVps39_8_5C")

ProtScav.Comparison.plotter(HOPS_fluxes, AAs = c(AAtoplot), expt_labels = HOPS_fluxes_labels, save=TRUE, 
                            directory="ProtScavFluxCalculation/Figures/HOPS complex/", filename="HOPS_comparison")

Compare_Timecourses_and_Fluxes(HOPS_fluxes, AAs = c(AAtoplot), expt_labels=HOPS_fluxes_labels, normTo1 = FALSE,
                               save=TRUE, directory="ProtScavFluxCalculation/Figures/HOPS complex/", filename="HOPS_comparison")

HOPS_summ <- HOPS_fluxes %>% filter(compound != "tyrosine") %>% group_by(data) %>% summarize(median = median(umol.protein.flux))


ProtScav.Comparison.plotter(HOPS_fluxes, AAs = c(AAtoplot), expt_labels = HOPS_fluxes_labels, save=FALSE)
ggplot(HOPS_summ, aes(x=data,y=median)) + geom_bar(stat="identity")



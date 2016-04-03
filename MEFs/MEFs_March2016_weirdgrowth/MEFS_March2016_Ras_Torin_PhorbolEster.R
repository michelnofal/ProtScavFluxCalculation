options(stringsAsFactors = FALSE)
setwd("/Users/Michel/Desktop/Research/code/")

source("mn_common/libraries_and_themes_jan15.R")
theme_set(figure_theme)
source("ProtScavFluxCalculation/data.processing.functions/calculate.protein.scavenging.flux.R")

source("ProtScavFluxCalculation/comparison_functions.R")
source("ProtScavFluxCalculation/protein_influx_clean_plots_v2.R") # just for scav.plot.colors

# importing growth data
kmef.growth.dat <- read.table("ProtScavFluxCalculation/MEFs/MEFs_March2016_weirdgrowth/MEFs_growth.csv", header=TRUE, sep=",")

# any effect of any drug on cell size?
kmef.cellsize <- kmef.growth.dat %>% filter(time.pt == 24) %>% mutate(pL = pcv/((cc1*1e6 + cc2*1e6)/2)*1e6)
kmef.cellsize.summ <- kmef.cellsize %>% group_by(cell.line, condition) %>% summarize(mean = mean(pL), stderr = sqrt(var(pL)/length(pL)))

kmef.cellsize.plot <- ggplot(kmef.cellsize.summ, aes(x=condition, y=mean)) + 
  geom_bar(stat="identity", color="black", fill="grey") + geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin=mean-1.96*stderr, ymax=mean+1.96*stderr), width=0.5) +
  xlab("Condition") + ylab("Volume per cell (pL)") + facet_wrap(~ cell.line)
kmef.cellsize.plot
# maybe... messy data

# processing pcv data
kmef.growth.clean <- kmef.growth.dat %>% dplyr::select(-cc1, -cc2)

wt.notrx.pcv <- kmef.growth.clean %>% filter(cell.line == "WT MEFs", condition %in% c("common","noTrx"))
wt.torin.pcv <- kmef.growth.clean %>% filter(cell.line == "WT MEFs", condition %in% c("common","Torin1"))
wt.pma.pcv <- kmef.growth.clean %>% filter(cell.line == "WT MEFs", condition %in% c("common","PMA"))
wt.pdbu.pcv <- kmef.growth.clean %>% filter(cell.line == "WT MEFs", condition %in% c("common","PDBu"))
kras.notrx.pcv <- kmef.growth.clean %>% filter(cell.line == "KRas MEFs", condition %in% c("common","noTrx"))
kras.torin.pcv <- kmef.growth.clean %>% filter(cell.line == "KRas MEFs", condition %in% c("common","Torin1"))
kras.pma.pcv <- kmef.growth.clean %>% filter(cell.line == "KRas MEFs", condition %in% c("common","PMA"))
kras.pdbu.pcv <- kmef.growth.clean %>% filter(cell.line == "KRas MEFs", condition %in% c("common","PDBu"))

wt.notrx.proc <- growth.func.new(wt.notrx.pcv, name = "WT noTrx")
wt.torin.proc <- growth.func.new(wt.torin.pcv, name = "WT Torin")
wt.pma.proc <- growth.func.new(wt.pma.pcv, name = "WT PMA")
wt.pdbu.proc <- growth.func.new(wt.pdbu.pcv, name = "WT PDBu")
kras.notrx.proc <- growth.func.new(kras.notrx.pcv, name = "KRas noTrx")
kras.torin.proc <- growth.func.new(kras.torin.pcv, name = "KRas Torin")
kras.pma.proc <- growth.func.new(kras.pma.pcv, name = "KRas PMA")
kras.pdbu.proc <- growth.func.new(kras.pdbu.pcv, name = "KRas PDBu")

wtmef.growth.comp <- growth.comparer.new(wt.notrx.proc, wt.torin.proc, wt.pma.proc, wt.pdbu.proc)
wtmef.growth.comp$plot + theme(legend.position = "right")
wtmef.growth.comp$meta
krasmef.growth.comp <- growth.comparer.new(kras.notrx.proc, kras.torin.proc, kras.pma.proc, kras.pdbu.proc)
krasmef.growth.comp$plot + theme(legend.position = "right")
krasmef.growth.comp$meta

# importing ms data
kmef.ms.dat <- tbl_df(read.table("ProtScavFluxCalculation/MEFs/MEFs_March2016_weirdgrowth/KRas MEFs Peaks.csv", header=TRUE, sep=","))
kmef.ms.dat <- kmef.ms.dat %>% dplyr::select(-X, -sample)

kmef.ms.dat %>% filter(compound == "lysine", !(time.pt %in% c(16,24)))

unique(kmef.ms.dat$cell.line)
unique(kmef.ms.dat$drug)

lab.aa <- c("lysine", "phenylalanine","threonine","tyrosine","valine","histidine")

kmef.ms.gather <- kmef.ms.dat %>% gather(compound, ion.count, -(sample:prelabeled))
kmef.ms.gather <- kmef.ms.gather %>% group_by(compound) %>% mutate(isotope = strsplit(as.character(compound), "_")[[1]][1],
                                                                   amino.acid = strsplit(as.character(compound), "_")[[1]][2])
kmef.ms.spread <- kmef.ms.gather %>% ungroup %>% select(-compound) %>% spread(isotope, ion.count, fill=0)
kmef.ms.spread <- kmef.ms.spread %>% mutate(compound = amino.acid) %>% select(-lab1, -lab5, -lab6, -unlab, -amino.acid)
kmef.ms.5aa <- kmef.ms.spread %>% filter(compound %in% c("histidine","lysine","phenylalanine","threonine","valine"))
kmef.ms.5aa$time.pt <- as.numeric(kmef.ms.5aa$time.pt)
kmef.ms.5aa %>% filter(compound == "threonine", condition %in% c("noDFBS_common","noDFBS"))

kmef.one.noDFBS.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","one"), condition %in% c("noDFBS","noDFBS_common"))
kmef.one.noDT.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","one"), condition %in% c("noDFBS+Torin","noDFBS_common"))
kmef.one.DFBS.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","one"), condition %in% c("DFBS","DFBS_common"))
kmef.one.DT.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","one"), condition %in% c("DFBS+Torin","DFBS_common"))
kmef.five.noDFBS.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","five"), condition %in% c("noDFBS","noDFBS_common"))
kmef.five.noDT.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","five"), condition %in% c("noDFBS+Torin","noDFBS_common"))
kmef.five.DFBS.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","five"), condition %in% c("DFBS","DFBS_common"))
kmef.five.DT.ms <- kmef.ms.5aa %>% filter(prelabeled %in% c("common","five"), condition %in% c("DFBS+Torin","DFBS_common"))
kmef.one.noDFBS.ms %>% filter(compound == "threonine")

noDFBS.aatable <- generate.aa.table(perc.fbs.added = 0)

kmef.one.noDFBS.scav <- calculate.protein.scavenging.flux(kmef.one.noDFBS.pcv, kmef.one.noDFBS.ms, kmef.one.noDFBS.ms, 
                                                          "kmef.one.noDFBS", aa.table = noDFBS.aatable)
kmef.one.noDT.scav <- calculate.protein.scavenging.flux(kmef.one.noDT.pcv, kmef.one.noDT.ms, kmef.one.noDT.ms, 
                                                        "kmef.one.noDT", aa.table = noDFBS.aatable)
kmef.one.DFBS.scav <- calculate.protein.scavenging.flux(kmef.one.DFBS.pcv, kmef.one.DFBS.ms, kmef.one.DFBS.ms, "kmef.one.DFBS")
kmef.one.DT.scav <- calculate.protein.scavenging.flux(kmef.one.DT.pcv, kmef.one.DT.ms, kmef.one.DT.ms, "kmef.one.DT")

kmef.five.noDFBS.scav <- calculate.protein.scavenging.flux(kmef.five.noDFBS.pcv, kmef.five.noDFBS.ms, kmef.five.noDFBS.ms, 
                                                           "kmef.five.noDFBS", aa.table = noDFBS.aatable)
kmef.five.noDT.scav <- calculate.protein.scavenging.flux(kmef.five.noDT.pcv, kmef.five.noDT.ms, kmef.five.noDT.ms, 
                                                         "kmef.five.noDT", aa.table = noDFBS.aatable)
kmef.five.DFBS.scav <- calculate.protein.scavenging.flux(kmef.five.DFBS.pcv, kmef.five.DFBS.ms, kmef.five.DFBS.ms, "kmef.five.DFBS")
kmef.five.DT.scav <- calculate.protein.scavenging.flux(kmef.five.DT.pcv, kmef.five.DT.ms, kmef.five.DT.ms, "kmef.five.DT")

kmef.compare <- rbind(kmef.one.noDFBS.scav$fluxes$fluxes, kmef.one.noDT.scav$fluxes$fluxes, 
                      kmef.one.DFBS.scav$fluxes$fluxes, kmef.one.DT.scav$fluxes$fluxes, 
                      kmef.five.noDFBS.scav$fluxes$fluxes, kmef.five.noDT.scav$fluxes$fluxes, 
                      kmef.five.DFBS.scav$fluxes$fluxes, kmef.five.DT.scav$fluxes$fluxes)
kmef.compare$data <- factor(kmef.compare$data, levels=c("kmef.one.noDFBS","kmef.one.noDT","kmef.one.DFBS","kmef.one.DT",
                                                        "kmef.five.noDFBS","kmef.five.noDT","kmef.five.DFBS","kmef.five.DT")) 
kmef.compare.labels <- c("kmef.one.noDFBS","kmef.one.noDT","kmef.one.DFBS","kmef.one.DT",
                         "kmef.five.noDFBS","kmef.five.noDT","kmef.five.DFBS","kmef.five.DT")

lab.aa <- c("histidine","phenylalanine","threonine","valine")
ProtScav.Comparison.plotter(kmef.compare, AAs = lab.aa, expt_labels = kmef.compare.labels, save=T, plot_width=15, plot_height=12, 
                            directory="ProtScavFluxCalculation/Figures/KRas_MEF/", filename="KMef_Absolute_fluxes")

Compare_Timecourses_Stagger(kmef.compare, AAs = lab.aa, expt_labels = kmef.compare.labels, save=F, normTo1 = F,
                            plot_width=10, plot_height=10, 
                            directory="ProtScavFluxCalculation/Figures/KRas_MEF/", filename="KMef")

abs.uptake.fluxes <- kmef.compare %>% mutate(total.uptake = uptake.rate + resulting.aa.flux)
total.uptake.flux.plot <- ggplot(abs.uptake.fluxes, aes(x=data, y=total.uptake)) + geom_abline(slope=0) +
  geom_bar(stat="identity", position="dodge",color="black",fill="tomato") + facet_wrap( ~ compound) + 
  xlab("") + ylab("Total uptake flux of AA") + theme(axis.text.x = element_text(angle=90))
total.uptake.flux.plot


ProtScav.Comparison.plotter(kmef.compare, AAs = lab.aa, expt_labels = kmef.compare.labels, save=F, plot_width=15, plot_height=12, 
                            directory="ProtScavFluxCalculation/Figures/KRas_MEF/", filename="KMef_Absolute_fluxes")

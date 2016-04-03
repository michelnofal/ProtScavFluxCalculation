options(stringsAsFactors = FALSE)
setwd("/Users/Michel/Desktop/Research/code/")

source("mn_common/libraries_and_themes_jan15.R")
theme_set(figure_theme)
source("ProtScavFluxCalculation/data.processing.functions/calculate.protein.scavenging.flux.R")

source("ProtScavFluxCalculation/comparison_functions.R")
source("ProtScavFluxCalculation/protein_influx_clean_plots_v2.R") # just for scav.plot.colors

# importing growth data
kmef.growth.dat <- read.table("ProtScavFluxCalculation/Data/KRas_MEF/KRasMefs_growthdat.csv", header=TRUE, sep=",")

kmef.growth.gather <- kmef.growth.dat %>% gather(measurement, value, -prelabeled, -condition, -time.pt)
kmef.growth.gather <- kmef.growth.gather %>% filter(!is.na(value))
kmef.growth.split <- kmef.growth.gather %>% group_by(measurement) %>% mutate(meas = strsplit(as.character(measurement), "_")[[1]][1],
                                                                             rep = strsplit(as.character(measurement), "_")[[1]][2])
kmef.pcv <- kmef.growth.split %>% filter(meas == "PCV") %>% mutate(pcv = value) %>% 
  ungroup %>% dplyr::select(-measurement, -value, -meas, -rep)

kmef.cc <- kmef.growth.split %>% filter(meas == "cellnum") %>% group_by(rep) %>% mutate(num = strsplit(as.character(rep), "\\.")[[1]][1],
                                                                                        lett = strsplit(as.character(rep), "\\.")[[1]][2])
kmef.cellnum <- kmef.cc %>% group_by(prelabeled, condition, time.pt, num) %>% summarize(cellnum = mean(value)) %>% dplyr::select(-num)

# krpca.pcv <- nov2015.pcv.dat %>% filter(cell.line == "krpca") %>% select(-cell.line)
# krpcu.pcv <- nov2015.pcv.dat %>% filter(cell.line == "krpcu") %>% select(-cell.line)
# atg7_2.pcv <- nov2015.pcv.dat %>% filter(cell.line == "atg7_2") %>% select(-cell.line)

kmef.one.noDFBS.pcv <- kmef.pcv %>% filter(prelabeled == "one", condition %in% c("common", "noDFBS"))
kmef.one.noDFBS.proc <- growth.func.new(kmef.one.noDFBS.pcv, name = "kmef.one.noDFBS")
kmef.one.noDT.pcv <- kmef.pcv %>% filter(prelabeled == "one", condition %in% c("common", "noDFBS+Torin"))
kmef.one.noDT.proc <- growth.func.new(kmef.one.noDT.pcv, name = "kmef.one.noDFBS+Torin")
kmef.one.DFBS.pcv <- kmef.pcv %>% filter(prelabeled == "one", condition %in% c("common", "DFBS"))
kmef.one.DFBS.proc <- growth.func.new(kmef.one.DFBS.pcv, name = "kmef.one.DFBS")
kmef.one.DT.pcv <- kmef.pcv %>% filter(prelabeled == "one", condition %in% c("common", "DFBS+Torin"))
kmef.one.DT.proc <- growth.func.new(kmef.one.DT.pcv, name = "kmef.one.DFBS+Torin")

# growth.curve.plotter(kmef.one.noDT.proc) 
kmef.one.growth.comp <- growth.comparer.new(kmef.one.noDFBS.proc, kmef.one.noDT.proc, kmef.one.DFBS.proc, kmef.one.DT.proc)
kmef.one.growth.comp$plot
kmef.one.growth.comp$meta

kmef.five.noDFBS.pcv <- kmef.pcv %>% filter(prelabeled == "five", condition %in% c("common", "noDFBS"))
kmef.five.noDFBS.proc <- growth.func.new(kmef.five.noDFBS.pcv, name = "kmef.five.noDFBS")
kmef.five.noDT.pcv <- kmef.pcv %>% filter(prelabeled == "five", condition %in% c("common", "noDFBS+Torin"))
kmef.five.noDT.proc <- growth.func.new(kmef.five.noDT.pcv, name = "kmef.five.noDFBS+Torin")
kmef.five.DFBS.pcv <- kmef.pcv %>% filter(prelabeled == "five", condition %in% c("common", "DFBS"))
kmef.five.DFBS.proc <- growth.func.new(kmef.five.DFBS.pcv, name = "kmef.five.DFBS")
kmef.five.DT.pcv <- kmef.pcv %>% filter(prelabeled == "five", condition %in% c("common", "DFBS+Torin"))
kmef.five.DT.proc <- growth.func.new(kmef.five.DT.pcv, name = "kmef.five.DFBS+Torin")

# growth.curve.plotter(kmef.five.noDT.proc) 
kmef.five.growth.comp <- growth.comparer.new(kmef.five.noDFBS.proc, kmef.five.noDT.proc, kmef.five.DFBS.proc, kmef.five.DT.proc)
kmef.five.growth.comp$plot
kmef.five.growth.comp$meta
kmef.all.growth.meta <- rbind(kmef.one.growth.comp$meta, kmef.five.growth.comp$meta)

# importing ms data
kmef.ms.dat <- tbl_df(read.table("ProtScavFluxCalculation/Data/KRas_MEF/KRasMefs4_toImport.csv", header=TRUE, sep=","))

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

setwd("/Users/Michel/Desktop/Research/CRISPR screen/")
theme_set(figure_theme)
setwd("/Users/Michel/Desktop/Research/MFA/code/")

# code to run for examples
# load(file="Akt1_fluxes.Rda")
# load(file="D3_1_fluxes.Rda")
# load(file="Ras1_fluxes.Rda")

# AAtoplot <- c("lysine","phenylalanine","threonine")

# D3_Akt_Ras_fluxes <- rbind(D3_1_fluxes, Akt1_fluxes, Ras1_fluxes)
# D3_Akt_Ras_fluxes$data <- factor(D3_Akt_Ras_fluxes$data, levels=c("BMK_D3_1","BMK_Ras_1","BMK_Akt_1"))
# BMK_labels <- c("Parental", "myrAkt", expression(paste("HRas"^"G12V")))

ProtScav.Comparison.plotter <- function(multi_expt_fluxes, AAs, expt_labels, save=FALSE, 
                                        directory="ProtScavFluxCalculation/Figures/", filename="choose_filename") {
  col.df <- scav.plot.colors %>% filter(aa %in% fluxes.toPlot$compound)
  
  multi_expt_fluxes <- multi_expt_fluxes %>% filter(compound %in% AAs)
  
  comp.umol.uptakePlot <- ggplot(multi_expt_fluxes, aes(x=data, y=umol.protein.flux*10^6, fill=compound)) + 
    geom_bar(stat="identity", color="black", width=0.7, position = position_dodge(width=0.8)) + # position="dodge", 
    scale_fill_manual(values=col.df$line_colors) +
    scale_x_discrete(labels=expt_labels) + 
    geom_abline(height=0, slope=0) + 
    theme(axis.text.x=element_text(size=8, angle=45,hjust=1),
          axis.title.y=element_text(size=8)) + 
    xlab("") + ylab("Scavenging Flux\n(pmol / uL cell / hr)") # + 
    # facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
  
  if (save) {
    ProtScavFlux.filename <- paste(directory, filename, ".Compare.ProteinScavFlux.pdf", sep="")
    ggsave(comp.umol.uptakePlot, file=ProtScavFlux.filename, width=7, height=6, unit="cm")
  }
  comp.umol.uptakePlot
}

# EXAMPLES:
# ProtScav.Comparison.plotter(D3_Akt_Ras_fluxes, AAs = c(AAtoplot), expt_labels = BMK_labels)
# ProtScav.Comparison.plotter(D3_Akt_Ras_fluxes, AAs = c(AAtoplot), expt_labels = BMK_labels, save=TRUE, 
#                             directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")

Compare_Timecourses_and_Fluxes <- function(multi_expt_fluxes, AAs, expt_labels, normTo1=TRUE, 
                                   save=FALSE, directory="ProtScavFluxCalculation/Figures/", filename="choose_filename") {
  col.df <- scav.plot.colors %>% filter(aa %in% multi_expt_fluxes$compound)
  
  multi_expt_fluxes <- D3_Akt_Ras_fluxes
  multi_expt_fluxes <- multi_expt_fluxes %>% filter(compound %in% AAs)
  
  multi.fluxes.toCompare <- multi_expt_fluxes %>% group_by(data) %>% 
    gather(entryRoute, aa.flux, resulting.aa.flux, uptake.rate) %>% select(compound, data, entryRoute, aa.flux) %>% 
    mutate(comp_entryRoute = paste(compound, entryRoute, sep="_")) %>%
    group_by(compound, data) %>% mutate(norm.flux = aa.flux/sum(aa.flux))
  if (normTo1) {
    compareFlux.plot <- ggplot(multi.fluxes.toCompare, aes(x=data, y=norm.flux, color=compound, fill=comp_entryRoute)) + 
      geom_bar(stat="identity", width=0.8) + geom_abline(height=0, slope=0) + 
      scale_color_manual(values=col.df$line_colors) +
      scale_fill_manual(values=alternate.white.and.fill(col.df$fill_colors)) +
      scale_x_discrete(labels=expt_labels) + 
      theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
      xlab("") + ylab("Fractional\nAmino Acid Uptake") +
      facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
    
    if (save == TRUE) {
      compFlux.filename <- paste(directory, filename, ".Comp.Timecourses.AbsoluteFluxes.pdf", sep="")
      ggsave(compareFlux.plot, file=compFlux.filename, width=7.3, height=6, unit="cm")
    }
  } else {
    compareFlux.plot <- ggplot(multi.fluxes.toCompare, aes(x=data, y=aa.flux*1000, color=compound, fill=comp_entryRoute)) + 
      geom_bar(stat="identity", width=0.8) + geom_abline(height=0, slope=0) + 
      scale_color_manual(values=col.df$line_colors) +
      scale_fill_manual(values=alternate.white.and.fill(col.df$fill_colors)) +
      scale_x_discrete(labels=expt_labels) + 
      theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
      xlab("") + ylab("Amino Acid Uptake\n(nmol / uL cell / hr)") +
      facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
    
    if (save == TRUE) {
      compFlux.filename <- paste(directory, filename, ".Comp.Timecourses.AbsoluteFluxes.pdf", sep="")
      ggsave(compareFlux.plot, file=compFlux.filename, width=7, height=6, unit="cm")
    }
  }
  compareFlux.plot
}

# EXAMPLES:
# Compare_Timecourses_and_Fluxes(fluxes.toPlot, AAs = c(AAtoplot), expt_labels=BMK_labels, normTo1 = FALSE,
#                                save=TRUE, directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")



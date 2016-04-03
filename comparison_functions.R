# code to run for examples
# setwd("/Users/Michel/Desktop/Research/code/ProtScavFluxCalculation/Data/BMK_White_Lab/R_data/")

ProtScav.Comparison.plotter <- function(multi_expt_fluxes, AAs, expt_labels, save=FALSE, 
                                        directory="ProtScavFluxCalculation/Figures/", filename="choose_filename",
                                        plot_width=7, plot_height=6) {
  col.df <- scav.plot.colors %>% filter(amino.acid %in% AAs)
  
  multi_expt_fluxes <- multi_expt_fluxes %>% filter(compound %in% AAs)
  multi_expt_fluxes$compound <- factor(multi_expt_fluxes$compound, 
                                       levels=scav.plot.colors$amino.acid[scav.plot.colors$amino.acid %in% lab.aa])
  
  comp.umol.uptakePlot <- ggplot(multi_expt_fluxes, aes(x=data, y=umol.protein.flux*10^6, fill=compound)) + 
    geom_bar(stat="identity", color="black", width=0.7, position = position_dodge(width=0.8)) + # position="dodge", 
    scale_fill_manual(values=col.df$line_colors) +
    scale_x_discrete(labels=expt_labels) + 
    geom_abline(slope=0) + 
    theme(axis.text.x=element_text(size=8, angle=45,hjust=1),
          axis.title.y=element_text(size=8)) + 
    xlab("") + ylab("Scavenging Flux\n(pmol / uL cell / hr)") # + 
  # facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
  comp.umol.uptakePlot
  
  if (save) {
    ProtScavFlux.filename <- paste(directory, filename, ".Compare.ProteinScavFlux.pdf", sep="")
    ggsave(comp.umol.uptakePlot, file=ProtScavFlux.filename, width=plot_width, height=plot_height, unit="cm")
  }
  comp.umol.uptakePlot
}


# EXAMPLES:
# ProtScav.Comparison.plotter(D3_Akt_Ras_fluxes, AAs = c(AAtoplot), expt_labels = BMK_labels)
# ProtScav.Comparison.plotter(D3_Akt_Ras_fluxes, AAs = c(AAtoplot), expt_labels = BMK_labels, save=TRUE, 
#                             directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")

# need this list to order fill colors of plots
comp_entry_levels <- c("lysine_resulting.aa.flux", "lysine_uptake.rate", 
                       "phenylalanine_resulting.aa.flux", "phenylalanine_uptake.rate", 
                       "threonine_resulting.aa.flux", "threonine_uptake.rate",
                       "tyrosine_resulting.aa.flux", "tyrosine_uptake.rate",
                       "valine_resulting.aa.flux", "valine_uptake.rate",
                       "histidine_resulting.aa.flux", "histidine_uptake.rate")

Compare_Timecourses_and_Fluxes <- function(multi_expt_fluxes, AAs, expt_labels, normTo1=TRUE, 
                                   save=FALSE, directory="ProtScavFluxCalculation/Figures/", filename="choose_filename",
                                   plot_width=7, plot_height=6) {
  col.df <- scav.plot.colors %>% filter(amino.acid %in% AAs)
  
  multi_expt_fluxes <- multi_expt_fluxes %>% filter(compound %in% AAs)
  
  multi.fluxes.toCompare <- multi_expt_fluxes %>% group_by(data) %>% 
    gather(entryRoute, aa.flux, resulting.aa.flux, uptake.rate) %>% select(compound, data, entryRoute, aa.flux) %>% 
    mutate(comp_entryRoute = paste(compound, entryRoute, sep="_")) %>%
    group_by(compound, data) %>% mutate(norm.flux = aa.flux/sum(aa.flux))
  
  multi.fluxes.toCompare$compound <- factor(multi.fluxes.toCompare$compound, 
                                            levels=scav.plot.colors$amino.acid[scav.plot.colors$amino.acid %in% AAs])
  multi.fluxes.toCompare$comp_entryRoute <- factor(multi.fluxes.toCompare$comp_entryRoute, 
                                                   levels=comp_entry_levels[which(comp_entry_levels %in% unique(multi.fluxes.toCompare$comp_entryRoute))])
  
  if (normTo1) {
    compareFlux.plot <- ggplot(multi.fluxes.toCompare, aes(x=data, y=norm.flux, color=compound, fill=comp_entryRoute)) + 
      geom_bar(stat="identity", width=0.8) + geom_abline(slope=0) + 
      scale_color_manual(values=col.df$line_colors) +
      scale_fill_manual(values=alternate.white.and.fill(col.df$fill_colors)) +
      scale_x_discrete(labels=expt_labels) + 
      theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
      xlab("") + ylab("Fractional\nAmino Acid Uptake") +
      facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
    
    if (save == TRUE) {
      compFlux.filename <- paste(directory, filename, ".Comp.Timecourses.AbsoluteFluxes.pdf", sep="")
      ggsave(compareFlux.plot, file=compFlux.filename, width=plot_width, height=plot_height, unit="cm")
    }
  } else {
    compareFlux.plot <- ggplot(multi.fluxes.toCompare, aes(x=data, y=aa.flux*1000, color=compound, fill=comp_entryRoute)) + 
      geom_bar(stat="identity", width=0.8) + geom_abline(slope=0) + 
      scale_color_manual(values=col.df$line_colors) +
      scale_fill_manual(values=alternate.white.and.fill(col.df$fill_colors)) +
      scale_x_discrete(labels=expt_labels) + 
      theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
      xlab("") + ylab("Amino Acid Uptake\n(nmol / uL cell / hr)") +
      facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
    
    if (save == TRUE) {
      compFlux.filename <- paste(directory, filename, ".Comp.Timecourses.AbsoluteFluxes.pdf", sep="")
      ggsave(compareFlux.plot, file=compFlux.filename, width=plot_width, height=plot_height, unit="cm")
    }
  }
  compareFlux.plot
}

# EXAMPLES:
# Compare_Timecourses_and_Fluxes(D3_Akt_Ras_fluxes, AAs = c(AAtoplot), expt_labels=BMK_labels, normTo1 = FALSE,
#                                save=TRUE, directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")

Compare_Timecourses_Stagger <- function(multi_expt_fluxes, AAs, expt_labels, normTo1=TRUE, 
                                        save=FALSE, directory="ProtScavFluxCalculation/Figures/", filename="choose_filename",
                                        plot_width=7, plot_height=6) {
  col.df <- scav.plot.colors %>% filter(amino.acid %in% AAs)
  
  multi_expt_fluxes <- multi_expt_fluxes %>% filter(compound %in% AAs)
  
  multi.fluxes.toCompare <- multi_expt_fluxes %>% group_by(data) %>% 
    gather(entryRoute, aa.flux, resulting.aa.flux, uptake.rate) %>% select(compound, data, entryRoute, aa.flux) %>% 
    mutate(comp_entryRoute = paste(compound, entryRoute, sep="_")) %>%
    group_by(compound, data) %>% mutate(norm.flux = aa.flux/sum(aa.flux))
  
  multi.fluxes.toCompare$compound <- factor(multi.fluxes.toCompare$compound, 
                                            levels=scav.plot.colors$amino.acid[scav.plot.colors$amino.acid %in% AAs])
  multi.fluxes.toCompare$comp_entryRoute <- factor(multi.fluxes.toCompare$comp_entryRoute, 
                                                   levels=comp_entry_levels[which(comp_entry_levels %in% unique(multi.fluxes.toCompare$comp_entryRoute))])
  
  compareFlux.plot <- ggplot(multi.fluxes.toCompare, aes(x=data, y=aa.flux*1000, color=compound, fill=comp_entryRoute)) + 
    geom_bar(stat="identity", position="dodge", width=0.8) + geom_abline(slope=0) + 
    scale_color_manual(values=col.df$line_colors) +
    scale_fill_manual(values=alternate.white.and.fill(col.df$fill_colors)) +
    scale_x_discrete(labels=expt_labels) + 
    theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
    xlab("") + ylab("Amino Acid Uptake\n(nmol / uL cell / hr)") +
    facet_grid( ~ compound) + theme(strip.text.x = element_text(size=0))
  
  if (save == TRUE) {
    compFlux.filename <- paste(directory, filename, ".Comp.Timecourses.AbsoluteFluxes.pdf", sep="")
    ggsave(compareFlux.plot, file=compFlux.filename, width=plot_width, height=plot_height, unit="cm")
  }
  
  compareFlux.plot
}

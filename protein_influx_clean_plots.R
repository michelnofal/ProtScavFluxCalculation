# setwd("/Users/Michel/Desktop/Research/code/")

# inputs needed: aa.lab.list
# dat <- aa.lab.list$data # %>% filter(compound %in% aas.toPlot)
# preds <- aa.lab.list$preds # %>% filter(compound %in% aas.toPlot)

aa.lab.plotter <- function (aa, dat, preds, save=FALSE, directory="ProtScavFluxCalculation/Figures/", filename="choose_filename",
                            max.y=0.35, col = "orange3", shape.num = 25) {
  dat.toplot <- dat %>% filter(compound == aa)
  preds.toplot <- preds %>% filter(compound == aa)
  
  aa.lab.intra <- dat.toplot %>% filter(loc == "Intracellular")
  aa.preds.intra <- preds.toplot %>% filter(loc == "Intracellular") 
  
  aa.lab.media <- dat.toplot %>% filter(loc == "Extracellular")
  aa.preds.media <- preds.toplot %>% filter(loc == "Extracellular") 
  
  if (aa == "lysine") {
    col = green_outside; shape.num = 22;
  } else if (aa == "phenylalanine") {
    col = red_outside; shape.num = 21;
  } else if (aa == "threonine") {
    col = blue_outside; shape.num = 24;
  } else if (aa == "valine") {
    col = purple_outside; shape.num = 23;
  } 
  
  ret <- list()
  aa.lab.intra.plot <- ggplot(aa.lab.intra %>% filter(compound == aa), aes(x=time.pt, y=frac)) + 
    geom_point(size=1.75, color=col, shape=shape.num) + 
    geom_line(data=aa.preds.intra %>% filter(compound == aa), aes(x=time.pt, y=predicted), color=col, lty=2) + 
    scale_x_continuous(breaks=c(0,12,16,20,24), labels=c(0,12,"","",24)) + scale_y_continuous(limits=c(0,max.y)) +
    xlab("Time (hr)") + theme(axis.title.x = element_text(size=8), axis.title.y = element_blank(), axis.text = element_text(size=7))
  ret$intra.plot <- aa.lab.intra.plot
  
  aa.lab.media.plot <- ggplot(aa.lab.media %>% filter(compound == aa), aes(x=time.pt, y=frac)) + 
    geom_point(size=1.75, color=col, shape=shape.num) + 
    geom_line(data=aa.preds.media %>% filter(compound == aa), aes(x=time.pt, y=predicted), color=col, lty=2) + 
    scale_x_continuous(breaks=c(0,12,16,20,24), labels=c(0,12,"","",24)) + scale_y_continuous(limits=c(0,max.y)) +
    xlab("Time (hr)") + theme(axis.title.x = element_text(size=8), axis.title.y = element_blank(), axis.text = element_text(size=7))
  ret$media.plot <- aa.lab.media.plot
  
  if(save == TRUE) {
    intra.filename <- paste(directory, filename, ".", aa, ".intra.pdf", sep="")
    ggsave(ret$intra.plot, file=intra.filename, width=3.5, height=3.5, unit="cm")
    
    media.filename <- paste(directory, filename, ".", aa, ".media.pdf", sep="")
    ggsave(ret$media.plot, file=media.filename, width=3.5, height=3.5, unit="cm")
  }
  
  ret
}

# EXAMPLE: 
# aa.lab.plotter("threonine", dat, preds, save=TRUE,
#                 directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")


# TO SAVE LEGEND
aa.lab.sample.plot <- ggplot(aa.lab.intra, aes(x=time.pt, y=frac, color=compound, shape=compound)) + 
  geom_point(size=1.75) + scale_color_manual(values=c(green_outside, red_outside, blue_outside, purple_outside, "orange3")) +
  scale_shape_manual(values=c(22,21,24,23,25))

# pdf(file="../figures/Ras_timecourse_cleanFigs/geom_point_legend.pdf")
# legend_plotter(aa.lab.sample.plot)
# dev.off()

# inputs needed: fluxes.list
# fluxes.toPlot <- fluxes.list$fluxes

ProtScavFlux.plotter <- function(fluxes.toPlot, 
                                 save=FALSE, directory="ProtScavFluxCalculation/Figures/", filename="choose_filename") {
  umol.uptakePlot <- ggplot(fluxes.toPlot, aes(x=compound, y=umol.protein.flux*10^6, fill=compound)) + 
    geom_bar(stat="identity", position="dodge", color="black", width=0.8) + 
    scale_fill_manual(values=c(green_outside, red_outside, blue_outside, "orange3", purple_outside)) +
    scale_x_discrete(labels=c("Lysine","Phenylalanine","Threonine","Tyrosine","Valine")) + 
    geom_abline(height=0, slope=0) + 
    theme(axis.text.x=element_text(size=8, angle=45,hjust=1),
          axis.title.y=element_text(size=8)) + 
    xlab("") + ylab("Scavenging Flux\n(pmol / uL cell / hr)")
  
  if (save == TRUE) {
    ProtScavFlux.filename <- paste(directory, filename, ".ProteinScavFlux.pdf", sep="")
    ggsave(umol.uptakePlot, file=ProtScavFlux.filename, width=4.5, height=6, unit="cm")
  }
  umol.uptakePlot
}

# EXAMPLE: 
# ProtScavFlux.plotter(fluxes.toPlot, save=TRUE,
#                directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")

# inputs needed: fluxes.list
# fluxes.toPlot <- fluxes.list$fluxes

FluxComparison.Plotter <- function(fluxes.toPlot, normTo1=TRUE,
                                   save=FALSE, directory="ProtScavFluxCalculation/Figures/", filename="choose_filename") {
  fluxes.toCompare <- fluxes.toPlot %>% gather(entryRoute, aa.flux, resulting.aa.flux, uptake.rate) %>% 
    select(compound, entryRoute, aa.flux) %>% 
    mutate(comp_entryRoute = paste(compound, entryRoute, sep="_")) %>%
    group_by(compound) %>%
    mutate(norm.flux = aa.flux/sum(aa.flux))  
  if (normTo1) {
    compareFlux.plot <- ggplot(fluxes.toCompare, aes(x=compound, y=norm.flux, color=compound, fill=comp_entryRoute)) + 
      geom_bar(stat="identity", width=0.8) + geom_abline(height=0, slope=0) + 
      scale_color_manual(values=c(green_outside, red_outside, blue_outside, "orange3", purple_outside)) +
      scale_fill_manual(values=c("white",green_inside,"white",red_inside,"white",blue_inside,"white",
                                 "orange1","white",purple_inside)) +
      scale_x_discrete(labels=c("Lysine","Phenylalanine","Threonine","Tyrosine","Valine")) + 
      theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
      xlab("") + ylab("Fractional\nAmino Acid Uptake")
    
    if (save == TRUE) {
      compFlux.filename <- paste(directory, filename, ".CompFractionalFluxes.pdf", sep="")
      ggsave(compareFlux.plot, file=compFlux.filename, width=4.8, height=6, unit="cm")
    }
  } else {
    compareFlux.plot <- ggplot(fluxes.toCompare, aes(x=compound, y=aa.flux*1000, color=compound, fill=comp_entryRoute)) + 
      geom_bar(stat="identity", width=0.8) + geom_abline(height=0, slope=0) + 
      scale_color_manual(values=c(green_outside, red_outside, blue_outside, "orange3", purple_outside)) +
      scale_fill_manual(values=c("white",green_inside,"white",red_inside,"white",blue_inside,"white",
                                 "orange1","white",purple_inside)) +
      scale_x_discrete(labels=c("Lysine","Phenylalanine","Threonine","Tyrosine","Valine")) + 
      theme(axis.text.x=element_text(size=8, angle=45,hjust=1), axis.title.y=element_text(size=8)) + 
      xlab("") + ylab("Amino Acid Uptake\n(nmol / uL cell / hr)")
    
    if (save == TRUE) {
      compFlux.filename <- paste(directory, filename, ".CompAbsoluteFluxes.pdf", sep="")
      ggsave(compareFlux.plot, file=compFlux.filename, width=4.5, height=6, unit="cm")
    }
  }
  compareFlux.plot
}

# EXAMPLE: 
# FluxComparison.Plotter(fluxes.toPlot, normTo1=FALSE, save=TRUE,
#                        directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")
# FluxComparison.Plotter(fluxes.toPlot, normTo1=TRUE, save=TRUE,
#                        directory="ProtScavFluxCalculation/Figures/BMK_White_Lab/", filename="BMK_Ras_LabPlot")


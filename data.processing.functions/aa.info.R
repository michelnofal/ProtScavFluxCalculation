# DMEM concentration of each AA
generate.aa.table <- function(essAA.perc = 40, gln.perc = 150, ser.perc = 100, perc.fbs.added = 5) {
  fbs.ps.dil = 100/(100+perc.fbs.added+1)
  aa_table <- data.table(compound = c("lysine","phenylalanine","threonine","tyrosine","valine","histidine"), 
                         DMEM.conc = c(800,400,800,500,800,200)*(essAA.perc/100)*fbs.ps.dil)
  unlab_aa_table <- data.table(compound = c("glutamine","arginine","methionine","tryptophan"),
                               DMEM.conc = c(4000*(gln.perc/100),400,200,78)*fbs.ps.dil)
  more_unlab_aa_table <- data.table(compound = c("serine","glycine"),
                                    DMEM.conc = c(400*(ser.perc/100), 400)*fbs.ps.dil)
  rbind(aa_table, unlab_aa_table, more_unlab_aa_table)
}
default.aa.table <- generate.aa.table()

# number of each AA in bovine serum albumin
alpha.df <- data.frame(compound = c("lysine","phenylalanine","threonine","tyrosine","valine","histidine",
                                    "arginine","glutamine","methionine","tryptophan","glycine","serine","alanine"), 
                       alpha = c(59, 27, 33, 20, 36, 17, 23, 20, 4, 2, 16, 28, 47))

# ... for labeling of facet titles in ggplot2
aa.labeller.complete <- function(variable, value) { 
  aa.names <- list('lysine'="Lysine",
                   'phenylalanine'="Phenylalanine",
                   'threonine'="Threonine",
                   'tyrosine'="Tyrosine",
                   'valine'="Valine",
                   'histidine'="Histidine",
                   'glutamine'="Glutamine",
                   'arginine'="Arginine",
                   'tryptophan'="Tryptophan",
                   'methionine'="Methionine",
                   'serine'="Serine",
                   'glycine'="Glycine")
  return(aa.names[value]) 
}
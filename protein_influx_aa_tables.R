# data on DMEM concentration

aa_table <- data.table(compound = c("lysine","phenylalanine","threonine","tyrosine","valine","histidine"), 
                       DMEM.conc = c(800,400,800,500,800,200)/1.06*(media_perc/100))
unlab_aa_table <- data.table(compound = c("glutamine","arginine","methionine","tryptophan"),
                             DMEM.conc = c(4000*gln_conc_factor,400,200,78)/1.06)
more_unlab_aa_table <- data.table(compound = c("serine","glycine"),
                                  DMEM.conc = c(400, 400)/1.06)

# alpha.df <- data.frame(compound = c("lysine","phenylalanine","threonine","tyrosine","valine","histidine"), 
#                        alpha = c(59, 27, 33, 20, 36, 17))

alpha.df <- data.frame(compound = c("lysine","phenylalanine","threonine","tyrosine","valine","histidine",
                                    "arginine","glutamine","methionine","tryptophan"), 
                       alpha = c(59, 27, 33, 20, 36, 17, 23, 20, 4, 2))

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
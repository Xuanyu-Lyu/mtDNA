# Create a function as a wrapper of the plot.pedigree function in kinship2 package
# Reference: kinship2 package by XX
# Parameters:
# ped: a pedigree simulated from SimPed function
# cex: the size of the text in the graph




PlotPedigree <- function(ped, cex = .5){
      library(kinship2)
      p <- ped[,-c(3,6)]
      colnames(p) <- c("ped","id","father","mother","sex")
      p[is.na(p)] <- 0 
      p$ped <- 1
      p$affected <- 0
      p$avail <- 0
      p2 <- pedigree(id=p$id, 
                     dadid=p$father, 
                     momid=p$mother, 
                     sex=p$sex, 
                     famid=p$ped)
      p3 <- p2['1']
      print(p3)
      return(plot.pedigree(p3,
                           cex = cex))
}

#PlotPedigree(SimPed(kpc = 2, Ngen = 6, marR = .8))

### prepare the fixed pedigrees and the dataframe for variance components

source("~/R-Project/mtDNA_mt2/Functions/evenInsert.R")
source("~/R-Project/mtDNA_mt2/Functions/famSizeCal.R")
source("~/R-Project/mtDNA_mt2/Functions/SimPed.R")

df_ped <- read.csv("~/R-Project/mtDNA_mt2/InitialData/Pedigrees_a2mt2.csv")
df_var <- read.csv("~/R-Project/mtDNA_mt2/InitialData/VarianceComb_a2mt2.csv")

set.seed(62)

l_ped <- list()
for(i in 1: nrow(df_ped)){
    ped_temp <- SimPed(kpc = df_ped$k[i],
                       Ngen = df_ped$G[i],
                       sexR = df_ped$p[i],
                       marR = df_ped$r[i])
    assign(paste0("ped",i),
           ped_temp)
    l_ped[[i]] <- ped_temp
    names(l_ped)[i] <- paste0("ped",i)
}


save.image("~/R-Project/mtDNA_mt2/InitialData/FixedPedVar_a2mt2.RData")

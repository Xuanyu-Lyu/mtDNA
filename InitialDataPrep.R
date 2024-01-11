### prepare the fixed pedigrees and the dataframe for variance components
#path_to_scripts <- "C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Functions/"
path_base = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/"
path_to_scripts <- paste0(path_base,"Functions/")
# Get a list of all .R files in the directory
r_scripts <- list.files(path = path_to_scripts, pattern = "\\.R$")
# Loop through the list of scripts and source each one
for(script in r_scripts) {
  source(paste0(path_to_scripts, script))
}

df_ped <- read.csv(paste0(path_base, "Pedigrees_a2mt2.csv"))
df_var <- read.csv(paste0(path_base, "VarianceComb_a2mt2.csv"))

# set.seed(62)

# l_ped <- list()
# for(i in 1: nrow(df_ped)){
# ped_temp <- SimPed(kpc = df_ped$k[i],
#                     Ngen = df_ped$G[i],
#                     sexR = df_ped$p[i],
#                     marR = df_ped$r[i])
# assign(paste0("ped",i),
#         ped_temp)
# l_ped[[i]] <- ped_temp
# names(l_ped)[i] <- paste0("ped",i)
# }


#save.image("~/R-Project/mtDNA_mt2/InitialData/FixedPedVar.RData")


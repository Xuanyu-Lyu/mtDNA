## Run the simulations
#source("C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/InitialDataPrep.R")
source("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/InitialDataPrep.R")

for (i in 1:1){
    for(j in 7:11){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/1214CheckSamePed","/","c",i,"p",j) 
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        }
        # do.call(RunSim, as.list(df_var[i,],
        #                         paste0("ped",j),
        #                         target_folder))
        RunSim_rd(Var = df_var[i,],
                  kpc = df_ped$k[j],
                  Ngen = df_ped$G[j],
                  sexR = df_ped$p[j],
                  marR = df_ped$r[j],
                  n = 5000,
                  path_results = target_folder)
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}

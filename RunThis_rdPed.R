## Run the simulations

source("C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/InitialDataPrep.R")

for(i in 1:5){
    
    for(j in 1:5){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/Data/1130Prelim","/","c",i,"p",j) 
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
                  n = 1000,
                  path_results = target_folder)
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}

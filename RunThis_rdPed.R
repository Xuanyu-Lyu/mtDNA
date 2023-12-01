## Run the simulations

source("~/R-Project/mtDNA_mt2/InitialData/InitialDataPrep.R")
source("~/R-Project/mtDNA_mt2/Functions/RunSim_rdped_full.R")

for(i in 2){
    
    for(j in 8){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("~/R-Project/mtDNA_mt2/Result_rdped_full","/","c",i,"p",j) 
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        }
        # do.call(RunSim, as.list(df_var[i,],
        #                         paste0("ped",j),
        #                         target_folder))
        RunSim_rd_full(Var = df_var[i,],
                  kpc = df_ped$k[j],
                  Ngen = df_ped$G[j],
                  sexR = df_ped$p[j],
                  marR = df_ped$r[j],
                  n = 10000,
                  path_results = target_folder)
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}

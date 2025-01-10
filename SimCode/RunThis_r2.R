## Run the simulations

source("~/R-Project/mtDNA_mt2/InitialData/InitialDataPrep.R")
source("~/R-Project/mtDNA_mt2/Functions/RunSim_20k.R")

for(i in 1: 5){
    
    for(j in 1: 10){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("~/R-Project/mtDNA_mt2/Result_20k","/","c",i,"p",j)
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        }
        # do.call(RunSim, as.list(df_var[i,],
        #                         paste0("ped",j),
        #                         target_folder))
        RunSim_20k(Var = df_var[i,],
               Ped = l_ped[[j]],
               n = 20000,
               path_results = target_folder )
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}




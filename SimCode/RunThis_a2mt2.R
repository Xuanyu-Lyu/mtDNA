## Run the simulations

load("~/R-Project/mtDNA_mt2/InitialData/FixedPedVar_a2mt2.RData")
source("~/R-Project/mtDNA_mt2/Functions/RunSim_a2mt2e2only.R")
options(expressions = 5e5)
# Run the code in terminal
# set R_CStackLimit=10000000

for(i in 1: 5){
    
    for(j in 1: 10){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("~/R-Project/mtDNA_mt2/ResultA2Mt2_20k_2","/","c",i,"p",j)
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        }
        # do.call(RunSim, as.list(df_var[i,],
        #                         paste0("ped",j),
        #                         target_folder))
        RunSim_a2mt2e2only(Var = df_var[i,],
                           Ped = l_ped[[j]],
                           n = 20000,
                           path_results = target_folder )
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}




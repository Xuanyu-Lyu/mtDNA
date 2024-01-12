## Run the simulations
#source("C:/Users/lxy75/OneDrive/Documents/R-Project/mtDNA/InitialDataPrep.R")
source("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/InitialDataPrep.R")


# # impact of mutation
# for (i in 2){
#     for(j in 12){
#         print(Sys.time())
#         cat(paste("start c",i,"p",j, "\n"))
#         target_folder <- paste0("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/Mutation","/","c",i,"p",j) 
#         if (!dir.exists(target_folder)){
#             dir.create(target_folder)}
#         for(loop in 1:100){
#             target_folder_loop <- paste0(target_folder,"/","loop",loop)
#             if (!dir.exists(target_folder_loop)){
#             dir.create(target_folder_loop)}
#             cat(paste("running loop",loop, "\n"))
#             RunSim_Mutation(
#                   Var = df_var_full[i,],
#                   kpc = df_ped$k[j],
#                   Ngen = df_ped$G[j],
#                   sexR = df_ped$p[j],
#                   marR = df_ped$r[j],
#                   n = 20000,
#                   path_results = target_folder_loop,
#                   gen_drop = 2,
#                   sex_drop = "F",
#                   n_drop = 1)
#         }      
#         print(Sys.time())
#         cat(paste("end c",i,"p",j, "\n"))
#     }
# }

# impact of missing links
for (i in 2){
    for(j in 12){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/MissingLink","/","c",i,"p",j) 
        if (!dir.exists(target_folder)){
            dir.create(target_folder)}
        for(loop in 1:50){
            target_folder_loop <- paste0(target_folder,"/","loop",loop)
            if (!dir.exists(target_folder_loop)){
            dir.create(target_folder_loop)}
            cat(paste("running loop",loop, "\n"))
            RunSim_MissingLink(
                  Var = df_var_full[i,],
                  kpc = df_ped$k[j],
                  Ngen = df_ped$G[j],
                  sexR = df_ped$p[j],
                  marR = df_ped$r[j],
                  n = 20000,
                  path_results = target_folder_loop,
                  gen_drop = 2,
                  sex_drop = "F",
                  n_drop = 1)
        }      
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}

# Power Run
for (i in 1:3){
    for(j in 1:5){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/20240112PowerRun","/","c",i,"p",j) 
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
                  n = 20000,
                  path_results = target_folder)
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}
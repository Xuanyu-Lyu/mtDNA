## test 

files <- list.files(path = "~/R-Project/mtDNA_mt2/Functions", pattern = "*.R$", full.names = TRUE)
for (file in files){
    source(file)
}

set.seed(1)
ped1 <- SimPed(3,5)
set.seed(2)
ped2 <- SimPed(3,5)
set.seed(3)
ped3 <- SimPed(3,5)
set.seed(4)
ped4 <- SimPed(3,5)
set.seed(5)
ped5 <- SimPed(3,5)
l_ped <- list(ped1, ped2, ped3, ped4, ped5)
df_var <- read.csv("~/R-Project/mtDNA_mt2/InitialData/VarianceComb_a2mt2.csv")

PlotPedigree(ped1)
Wide2Long(ped1)
write.csv(Wide2Long(ped1),"~/R-Project/mtDNA_mt2/Result_Check/Analysis_Check/ped1_long.csv" )

PlotPedigree(ped2)
Wide2Long(ped2)
write.csv(Wide2Long(ped2),"~/R-Project/mtDNA_mt2/Result_Check/Analysis_Check/ped2_long.csv" )

PlotPedigree(ped3)
Wide2Long(ped3)
write.csv(Wide2Long(ped3),"~/R-Project/mtDNA_mt2/Result_Check/Analysis_Check/ped3_long.csv" )

PlotPedigree(ped4)
Wide2Long(ped4)
write.csv(Wide2Long(ped4),"~/R-Project/mtDNA_mt2/Result_Check/Analysis_Check/ped4_long.csv" )

PlotPedigree(ped5)
Wide2Long(ped5)
write.csv(Wide2Long(ped5),"~/R-Project/mtDNA_mt2/Result_Check/Analysis_Check/ped5_long.csv" )

for(i in 2){
    
    for(j in 1: 5){
        print(Sys.time())
        cat(paste("start c",i,"p",j, "\n"))
        target_folder <- paste0("~/R-Project/mtDNA_mt2/Result_Check","/","c",i,"p",j)
        if (!dir.exists(target_folder)){
            dir.create(target_folder)
        }
        # do.call(RunSim, as.list(df_var[i,],
        #                         paste0("ped",j),
        #                         target_folder))
        RunSim_a2mt2e2only(Var = df_var[i,],
               Ped = l_ped[[j]],
               path_results = target_folder )
        print(Sys.time())
        cat(paste("end c",i,"p",j, "\n"))
    }
}


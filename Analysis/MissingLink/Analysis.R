setwd("~/R-Project/mtDNA_mt2/MissLink")
path_results <- "~/R-Project/mtDNA_mt2/MissLink"
#load(paste0(path_results,paste0("/modelSmr_",1,".RData")))
df_result <- data.frame("a2"=as.numeric(rep(NA,10)), "mt2"=as.numeric(rep(NA,10)),
                        "a2_miss"=as.numeric(rep(NA,10)),
                        "mt2_miss"=as.numeric(rep(NA,10))) 

for (count in 1: 100){
    load(paste0(path_results,paste0("/modelSmr_",count,".RData")))
    df_result[count,1] <- smr1$parameters$Estimate[1]
    df_result[count,2] <- smr1$parameters$Estimate[2]    
    df_result[count,3] <- smr2$parameters$Estimate[1]
    df_result[count,4] <- smr2$parameters$Estimate[2] 
}
apply(df_result,2,mean)
apply(df_result,2,sd)

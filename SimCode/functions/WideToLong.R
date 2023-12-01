### a function to convert relatedness matrix into long-form tables


Wide2Long <- function(Ped){
    library(reshape2)
    library(dplyr)
    df_melted_a2 <- melt(as.matrix(ped2add(Ped)), id.vars = NULL)
    colnames(df_melted_a2) <- c("ID1","ID2","a2")
    df_melted_mt2 <- melt(ped2mt_v3(Ped), id.vars = NULL)
    colnames(df_melted_mt2) <- c("ID1","ID2","mt2")
    df_long <- cbind(df_melted_a2, "mt2" = df_melted_mt2$mt2)
    df_count <- df_long %>% group_by(a2, mt2) %>% summarise(n = n())
    return(df_count)
    
}

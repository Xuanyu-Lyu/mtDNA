setwd("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/MissingLink/c2p12")
path_results2 <- "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/MissingLink/c2p12"
#load(paste0(path_results,paste0("/modelSmr_",1,".RData")))
df_result <- data.frame("a2_miss"=as.numeric(rep(NA,50)), 
                        "mt2_miss"=as.numeric(rep(NA,50)),
                        "a2"=as.numeric(rep(NA,50)),
                        "mt2"=as.numeric(rep(NA,50))) 



# collect the results of the 50 loops in one data frame
for (loop in 1:50){
    print(loop)
    load(paste0(path_results2,"/loop",loop,"/modelSmr.RData"))
    df_result[loop,1] <- smr1$parameters[1,5]
    df_result[loop,2] <- smr1$parameters[4,5]
    df_result[loop,3] <- smr2$parameters[1,5]
    df_result[loop,4] <- smr2$parameters[4,5]
}

# plot the distribution of the mt2_miss and mt2
library(ggplot2)
#library(Cairo)
# Reshaping the dataframe for plotting
df_long <- reshape2::melt(df_result, measure.vars = c("mt2_miss", "mt2"))
# Combined histogram
g1 =ggplot(df_long, aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5) +
    labs(x = "Value", y = "Density") +
    ggtitle("Density Plots of mt2_miss and mt2") +
    scale_fill_discrete(name = "Variable")

save_path = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/MissingLink"

ggsave( paste0(save_path,"/graph_distribution.png"), g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)



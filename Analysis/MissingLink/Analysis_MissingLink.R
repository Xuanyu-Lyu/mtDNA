setwd("/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/MissingLink/c2p12")
path_results2 <- "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/MissingLink/c2p12"
#load(paste0(path_results,paste0("/modelSmr_",1,".RData")))
df_result <- data.frame("j2_miss"=as.numeric(rep(NA,100)), 
                        "mt2_miss"=as.numeric(rep(NA,100)),
                        "j2"=as.numeric(rep(NA,100)),
                        "mt2"=as.numeric(rep(NA,100))) 



for (loop in 1:500){
    print(loop)
    load(paste0(path_results2,"/loop",loop,"/modelSmr.RData"))
    df_result[loop,1] <- smr1$parameters[5,5]
    df_result[loop,2] <- smr1$parameters[4,5]
    df_result[loop,3] <- smr2$parameters[5,5]
    df_result[loop,4] <- smr2$parameters[4,5]
}

# check the pedigree size and number of pedigrees
ncol(smr1[["dataSummary"]][["fam1"]])*length(smr1[["dataSummary"]])
ncol(smr1[["dataSummary"]][["fam1"]])
length(smr1[["dataSummary"]])
psych::describe(df_result, trim = 0) |> print(digits = 4)
# percentage of mt2 overestimation
(psych::describe(df_result, trim = 0)[4,3] - psych::describe(df_result, trim = 0)[2,3])/psych::describe(df_result, trim = 0)[4,3]
# percentage of am2 underestimation
(psych::describe(df_result, trim = 0)[3,3] - psych::describe(df_result, trim = 0)[1,3])/psych::describe(df_result, trim = 0)[3,3]

# plot the distribution of the mt2_miss and mt2
library(ggplot2)
#library(Cairo)
# Reshaping the dataframe for plotting
df_long <- reshape2::melt(df_result, measure.vars = c("mt2_miss", "mt2"))
# Combined histogram
g1 =ggplot(df_long, aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5) +
    labs(x = "Value", y = "Density") +
    ggtitle("Density Plots of correcly specified M and misspecified M") +
    scale_fill_manual(values = c("mt2_miss" = "#332288", "mt2" = "#E69F00"), name = c("MtDNA Matrix"), 
                      labels = c("mt2_miss" = "Misspecified M", "mt2" = "Correct M")) +
    theme_bw() +  # Start with a minimal theme
     theme(panel.background = element_rect(fill = "white", colour = "white"),  # Set panel background to white
         panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank(),  # Remove minor grid lines
         axis.line = element_line(colour = "black"))  # Add axis lines for clarity


save_path = "/Users/xuly4739/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/MissingLink"

ggsave( paste0(save_path,"/MissingLinkPlot.png"), g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)



setwd("/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/20240219TypeIErr/c2p12")
path_results2 <- "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Data/20240219TypeIErr/c2p12"
#load(paste0(path_results,paste0("/modelSmr_",1,".RData")))
df_result <- data.frame("mt2_typeIerror"=as.numeric(rep(NA,50)), 
                        "mt2_typeIerrorJ"=as.numeric(rep(NA,50))) 



# collect the results of the 50 loops in one data frame
for (loop in 1:50){
    print(loop)
    load(paste0(path_results2,"/loop",loop,"/modelSmr.RData"))
    df_result[loop,1] <- smr1$parameters[4,5]
    df_result[loop,2] <- smr1$parameters[5,5]

}
for (loop in 101:227){
    print(loop)
    load(paste0(path_results2,"/loop",loop,"/modelSmr.RData"))
    df_result[loop,1] <- smr1$parameters[4,5]
    df_result[loop,2] <- smr1$parameters[5,5]

}
df_result <- df_result[-(51:100),]
# calculate the proportion of values that within the range of -.01 and .01 in a seperate vector
sum(df_result$mt2_typeIerror == 0)
sum(df_result$mt2_typeIerror > -.03 & df_result$mt2_typeIerror < .03)
freq01 <- sum(df_result$mt2_typeIerror > -.03 & df_result$mt2_typeIerror < .03)/length(df_result$mt2_typeIerror)
freq01
sum(df_result$mt2_typeIerror > -.005 & df_result$mt2_typeIerror < .005)
freq005 <- sum(df_result$mt2_typeIerror > -.005 & df_result$mt2_typeIerror < .005)/length(df_result$mt2_typeIerror)
freq005

freq005 <- sum(df_result$mt2_typeIerror > -.001 & df_result$mt2_typeIerror < .001)/length(df_result$mt2_typeIerror)
freq005

sum(df_result$mt2_typeIerrorJ > -.03 & df_result$mt2_typeIerrorJ < .03)/length(df_result$mt2_typeIerrorJ)

print(psych::describe(df_result, trim = 0),digits = 4)
# plot the distribution of the mt2_miss and mt2
library(ggplot2)

# Combined histogram
g1 =ggplot(df_result, aes(x = mt2_typeIerror)) +
    geom_histogram(fill = "#E41A1C") +
    geom_density(alpha = 0.4, fill = "#E69F00", color = "white") +
    labs(x = "Value", y = "Count") +
    theme_bw() +    # Start with a minimal theme
    theme(panel.background = element_rect(fill = "white", colour = "white"),  # Set panel background to white
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          axis.line.x = element_line(colour = "black"),  # Add x axis line
          axis.line.y = element_line(colour = "black"),  # Add y axis line
          panel.border = element_blank())  # Remove panel border

save_path = "/Users/lyux20/Library/CloudStorage/OneDrive-UCB-O365/Documents/mtDNA/mtDNA/Analysis/20240219TypeIErr"

ggsave( paste0(save_path,"/TypeIerror.png"), g1, width = 6, height = 4.5,  type = "cairo-png", dpi = 900)




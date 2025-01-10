library(ggplot2)

# Number of simulations
n_sims <- 100000

# Function to generate samples from the mixture distribution
generate_samples <- function(ncp=0, n=n_sims) {
    chi_samples <- rchisq(n, df=1, ncp=ncp)
    zero_one <- sample(c(0, 1), size=n, replace=TRUE)
    return(chi_samples * zero_one)
}

# Generate samples for null and alternative hypothesis
null_samples <- generate_samples()
alt_samples <- generate_samples(ncp=5)  # You can change this to a different lambda value

# Estimate critical value
critical_value <- quantile(null_samples, 0.95)

# Create data frame for ggplot
data <- data.frame(
    Value = c(null_samples, alt_samples),
    Distribution = factor(rep(c("Null", "Alternative"), each=n_sims))
)

# Create the ggplot
p <- ggplot(data, aes(x=Value, fill=Distribution)) +
    geom_density(alpha=0.4, aes(y=..scaled..), size =.8) +
    geom_vline(aes(xintercept=critical_value), linetype="dashed", color="darkred", size = 1) +
    geom_vline(aes(xintercept=5), linetype="dashed", color="black", size = 1) +
    labs(x="Value", y="Scaled Density") +
    scale_fill_manual(values=c("#E41A1C", "#332288")) +
    xlim(0,20)+
    annotate("text", x = critical_value - 1, y = -.05, label = "Critical Value", angle = 0, size = 4) +
    annotate("text", x = 5 + 0.4, y = -.05, label = expression(lambda), angle = 0, size = 4) +
    theme(panel.background = element_rect(fill = "transparent"),
                         panel.grid = element_line(color = "transparent"),
                         axis.line = element_line(size = 1, colour = "black"),
                         #axis.line.y = element_blank(),
                         axis.text = element_text( color = "black"),
                         #axis.text.y = element_blank(),
                         #axis.ticks.y = element_blank(),
                         text=element_text( family="Calibri",  size = 12),
                         legend.spacing = unit(-17,'pt'),
                         legend.margin = margin(t=0,b=0,unit='pt'),
                         legend.background = element_rect(),
                         legend.position=c(.8,.2))

print(p)
ggsave( "~/R-Project/mtDNA_mt2/Analysis/power_distribution.png",p, width = 4.5, height = 4.5,  type = "cairo-png", dpi = 900)
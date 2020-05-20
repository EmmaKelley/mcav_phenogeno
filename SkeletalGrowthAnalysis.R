# set working directory (where you want files to be saved)
setwd("/Volumes/wrightlab/MCav_VibrioChallenge/skeletalGrowth/")

# load any packages you want to use
# install.packages("tidyverse") # you only have to do this once!!
library(tidyverse)

# load files (this can be from your working directory or any other directory, if you know the path)
weights <- read_csv("buoyantWeight.csv")
summary(weights)

# recast "genet" as a factor, not an integer
weights <- weights %>% mutate(genet = as.factor(genet), sam = as.factor(sam), 
                              treat = as.factor(treat), bank = as.factor(bank), pheno = as.factor(pheno))

# look at the file
head(weights)
tail(weights)

# summarize your data
# look for missing data or strange outliers
summary(weights)
hist(weights$weight_4_15)
hist(weights$weight_2_16)
hist(weights$Area)

# start analyzing! ------
# notice that I just made my first "chunk" by adding >=4 symbols after the title

# calculate the difference between the first and second time points
# calculate the difference as a percentage of starting weight
weights <- weights %>% 
        mutate(weightdiff = weight_2_16-weight_4_15) %>% 
        mutate(percDiff = (weightdiff/weight_4_15)*100)

hist(weights$weightdiff)
hist(weights$percDiff)


# Normalize weightdiff w/ regard to area- 

weights <- weights %>% 
        mutate(normweight = weightdiff/Area)

summary(weights$normweight)

hist(weights$normweight)

weights <- weights %>% 
        mutate(calcification = normweight/392)
hist(weights$calcification)

# stats using MCMC generalized linear mixed model -----
library(MCMCglmm)
set.seed(1)
mcmc.weight <- MCMCglmm(calcification~bank+pheno+treat, 
                        random=~genet, data=weights,
                        family="gaussian")
summary(mcmc.weight)
pMCMC.bank <- round(summary(mcmc.weight)$solutions[2,5],3)
pMCMC.pheno <- round(summary(mcmc.weight)$solutions[3,5],3)
pMCMC.treat <- round(summary(mcmc.weight)$solutions[4,5],3)


# plot with ggplot -----

# by bank -----
plot.bank <- weights %>% ggplot(aes(x = bank, y = calcification)) +
        # Make the box plot
        geom_boxplot(aes(fill = bank)) +
        # Change the colors of the boxplot and the name in the legend
        scale_fill_manual(values=c("navy", "darkgreen"),
                          name = "Bank",
                          labels = c( "East", "West")) +
        # Plot individual data points
        geom_point() +
        # Add horizontal jitter
        geom_jitter(width = 0.05) +
        # Change the colors of the points (if you want)
        scale_color_manual(values=c("black", "black")) +
        # Change the title and axis labels
        labs(title="Calcification Rate by Bank", 
             x="Bank",
             y=expression(Calcification~Rate~(g/cm^2/d))) +
        scale_x_discrete(labels = c("East", "West")) + 
        # Add annotation with p-value
        annotate("text", x = 2, y = 0, size = 5,
                 label = paste("pMCMC = ", pMCMC.bank)) +
        # Remove the grey background
        theme_bw() +
        # Change the sizes and placements of text and symbols
        theme(plot.title = element_text(size=35, hjust=0.5),
              axis.title.x = element_text(size=25),
              axis.title.y = element_text(size=25),
              axis.text.x = element_text(size=25),
              axis.text.y = element_text(size=25),
              legend.title = element_text(size=25),
              legend.text = element_text(size=25),
              legend.key.size = unit(3,"line"))
plot.bank

# by treatment -----
plot.treat <- weights %>% ggplot(aes(x = treat, y = calcification)) +
        # Make the box plot
        geom_boxplot(aes(fill = treat)) +
        # Change the colors of the boxplot and the name in the legend
        scale_fill_manual(values=c("grey60", "red4"),
                          name = "Treatment",
                          labels = c( "Control", "Vibrio")) +
        # Plot individual data points
        geom_point() +
        # Add horizontal jitter
        geom_jitter(width = 0.05) +
        # Change the colors of the points (if you want)
        scale_color_manual(values=c("black", "black")) +
        # Change the title and axis labels
        labs(title="Calcification Rate by Treatment", 
             x="Treatment",
             y=expression(Calcification~Rate~(g/cm^2/d))) +
        scale_x_discrete(labels = c("Control", expression(italic("Vibrio")))) + 
        # Add annotation with p-value
        annotate("text", x = 2, y = 0, size = 5,
                 label = paste("pMCMC = ", pMCMC.treat)) +
        # Remove the grey background
        theme_bw() +
        # Change the sizes and placements of text and symbols
        theme(plot.title = element_text(size=35, hjust=0.5),
              axis.title.x = element_text(size=25),
              axis.title.y = element_text(size=25),
              axis.text.x = element_text(size=25),
              axis.text.y = element_text(size=25),
              legend.title = element_text(size=25),
              legend.text = element_text(size=25),
              legend.key.size = unit(3,"line"))
plot.treat

# by pheno -----
plot.pheno <- weights %>% ggplot(aes(x = pheno, y = calcification)) +
        # Make the box plot
        geom_boxplot(aes(fill = pheno)) +
        # Change the colors of the boxplot and the name in the legend
        scale_fill_manual(values=c("coral3", "azure3"),
                          name = "Phenotype",
                          labels = c( "Resistant", "Susceptible")) +
        # Plot individual data points
        geom_point() +
        # Add horizontal jitter
        geom_jitter(width = 0.05) +
        # Change the colors of the points (if you want)
        scale_color_manual(values=c("black", "black")) +
        # Change the title and axis labels
        labs(title="Calcification Rate by Phenotype", 
             x="Phenotype",
             y=expression(Calcification~Rate~(g/cm^2/d))) +
        scale_x_discrete(labels = c("Resistant", "Susceptible")) +
        # Add annotation with p-value
        annotate("text", x = 2, y = 0, size = 5,
                 label = paste("pMCMC = ", pMCMC.pheno)) +
        # Remove the grey background
        theme_bw() +
        # Change the sizes and placements of text and symbols
        theme(plot.title = element_text(size=35, hjust=0.5),
              axis.title.x = element_text(size=25),
              axis.title.y = element_text(size=25),
              axis.text.x = element_text(size=25),
              axis.text.y = element_text(size=25),
              legend.title = element_text(size=25),
              legend.text = element_text(size=25),
              legend.key.size = unit(3,"line"))

plot.pheno

# plot all together
library(gridExtra)
combinedplot <- grid.arrange(plot.bank, plot.treat, plot.pheno, ncol=3)

# Save plot to working directory
ggsave(filename="CalcificationRateBigPlot-adj.jpg", plot=combinedplot,
       width = 67, height = 33, 
       units = "cm", # other options are "in", "cm", "mm" 
       dpi = 200)

# SAVE ------
save(weights, file="calcificationData.Rdata")
load("calcificationData.Rdata")

# TO DO --------
# can we "bin" growth rates for DESeq2?  -----

# step 1: make a per-genet average calcification rate
head(weights)
 means <- weights %>% group_by(genet) %>%
        summarize(meanCalc = mean(calcification),
                  meanNorm = mean(normweight))
means

# step 2: see if any natural breaks fall out (low, medium, high calcification?)
# maybe try rounding?
hist(round(means$meanCalc,4))
hist(means$meanCalc)
hist(means$meanNorm,10)
# play around with this
# is there any separation in bins that makes sense?


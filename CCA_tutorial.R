# Canonical Correlation Analysis (CCA) for brain-behaviour associations. 
# bethany.little@ncl.ac.uk 2023


##### 0. load packages #####

# load packages
library(tidyverse)
library(crayon)
library(Hmisc)
library(MVN)
library(CCA)
library(CCP)
library(ggseg)


##### 1. set up the environment and load data #####

# clear the environment
remove(list = ls())

# set the working directory to the source file location
CCAdir <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/")
setwd(CCAdir)

# load in datasets using read.csv (row.names=1 sets the first 
# column (subject IDs) in the csv file as the rownames of the dataframe)
brain <- read.csv("NIMH-IHV_corticalthickness.csv", row.names=1)
behav <- read.csv("NIMH-IHV_NIHtoolbox.csv", row.names=1)

# remove mean thickness variables from brain data
brain <- subset(brain, select=-c(lh_MeanThickness_thickness,rh_MeanThickness_thickness))

# save the number of variables in each dataset
nbrain <- ncol(brain)
nbehav <- ncol(behav)
# save sample size
sampleN <- nrow(brain)
# save variable names
vars.behav <- colnames(behav)
vars.brain <- colnames(brain)

# change the working directory to a folder to the Results folder
setwd(paste0(CCAdir,"/Results"))


##### 2. PCA #####

# If the number of (total) variables is larger than the sample size, then data 
# reduction techniques should be used to reduce the number of variables in one 
# or both datasets. Usually, there are many variables in the brain data, for 
# example, in this tutorial we have cortical thickness data for 68 regions 
# according to the Desikan-Killiany parcellation atlas. Therefore, we will run 
# Principal Components Analysis (PCA) on the brain data only here, however the 
# same steps can be applied to the behavioral data if required. This script 
# saves the number of principal components (PCs) that explin 90% of the data. 


###### 2. a) check data meets assumptions for PCA ######

# To use PCA, data should: 
# 1. be multiple continuous variables (or ordinal)
# 2. have linear relationship between all variables
# 3. have sampling adequacy (rule-of-thumb are a min of 150 cases)
# 4. have adequate correlations between the variables 
# 5. have no significant outliers (component scores >3SDs away from the mean) 
# 6. have multivariate normality
# The rest of the script assumes that the data meets these requirements. 
# Some suggestions of how to test some of these assumptions are below. 

# can test linear relationships with scatterplots of pairs of variables, e.g.: 
ggplot(brain, aes(lh_caudalmiddlefrontal_thickness,lh_posteriorcingulate_thickness)) + geom_point() + geom_smooth(method="lm", alpha=0.1, fill="Blue")
ggplot(brain, aes(lh_transversetemporal_thickness,lh_parsorbitalis_thickness)) + geom_point() + geom_smooth(method="lm", alpha=0.1, fill="Blue")

# can test if there are adequate correlations between the variables by creating a 
# correlation matrix - check all variables have at at least 1 correlation above 0.3
testCors <- data.frame(cor(brain, method = c("pearson")))
testCors

# check which variables have outliers (scores >3 SDs away from the mean)
outs <- c() 
for (i in 1:nbrain){
  varname <- colnames(brain)[i]
  var <- brain[[varname]]
  Nhigh <- sum(var > (mean(var, na.rm=T) + 3*sd(var, na.rm=T)), na.rm=T)
  Nlow <- sum(var < (mean(var, na.rm=T) - 3*sd(var, na.rm=T)), na.rm=T)
  if (Nhigh + Nlow > 0) {
    outs<- c(outs, varname)
  }
}
outs

# test for multivariate normality using Mardia's test of skewness and kurtosis
Mardia.brain <- MVN::mvn(data = brain, mvnTest ="mardia", univariateTest = "SW")
Mardia.brain$multivariateNormality


###### 2. b) run PCA ######

# run the PCA
pca.brain <- prcomp(brain)

# visually check how many components explain 90% of the data and plot PCs
summary(pca.brain)
plot(pca.brain)

# save the first X components that explain 90% of data
pca.brain.vars <- apply(pca.brain$x, 2, var)            # use var() to get the sample variance of each column (PC)
pca.brain.prop <- pca.brain.vars / sum(pca.brain.vars)  # calculate the proportion of sample variance 
pca.brain.cs <- cumsum(pca.brain.prop)                  # get the cumulative sum of the proportion of sample variance
pca.brain.N <- which(pca.brain.cs > 0.9)[1]             # find the first case that explains at least 90% of the cumulative variance explained
brain.PCs <- data.frame(pca.brain$x[,1:pca.brain.N])    # save only the PCs that together explain 90% of the variance in the data



##### 3. Canonical Correlation Analysis (CCA) #####

###### 3. a) check data meets assumptions for CCA ######

# CCA assumes that the data are suitable for correlation-based analysis. 
# Check both datasets have multivariate normality and no extreme outliers. 


### brain data ###

# list any brain variables (PCs) that are not normally distributed (shapiro-wilk 
# test p<.05); if there are none, the output will be NULL 
notNorm.brain.PC <- c()
for (i in 1:length(brain.PCs)){
  if (shapiro.test(brain.PCs[[i]])$p <.05){
    notNorm.brain.PC<- c(notNorm.brain.PC, colnames(brain.PCs[i]))
  }
}
notNorm.brain.PC 

# check multivariate normality of brain PCs dataset
Mardia.brain.PCs <- MVN::mvn(data = brain.PCs, mvnTest ="mardia", univariateTest = "SW")
Mardia.brain.PCs$multivariateNormality

# get list of variables with outliers (+-3 SD away from mean)
Outbrain.PCs <- c()
for (i in 1:length(brain.PCs)){
  var = brain.PCs[[i]]
  if (sum(var > mean(var) + 3*sd(var) | var < mean(var) - 3*sd(var))){
    Outbrain.PCs<- c(Outbrain.PCs, colnames(brain.PCs[i]))
  }
}
Outbrain.PCs

# if any of the PCs are not normal or have outliers, can visualise the data:
hist(brain.PCs$PC3)
boxplot(brain.PCs$PC3)


### behavioural data ###

# list any behav variables that are not normally distributed (shapiro-wilk test p<.05) 
notNorm.behav <- c()
for (i in 1:length(behav)){
  if (shapiro.test(behav[[i]])$p <.05){
    notNorm.behav<- c(notNorm.behav, colnames(behav[i]))
  }
}
notNorm.behav

# check multivariate normality of behav PCs dataset
Mardia.behav <- MVN::mvn(data = behav, mvnTest ="mardia", univariateTest = "SW")
Mardia.behav$multivariateNormality

# get list of variables with outliers (+-3 SD away from mean)
Outbehav <- c()
for (i in 1:length(behav)){
  var = behav[[i]]
  if (sum(var > mean(var) + 3*sd(var) | var < mean(var) - 3*sd(var))){
    Outbehav<- c(Outbehav, colnames(behav[i]))
  }
}
Outbehav

# if any of the PCs are not normal or have outliers, can visualise the data:
hist(behav[,1])
boxplot(behav[,1])


# *** PLEASE NOTE BEFORE PROCEEDING ***
# The following script assumes that the datasets have multiariate normality and 
# don't have any extreme outliers. If data does not meet these requirement, 
# transform variables accordingly. 


###### 3. b) run CCA ######

# before running CCA, check participants are in the same order in each dataset
if ((sum(rownames(brain)!=rownames(behav)))==0){
  cat(green("all subjects aligned"))} else {
    cat(red("ERROR - subjects don't match in brain and behaviour datasets"))}

# run CCA analysis using cc() function (behaviour data is dataset U, brain data is dataset V)
CCA <- CCA::cc(behav, brain.PCs)

# save canonical correlation coefficients
CCA.cor <- CCA$cor

# save the canonical variates for each subject for each dataset
CCA.subXscore <- data.frame(CCA$scores$xscores)
CCA.subYscore <- data.frame(CCA$scores$yscores)


### get the canonical loadings and cross-loadings ###

# save the canonical loadings
CCA.Uloadings <- data.frame(CCA$scores$corr.X.xscores)
CCA.Vloadings <- data.frame(CCA$scores$corr.Y.yscores)
# save the canonical cross-loadings 
CCA.Ucrossload <- data.frame(CCA$scores$corr.X.yscores)
CCA.Vcrossload <- data.frame(CCA$scores$corr.Y.xscores)

# To get the canonical loadings for the original brain variables (rather than 
# the PCs) with canonical variates from V, correlate the original brain matrix 
# (subjects*ROIs) with the canonical variates from V (subjects*VcanonicalVariates). 
CCA.brain.load <- data.frame(cor(brain, CCA.subYscore, method="pearson"))
CCA.brain.load.r <- data.frame(rcorr(as.matrix(brain), as.matrix(CCA.subYscore), type="pearson")$r)
CCA.brain.load.p <- data.frame(rcorr(as.matrix(brain), as.matrix(CCA.subYscore), type="pearson")$P)
CCA.brain.load.n <- data.frame(rcorr(as.matrix(brain), as.matrix(CCA.subYscore), type="pearson")$n)
# merge into one data frame and save correlation results
CCA.brain.load.all <- cbind(rbind(rep("rho", nbehav), data.frame(CCA.brain.load.r[1:nbrain,(nbrain+1):(nbrain+nbehav)])), 
                            rbind(rep("p", nbehav), CCA.brain.load.p[1:nbrain,(nbrain+1):(nbrain+nbehav)]))

# To get the canonical cross-loadings for the original brain variables (rather 
# than the PCs) with behavioural canonical variates (U), correlate the original 
# brain matrix (subjects*ROIs) with the canonical variates from U (subjects*UcanonicalVariates). 
CCA.braincrossLoad <- data.frame(cor(brain, CCA.subXscore, method="pearson"))
CCA.braincrossLoad.r <- data.frame(rcorr(as.matrix(brain), as.matrix(CCA.subXscore), type="pearson")$r)
CCA.braincrossLoad.p <- data.frame(rcorr(as.matrix(brain), as.matrix(CCA.subXscore), type="pearson")$P)
CCA.braincrossLoad.n <- data.frame(rcorr(as.matrix(brain), as.matrix(CCA.subXscore), type="pearson")$n)
# merge into one data frame and s Save correlation results (without group)
CCA.braincrossLoad.all <- cbind(rbind(paste("r_", c(1:nbehav), sep = ""), data.frame(CCA.braincrossLoad.r[1:nbrain,(nbrain+1):(nbrain+nbehav)])), 
                                rbind(paste("p_", c(1:nbehav), sep = ""), CCA.braincrossLoad.p[1:nbrain,(nbrain+1):(nbrain+nbehav)]))
colnames(CCA.braincrossLoad.all) <- CCA.braincrossLoad.all[1,] # set column names
CCA.braincrossLoad.all <- CCA.braincrossLoad.all[-c(1),] # delete first row

# save all outputs as csv files
write.csv(CCA.Uloadings,"CCA.Uloadings.csv")
write.csv(CCA.Vloadings,"CCA.Vloadings.csv")
write.csv(CCA.Ucrossload,"CCA.Ucrossload.csv")
write.csv(CCA.Vcrossload,"CCA.Vcrossload.csv")
write.csv(CCA.subXscore,"CCA.subXscore.csv")
write.csv(CCA.subYscore,"CCA.subYscore.csv")
write.csv(CCA.cor,"CCA.cor.csv")
write.csv(CCA.brain.load.all,"CCA.Vloadings.ROIs.csv")
write.csv(CCA.braincrossLoad.all,"CCA.Vcrossload.ROIs.csv")


##### 4. Test significance of canonical correlations #####

# run significance tests and save output in a txt file (here we use Wilks's test 
# as a common test of significance, but also Pillai's trace as a robust test)
cca.wilks <- p.asym(rho=CCA$cor, N=sampleN, p=nbehav, q=nbrain, tstat="Wilks")
cca.pillai <- p.asym(rho=CCA$cor, N=sampleN, p=nbehav, q=nbrain, tstat="Pillai")

# see p-values for each canonical correlation
cca.wilks$p.value
cca.pillai$p.value

# save output in txt file
sink("CCA_sigTests.txt")
cca.wilks
cca.pillai
sink() 


##### 5. Perform permutation test for CCA results (if first CCA looks ok) #####

# run the permutation tests using p.perm() function from the CCP package (run 
# both Wilk's test and Pillai's test of significance)
# here we are permuting the data 10,000 times
CCA.perm.W <- p.perm(behav, brain.PCs, nboot = 10000, rhostart = 1, type = "Wilks")
CCA.perm.P <- p.perm(behav, brain.PCs, nboot = 10000, rhostart = 1, type = "Pillai")

# save plots of the permutation distributions as pdf files
pdf(file="CCA.perm.Wilks.pdf")
plt.perm(CCA.perm.W)
dev.off()
pdf(file="CCA.perm.Pillai.pdf")
plt.perm(CCA.perm.P)
dev.off()

# save results of permutation tests in a txt file 
sink("CCA.perm.txt")
cat("Results of permutation test for CCA (untransformed behav; brain data not z-scored; nboot=10000)", "\n")
cat("\n", "Permutation test - Wilks", "\n")
cat(paste0("original value of the statistic (without resampling): ",round(CCA.perm.W$stat0,3)), "\n")
cat("summary of N=nboot permuted statistics: ")
print(summary(CCA.perm.W$stat))
cat(paste0("mean of N=nboot permuted statistics=",round(mean(CCA.perm.W$stat),3)), "\n")
cat(paste0("number of permutation resamplings that resulted in a more extreme value of the statistic than stat0: ", CCA.perm.W$nexcess), "\n")
cat(paste0("p-value derived from nexcess: ", round(CCA.perm.W$p.value,3)), "\n")
cat("\n", "Permutation test - Pillai", "\n")
cat(paste0("original value of the statistic (without resampling): ",round(CCA.perm.P$stat0,3)), "\n")
cat("summary of N=nboot permuted statistics: ")
print(summary(CCA.perm.P$stat))
cat(paste0("mean of N=nboot permuted statistics=",round(mean(CCA.perm.P$stat),3)), "\n")
cat(paste0("number of permutation resamplings that resulted in a more extreme value of the statistic than stat0: ", CCA.perm.P$nexcess), "\n")
cat(paste0("p-value derived from nexcess: ", round(CCA.perm.P$p.value,3)), "\n")
sink() 


##### 6. Visualize results #####

###### 6. a) scatterplot of first canonical correlation ######

# make scatterplot using ggplot2
CCA.1stScatter <- ggplot(brain, aes(x=as.matrix(CCA.subXscore[1]), 
                                    y=as.matrix(CCA.subYscore[1]))) + 
  geom_point() +
  labs(title="Scatterplot of the first canonical correlation",
       x="U1 (cogntitive data)",
       y="V1 (brain data)",
       subtitle=paste0("(n=", nrow(brain),")")) +
  theme_light() +
  theme(text = element_text(size = 16))
CCA.1stScatter

# save scatterplot as pdf file
pdf(file="CCA_1stcanobehavair_scatter.pdf")
CCA.1stScatter
dev.off()


###### 6. b) plot cross-loadings for each ROI ######


# Plots the results for the first canonical correlation only. 

# This part of the scripts assumes that the brain data was parcellated using 
# the desikan-killianny atlas. 
# Mowinckel & Vidal-Piñeiro (2020) provide a tutorial on how to use ggseg.

# get first canonical variate statistics and set up to use in ggseg
CCA.ggseg <- subset(CCA.braincrossLoad.all, select=c("r_1", "p_1"))
CCA.ggseg$r_1 <- as.numeric(CCA.ggseg$r_1)
CCA.ggseg$p_1 <- as.numeric(CCA.ggseg$p_1)
CCA.ggseg$label <- row.names(CCA.ggseg)
CCA.ggseg <- CCA.ggseg %>% mutate_at("label", str_replace, "_thickness", "") %>%  as_tibble()

# set limits for plotting CCA rho (sets it to .1 above/below from the max/min)
rholim <- max(c(abs(max(CCA.ggseg$r_1)), abs(min(CCA.ggseg$r_1))))
rholim <- round(rholim,1)+.1

# plot results using ggseg
CCA.ROIplot <- CCA.ggseg %>% 
  ggseg(mapping = aes(fill=r_1), position="stacked") + 
  labs(title = paste0("Cross-loadings of each brain region"), 
       fill = "r") + 
  scale_fill_gradient2(low = "red3", mid="white", high = "mediumblue", 
                       limits=c(-rholim,rholim)) +
  guides(fill = guide_colourbar(barheight = 6)) +
  theme_minimal() + theme(text = element_text(size = 12))
CCA.ROIplot

# save as pdf file
pdf("CCA.brain.pdf")
CCA.ROIplot
dev.off()


##### END #####
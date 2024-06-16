knitr::opts_chunk$set(echo = TRUE)

# Step 2: Predictor selection, GLM & GAM for richness
options(warn = -1)
setwd("/Users/malatji/Downloads/Masegos project")

#all data
env<-read.csv("new_env_var2.csv") #, sep=";")

average_richness <- mean(env$ Richness)
print(average_richness)

env$Site <- as.factor(env$Site)
env$Prec <-as.numeric(env$Prec)
env$Prec.of.driest.quarter<-as.numeric(env$Prec.of.driest.quarter)
env$Ele <-as.numeric(env$Ele)
env$Richness <-as.numeric(env$Richness)
env$Abundance <-as.numeric(env$Abundance)
env$Date <- as.factor(env$Date)
env$Collector <- as.factor(env$Collector)
env$Habitat <- as.factor(env$Habitat)

sbs<-read.csv("site by species.csv") #, sep=";") #columns seperated by ; before said only have 1 variable, introduced sep=

sbs=sbs[,4:233]
# sites aligned
identical(sbs[,1],env[,1]) 

#library(PerformanceAnalytics)
#library(car)
#library(mgcv)
library(ggplot2)
library(tidyr)

# pairwise correlations after removing categorical variables
cortable <- cor(env[,c("Temp", "Prec", "Ele", "Prec.of.driest.quarter", "Temp.annual.range","Richness")],method="spearman",use="complete.obs") 
heatmap(cortable)

write.csv(cortable,"cortable.csv",row.names = TRUE)
# 0 correlations > 0.75

cortable[6,] # correlations against Richness

en <- env[,c("Temp",
             "Prec",
             "Prec.of.driest.quarter",
             "Temp.annual.range",
             "Ele",
             "Richness",
             "Habitat")]
          
env$Temp <- (env$Temp)
env$Prec <- (env$Prec)
env$Prec.of.driest.quarter <- (env$Prec.of.driest.quarter)
env$Temp.annual.range <- (env$Temp.annual.range)
env$Ele <- (env$Ele)
env$Date <- (env$Date)
env$Collector <- (env$Collector)
env$Habitat <- (env$Habitat)

# check for Concurvity to ensure <0.8 (worst)
model.3 <-gam(env$Richness ~
              +s(Temp)
              +s(Prec)+
              Date+
              Habitat,
              ,data=env,method="REML")
cat("Deviance explained =",1-model.3$deviance/model.3$null.deviance)

#plot(model.3, all.terms=TRUE, pages = 3)
boxplot(env$Richness ~ env$Date,las = 2, names= c(2009,2011, 2013,2016,2017), col="red")
boxplot(env$Richness ~ env$Habitat,las=2, names= c("Cor","For","Gra","Pla","Thi"), col= "deepskyblue2")
#boxplot(env$Richness ~ env$Collector, las=2, names= c("Pryke", "Crous","Rene","Benie"),col ="darkgreen")
summary(model.3) #edf near 1 or 1 is linear
#hist(env$Temp)
concurvity(model.3)

#Make the plot with residuals
plot(model.3, residuals = TRUE)

# Change shape of residuals
plot(model.3, residuals = TRUE, pch = 1, cex = 1)

# Plot all effects
plot(model.3, residuals = TRUE, pch = 1, cex = 1, all.terms=TRUE, page=1)

# Make another plot adding the intercept value and uncertainty
plot(model.3, residuals = TRUE, pch = 1, cex = 1, all.terms=TRUE, page=1, shade = TRUE, shade.col = "purple", shift = coef(model.3)[1], seWithMean = TRUE)

#gam.check() helps you understand whether you have enough basis functions to model the data.
#Smooths that have significant effects in the diagnostic test (p < 0.05, with asterisks) generally do not have enough basis functions. This indicates non-random patterns in residuals.
gam.check(model.3)

#changed by increasing k=15
set.seed(3)
model.3a <-gam(env$Richness ~
               +s(Temp, k=10)
               +s(Prec, k=15)
               + Date +
                 Collector+
                Habitat
               ,data=env,method="REML")

gam.check(model.3a)

summary(model.3a)

env$exp.Temp <- exp(env$Temp)

model.3e <-gam(env$Richness ~
               +s(Temp)
               +s(Prec)
               +Date+
                Habitat
               ,data=env,method="REML")
cat("Deviance explained =",1-model.3e$deviance/model.3e$null.deviance)

plot(model.3e, all.terms=TRUE, pages = 2)

summary(model.3e) #edf near 1 or 1 is linear

concurvity(model.3e)

# Make the plot with residuals
plot(model.3e, residuals = TRUE)

# Plot all effects
plot(model.3e, residuals = TRUE, pch = 1, cex = 1, all.terms=TRUE, page=1)

# Make another plot adding the intercept value and uncertainty
plot(model.3e, residuals = TRUE, pch = 1, cex = 1, all.terms=TRUE, page=1, shade = TRUE, shade.col = "hotpink", shift = coef(model.3e)[1], seWithMean = TRUE)

gam.check(model.3e)

res.gam <- residuals(model.3a)
plot(res.gam)

quantile(res.gam,probs = c(0,0.025,0.25,0.5,0.75,0.975,1))

qq <- qnorm(seq(0, 1, 1/(length(res.gam)+1)))
qq <- qq[-c(1,length(qq))]
plot(qq, sort(res.gam))
abline(lm(sort(res.gam) ~ qq), col = "red")

# outliers? 
cooksd <- cooks.distance(model.3a)
plot(cooksd, pch = 20, main = "Cook's Distance Plot")

pos <- cooksd < quantile(cooksd,0.95) 
pos[pos == FALSE] <- NA
rich <- ifelse(pos,env$Richness,NA)
write.csv(rich,"outlierRM.csv",row.names = FALSE)

model.3b<- gam(rich ~
               +s(Temp, k=10)
               +s(Prec, k=15)+
                 Date+
                 Habitat
               ,data=env,method="REML")

cooksdb <- cooks.distance(model.3b)
plot(cooksdb, pch = 20, main = "Cook's Distance Plot")
posb <- cooksd < quantile(cooksdb,0.95) 
posb[posb == FALSE] <- NA
richb <- ifelse(posb,env$Richness,NA)

coef(model.3b)

gam.check(model.3b)

summary(model.3b, full=TRUE)

concurvity(model.3b, full=TRUE) #overall

op <- par(mfrow=c(3,2))
plot(model.3b,select=1,scheme=1)
plot(model.3b,select=2,scheme=1)
plot(model.3b,select=3,scheme=1)
plot(model.3b,select=4,scheme=1)
plot(model.3b,select=5,scheme=1)
par(op)

saveRDS(model.3, file = "wio_rich_gam.rds")

# a function to extract predictors of a gam
get_predictors <- function(model.3b) {
  # Check if the model is a fitted gam object
  if (!inherits(model.3b, "gam")) {
    stop("The model must be a fitted gam object")
  }
  
  # Get the formula from the model
  formula <- model.3b$formula
  
  # Get the predictors from the formula
  predictors <- all.vars(formula)
  
  # Return the list of predictors
  return(predictors)
}

#tried glm 5% explained

model.3bc<- glm(rich ~
                Temp
                +Prec
                +Date+Habitat
                ,data=env,family=quasipoisson())

cat("Deviance explained =",1-model.3bc$deviance/model.3bc$null.deviance)

plot(model.3bc, all.terms=TRUE, pages = 2)
 


















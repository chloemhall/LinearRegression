#Script for Davide Warm - lmer ####
#load packages 
library(tidyverse)
library(data.table)
#library(nlme)
library(lme4)
#library(emmeans) #for estimated marginal means. Can also do post-hoc contrasts, but is very intensive for large datasets e.g. here.  
library(report)
#library(stargazer) # for LateX-friendly tables. 
#library(patchwork)
#load data
setwd("/Users/chloe/Documents/R/Davide_lme")
DW.data<- read.csv("lmdataWF.csv", header=F)

#now rename the column headings, to better descriptive uses. 
colnames(DW.data)[1]  <- "Value"    # change column name for x column
colnames(DW.data)[2] <- "n"
colnames(DW.data)[3] <- "Age"
colnames(DW.data)[4] <- "Region"
#view(DW.data)

#### now we have a nice df. ####
### lme
Model1<- lmer(Value ~Age*Region+
                (1|n), 
              na.action=na.exclude, data = DW.data)
summary(Model1) # display summary of results
report(Model1) # display word results
report_table(Model1) # display table results


fixef(Model1)
plot(Model1, type = c("p", "smooth")) # data looks heteroscedastic 

#repeat on logged data
Model2 <- lmer(log10(Value) ~Age*Region+
                          (1|n), 
                        na.action=na.exclude, data = DW.data)
summary(Model2) # display summary of results
report(Model2) # display word results
report_table(Model2) # display table results


fixef(Model2)
plot(Model2, type = c("p", "smooth")) # data looks homoscedastic after log10 transformation 

plot(Model2, sqrt(abs(resid(.))) ~ fitted(.), type = c("p", "smooth"))
qqnorm(resid(Model2)) # q-q plot for normality. Roughly normal; should be ok.  
plot(ranef(Model2)) #?? 
library(car)
vif(Model2) # ??

summary(lm(Model2))

#test for linearity. We want to see random distribution of data,
# no set patterns. 
plot(resid(Model2), DW.data$Value) # hmmmm looks a bit of a linear relationship- unideal.  


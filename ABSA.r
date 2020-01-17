# ABSA in practice

# Call needed libraries

library(data.table)
library(nnet)
library(ROCR)
library(caret)

## Set your working directory - Should contain the 5 .r files described below and the needed input files


## Call functions built to run ABSA

# Internal functions
source("~/AB_SA/R/ReadParseScoary.r") # function called by create_input_mnl.r
source("~/AB_SA/R/DefStrainCategory.r") # function called by create_input_mnl.r

# Functions for ABSA user
source("~/AB_SA/R/CreateInputMNL.r") # function to be used to create mnl inputs (for fitting and prediction)
source("~/AB_SA/R/MNLTrainTest.r") # Training and testing the model
source("~/AB_SA/R/MNLFit.r") # Fitting (full dataset)
source("~/AB_SA/R/MNLPredict.r") # Predicting membership probabilities of strains to attribute 

## Application of ABSA on a dataset (here 'Salmonella' Typhimurium COMPARE French dataset)
setwd("~/AB_SA/data")
CreateInputMNL("~/AB_SA/data/FR_scoary_trait.csv","~/AB_SA/data/gene_presence_absence.Rtab",4)

testedMNL<-MNLTrainTest("mnl_input_0.csv",0.70,100)

final.trained<-MNLFit("mnl_input_0.csv")

predict.unknown<-MNLPredict("predict_sporadic.csv",final.trained)
write.table(file="predicted_sources.csv",predict.unknown,sep=";")
barplot(t(predict.unknown),legend=row.names(t(predict.unknown)),args.legend = list(x='right',bty='n',inset=c(-0.1,0),xpd=TRUE),xlim = c(0,45),cex.names = 0.8,xlab="Environnemental strains",ylab = "Membership probabilities")


## Optimisation for French dataset
AIC<-c()
Accuracy<-matrix(c(0),10,3)
for (ng in 1:10)
{CreateInputMNL("~/AB_SA/data/FR_scoary_trait.csv","~/AB_SA/data/gene_presence_absence.Rtab",ng)
  testedMNL<-MNLTrainTest("mnl_input_0.csv",0.70,10)
  percentiles_accuracy<-testedMNL[[1]]
  Accuracy[ng,1:3]<-percentiles_accuracy[c(1,3,5)]
  final.trained<-MNLFit("mnl_input_0.csv")
  AIC[ng]<-final.trained$AIC
}

# Plot of coef 
library(questionr)
library(GGally)
mnl_input<-c("mnl_input_0.csv")
data<-read.table(mnl_input,sep=",",header = TRUE)
data2<-subset(data,select=-1)
data2$Source<-relevel(data2$Source,ref="Ruminant_FR")
multinomModel_full <- nnet::multinom(Source ~ ., data=data2,trace=FALSE) # multinom Model
ggcoef(multinomModel_full,exponentiate=FALSE,conf.int = FALSE)+facet_grid(~y.level)
tidy(multinomModel_full,exponentiate = TRUE,conf.int = FALSE)
odds.ratio(multinomModel_full)


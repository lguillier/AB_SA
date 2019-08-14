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

## Optimisation for French dataset
AIC<-c()
Accuracy<-matrix(c(0),20,3)
for (ng in 1:20)
{CreateInputMNL("~/AB_SA/data/FR_scoary_trait.csv","~/AB_SA/data/gene_presence_absence.Rtab",ng)
  testedMNL<-MNLTrainTest("mnl_input_0.csv",0.70,10)
  percentiles_accuracy<-testedMNL[[1]]
  Accuracy[ng,1:3]<-percentiles_accuracy[c(1,3,5)]
  final.trained<-MNLFit("mnl_input_0.csv")
  AIC[ng]<-final.trained$AIC
}
  

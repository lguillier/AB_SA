# ABSA in practice

# Call needed libraries

library(data.table)
library(nnet)
library(ROCR)
library(caret)

## Set your working directory - Should contain the 5 .r files described below and the needed input files


## Call functions built to run ABSA

# Internal functions
source("ReadParseScoary.r") # function called by create_input_mnl.r
source("DefStrainCategory.r") # function called by create_input_mnl.r

# Functions for ABSA user
source("CreateInputMNL.r") # function to be used to create mnl inputs (for fitting and prediction)
source("MNLTrainTest.r") # Training and testing the model
source("MNLFit.r") # Fitting (full dataset)
source("MNLPredict.r") # Predicting membership probabilities of strains to attribute 

## Applied ABSA on a dataset (here Salmonella Typhimurium )

CreateInputMNL("DE_scoary_trait.csv","gene_presence_absence.Rtab",3)

testedMNL<-MNLTrainTest("mnl_input_0.csv",0.70,10)

final.trained<-MNLFit("mnl_input_0.csv")

predict.unknown<-MNLPredict("predict_sporadic.csv",final_trained)

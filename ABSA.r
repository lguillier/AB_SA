# ABSA in practice

# Call needed libraries

library(data.table)
library(nnet)
library(ROCR)
library(caret)

## Set your working directory - Should contain the 5 .r files described below and the needed input files


## Call functions built to run ABSA

# Internal functions
source("read_parse_scoary.r") # function called by create_input_mnl.r
source("def_strain_category.r") # function called by create_input_mnl.r

# Functions for ABSA user
source("create_input_mnl.r") # function to be used to create mnl inputs (for fitting and prediction)
source("multinomial_train_test.r") # Training and testing the model
source("multinomial_fit.r") # Fitting (full dataset)
source("multinomial_predict.r") # Predicting membership probabilities of strains to attribute 

## Applied ABSA on DE dataset

create_input_mnl("DE_scoary_trait.csv","gene_presence_absence.Rtab",3)

multinomial_train_test("mnl_input_0.csv",0.70,10)

final_trained<-multinomial_fit("mnl_input_0.csv")

prediction<-multinomial_predict("predict_sporadic.csv",final_trained)

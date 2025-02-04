AB_SA - Source attribution based on accessory genes
========
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3507204.svg)](https://doi.org/10.5281/zenodo.3507204)

Predicts the source of bacterial strains based on their accessory genes.

## Contents
  * [Introduction](#introduction)
  * [Input data](#input-data)
    * [Pangenome](#pangenome)
    * [Genes enriched](#genes-enriched)
  * [Running AB_SA](#running-ab_sa)
    * [Preparing input files](#preparing-inpu-files)
    * [Training and testing a multinomial logistic model](#training-and-testing)
    * [Fitting the multinomial logistic model](#fitting-MNL)
    * [Predicting the source of other strains](#predict)
    * [Plotting the results](#plotting)
  * [Licence](#licence)
  
## Introduction

The [AB_SA](https://github.com/lguillier/AB_SA) ("Accessory-based source attribution") method aims at attributing the origin of bacterial strains associated to human cases or isolated in environment. 
AB_SA method relies on accessory genes and multinomial logistic model. AB_SA uses as input a list of enriched genes in the sources, it trains and tests the multinomial logistic model and can predict from genomes of unkown strains the membership probabilities to the sources.  


## Input data

### Pangenome
The pangenome of the strains (both strains from the sources and the strains to attribute) included in the study should be determined. The AB_SA method parses the .Rtab produced by [Roary](http://sanger-pathogens.github.io/Roary). 

### Genes enriched
In order to know which genes were enriched in each host groups, [Scoary](https://github.com/AdmiralenOla/Scoary) is used. Scoary takes the gene_presence_absence.csv file from Roary as well as a traits file reporting the source associated to each strains. The following code was applied to determine the genes that are enriched in which sources.

```
scoary -g ./DE_pangenome_default/gene_presence_absence.csv -t ./DE_scoary/DE_scoary_trait.csv -p 1E-2 -c I --no_pairwise --collapse
```

The --no_pairwise flag is used for enrichment. The --collapse flag was used to regroup the genes which have the same pattern of distribution in the sources (a single gene, the first of the groups is then taken into accoiunt by AB_SA method) . The naive p-value is used to show the genes most overrepresented in a specific source. A threshold p-value of 0.01 is used. To reduce or increase the number of genes to be used by AB_SA method threshold user can modify the threshold p-value. Association of sources and genes are then carried out by the multinomial model.

## Running AB_SA

### Preparing input files

A script (https://github.com/lguillier/AB_SA/blob/master/R/CreateInputMNL.r) was developped to prepare i) the input file for multinomial logistic model training/testing/fitting and ii) the input file for predicting the source of unknwon strains. The user has to provide the file name of the scoary trait file and the .Rtab file name of the gene presence/absence output from roary. The user defines the maximal number genes to include per source (maxGenes). The genes are chosen according to their rank (i.e. lower p-values).  

````
CreateInputMNL(traitfile="DE_scoary_trait.csv",roary_Rtab="gene_presence_absence.Rtab",maxGenes=10)
````
The function thus produces:
- a file (default name "mnl_input_n.csv" where n is the number of genes selected) for fitting multinomial logistic model,
- a file (default name "predict_sporadic.csv") used for determining the probability of association  the source of unknown strains (sporadic human strains or environmental strains) 

### Training and testing a multinomial logistic model
 In order to assess the performance of the multinomial model with the enriched genes, a random splitting of data into training set (default value percent_cross= 70% for building a predictive model) and test set (1-percent_cross for evaluating the model) is carried out. The radonm splitting is carried out nboot times (default nboot=100)
 
````
testedMNL<-MNLTrainTest("mnl_input_0.csv",percent_cross=0.70,nboot=100)
````
The MNLTrainTest() function returns the accuracy values of the model according to the different training/testing sets. More precisely the quantiles of bootstrap  values of global accuracy, the median of balanced accuracies (one per source) for both accounting to unbalanced dataset and for exploring model performance for each source. 

````
# Explore the outputs:
#[[1]] quantiles describing bootstrap accuracies, 
testedMNL[[1]]

#[[2]] median of balanced accuracies for each source
testedMNL[[2]]

#[[3]] ksdensity of bootstrap accuracies 
testedMNL[[3]]
````

### Fitting the multinomial logistic model

Once the most appropriate model (i.e. after testing different maxGenes and selecting the best model based on accuracy obtained on the test sets) the multinomial logistic model can be fitted on the full input dataset with MNLFit function.
````
final.trained<-MNLFit("mnl_input_0.csv")
````

The AIC value is used as a way to compare different models (in order to prevent overfitting).

````
final.trained$AIC
````
### Predicting the source of other strains

Finally, the source membership probabilities of strains with unknown origin are estimated with MNLPredict.
````
predict.unknown<-MNLPredict("predict_sporadic.csv",final.trained)
write.table(file="predicted_sources.csv",predict.unknown,sep=";")

barplot(t(predict.unknown),legend=row.names(t(predict.unknown)),args.legend = list(x='right',bty='n',inset=c(-0.1,0),xpd=TRUE),xlim = c(0,45),cex.names = 0.7,xlab="Environnemental strains",ylab = "Membership probabilities")
````


## Licence

AB_SA is freely available under a GPLv3 license.

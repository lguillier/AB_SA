AB_SA - Source attribution based on accessory genes
========
Predicts the source of bacterial strains based on their accessory genes.

## Contents
  * [Introduction](#introduction)
  * [Input data](#input-data)
    * [Pangenome](#pangenome)
    * [Genes enriched](#genes-enriched)
    * [Preparing input files](#preparing-inpu-files)
  * [Running AB_SA](#running-ab_sa)
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
In order to know which genes were enriched in each host groups, [Scoary](https://github.com/AdmiralenOla/Scoary) is used. Scoary takes the gene_presence_absence.csv file from Roary as well as a traits file reporting the source associated to each strains. The following code was applied  to determine the genes that are enriched in which sources.

```
scoary -g ./DE_pangenome_default/gene_presence_absence.csv -t ./DE_scoary/DE_scoary_trait.csv -p 1E-2 -c I --no_pairwise --collapse
```

the --no_pairwise flag is used for enrichment. The --collapse flag was used . The naive p-value is used to show the genes most overrepresented in a specific source. A threshold p-value of 0.01 is used. To reduce or increase the number of genes to be used by AB_SA method threshold as well as . Association of sources and genes are then carried out by the multinomial model.

### Preparing input files
A script (AB_SA/R/CreateInputMNL.r) was developped to prepare i) the input file for multinomial logistic model training/testing/fitting and ii) the input file for predicting the source of unknwon strains 

## Running AB_SA

### Training and testing a multinomial logistic model

### Fitting the multinomial logistic model

### Predicting the source of other strains

### Ploting the results

## Licence

AB_SA is freely available under a GPLv3 license.

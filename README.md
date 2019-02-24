AB_SA - Source attribution based on accessory genes
========
Predicts the source of bacterial strains based on their accessory genes.

## Contents
  * [Introduction](#introduction)
  * [Input data](#input-data)
    * [Pangenome](#pangenome)
    * [Genes enriched](#genes-enriched)

## Introduction

The AB_SA ("Accessory-based source attribution") method aims at attributing the origin of bacterial strains associated to human cases or isolated in environment. 
AB_SA method relies on accessory genes and multinomial logistic model. AB_SA uses as input a list of enriched genes in the sources, it trains and tests the multinomial logistic model and can predict  genomes of strains with unknown sources.  


## Input data

### Pangenome
The pangenome of the strains (both strains from the sources and the strains to attribute) included in the study should be determined. The AB_SA method parses the .Rtab produced by [Roary](http://sanger-pathogens.github.io/Roary). 

### Genes enriched
In order to know which genes were enriched in each host groups, [Scoary](https://github.com/AdmiralenOla/Scoary) is used. Scoary takes the gene_presence_absence.csv file from Roary as well as a traits file reporting the source associated to each strains. The following code was applied  to determine the genes that are enriched in which sources.


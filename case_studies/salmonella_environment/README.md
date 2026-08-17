# AB_SA in practice — Salmonella case study

This folder contains the first case study developed with **AB_SA (Attribution Based on Source Association)** and corresponds to the application presented in the original *Microbial Genomics* article describing the approach.

The case study applies AB_SA to a French *Salmonella* Typhimurium dataset from the COMPARE project. It illustrates the main steps of the workflow:

1. creation of the multinomial logistic regression input from Scoary and Roary outputs;
2. assessment of model predictive performance;
3. fitting of the final multinomial model;
4. prediction of the source of strains of unknown origin;
5. optimisation of the number of source-associated genes used in the model.

## Repository structure

The scripts used by AB_SA are located in:

```text
R/
```

The input data required for this case study are located in:

```text
data/salmonella/
```

The present folder contains the scripts specific to this published case study.

## Running the case study

The example below assumes that the R working directory is set to the root of the AB_SA repository.

### Load required libraries

```r
library(data.table)
library(nnet)
library(ROCR)
library(caret)
```

### Load AB_SA functions

```r
# Internal functions
source("R/ReadParseScoary.r")
source("R/DefStrainCategory.r")

# Main AB_SA functions
source("R/CreateInputMNL.r")
source("R/MNLTrainTest.r")
source("R/MNLFit.r")
source("R/MNLPredict.r")
```

### Create the multinomial logistic model input

For a chosen maximum number of enriched genes per source (`maxGenes`):

```r
CreateInputMNL(
  "data/salmonella/FR_scoary_trait.csv",
  "data/salmonella/gene_presence_absence.Rtab",
  maxGenes = 4
)
```

This creates the input file used for multinomial logistic regression.

### Assess model performance

```r
testedMNL <- MNLTrainTest(
  "mnl_input_0.csv",
  percent_cross = 0.70,
  nboot = 100
)
```

The main outputs are:

```r
# Accuracy distribution summary
testedMNL[[1]]

# Median balanced accuracy for each source
testedMNL[[2]]

# Distribution of accuracies
testedMNL[[3]]
```

### Fit the final model

The complete dataset can then be used to fit the final multinomial model:

```r
final.trained <- MNLFit("mnl_input_0.csv")
```

### Predict the source of strains of unknown origin

```r
predict.unknown <- MNLPredict(
  "data/salmonella/predict_sporadic.csv",
  final.trained
)

write.table(
  predict.unknown,
  file = "predicted_sources.csv",
  sep = ";"
)
```

Membership probabilities can for example be visualised using:

```r
barplot(
  t(predict.unknown),
  legend = row.names(t(predict.unknown)),
  args.legend = list(
    x = "right",
    bty = "n",
    inset = c(-0.1, 0),
    xpd = TRUE
  ),
  xlim = c(0, 45),
  cex.names = 0.7,
  xlab = "Environmental strains",
  ylab = "Membership probabilities"
)
```

## Optimisation performed for the published application

For the application presented in the original AB_SA publication, different maximum numbers of enriched genes per source were compared.

```r
AIC <- c()
coefnames <- c()

Accuracy <- matrix(0, 10, 3)
Balanced_accuracies <- matrix(0, 10, 3)

for (ng in 1:5) {

  CreateInputMNL(
    "data/salmonella/FR_scoary_trait.csv",
    "data/salmonella/gene_presence_absence.Rtab",
    ng
  )

  testedMNL <- MNLTrainTest(
    "mnl_input_0.csv",
    0.70,
    100
  )

  percentiles_accuracy <- testedMNL[[1]]
  balanced_accuracies <- testedMNL[[2]]

  Accuracy[ng, 1:3] <- percentiles_accuracy[c(1, 3, 5)]
  Balanced_accuracies[ng, 1:3] <- balanced_accuracies

  final.trained <- MNLFit("mnl_input_0.csv")

  AIC[ng] <- final.trained$AIC
  coefnames <- c(coefnames, final.trained$coefnames)
}
```

The resulting models can be compared using their predictive performance and AIC.

## Exploring model coefficients

```r
library(questionr)
library(GGally)

data <- read.table(
  "mnl_input_0.csv",
  sep = ",",
  header = TRUE
)

data2 <- subset(data, select = -1)
data2$Source <- relevel(data2$Source, ref = "Ruminant_FR")

multinomModel_full <- nnet::multinom(
  Source ~ .,
  data = data2,
  trace = FALSE
)

ggcoef(
  multinomModel_full,
  exponentiate = FALSE,
  conf.int = FALSE
) +
  facet_grid(~y.level)

odds.ratio(multinomModel_full)
```

## Note

This case study corresponds to the original implementation of AB_SA used for the first publication describing the approach. Some AB_SA functions may have subsequently been updated or improved; therefore, numerical results obtained with the current version of the code may differ slightly from those obtained with the original release.

#Code to run the compositional regression models for DTU
library(CompDTUReg)

#Specify the same directory where the gene level files to be used for the analysis were saved in (3)
save_dir <- "/pine/scr/s/k/skvanbur/SQCCDataReproduceOldResBeforeCommonCodeTest/"
direc_to_save <- paste0(save_dir, "GeneLevelFiles/")

#Generate list of all gene level files
GeneFiles <- list.files(direc_to_save, full.names=T)

#Create fake additional predictors just for illustration of how they can be included
pred1 <- c(54,23,45,26,78,33,22,44,55,66)
pred2 <- c(5,2,5,2,5,2,5,2,5,2)

#Creating extraPredictors as data frame then convert to matrix to result in a matrix with each column
  #being a predictor with the names pred1, pred2, as desired
extraPredictors <- as.matrix(data.frame(pred1, pred2))

#CompDTUReg is run separately for each gene, so easiest way to run is using lapply with the startCompDTUReg
  #function
#Condition variable is included automatically, and the returned pvalue will be the pvalue for the significance
  #test of condition (controlling for any extra predictors if included)

CompDTUResultsList <- lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, extraPredictors = extraPredictors)
CompDTUmeResultsList <- lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, extraPredictors = extraPredictors)

#Aggregate results from all genes into one data.table object
CompDTUResults <- rbindlist(CompDTUResultsList)
CompDTUmeResults <- rbindlist(CompDTUmeResultsList)


#Sample to run the compositional regression models for DTU analysis

library(CompDTUReg)

#Specify the same directory where the gene level files to be used for the analysis were saved in file (3)
def_wd <- "/Users/Scott/Documents/Dissertation Data/CompDTURegData/"

#Directory where the gene level files were saved in file (3)SaveNecessaryDatasetsForCompDTUReg.R  (or whatever directory the GeneLevelFiles to be used in the analysis are saved in)
GeneLevelFilesSaveDir <- paste0(def_wd, "GeneLevelFiles/")

#Generate list of all gene level files
GeneFiles <- list.files(GeneLevelFilesSaveDir, full.names = TRUE)

#Create fake additional predictors for illustration of how they can be included
pred1 <- c(54,23,45,26,78,33,22,44,55,66)
pred2 <- c(5,2,5,2,5,2,5,2,5,2)

#Creating extraPredictors as data frame then convert to matrix to result in a matrix with each column
  #being a predictor with the names pred1, pred2, as desired
  #Can leave at the default value of NULL if not needed
extraPredictors <- as.matrix(data.frame(pred1, pred2))



#CompDTUReg is run separately for each gene, so easiest way to run is using lapply with the startCompDTUReg
  #function
#Condition variable is loaded by the function and included automatically, and the returned pvalue will be the pvalue for the significance
  #test of condition (controlling for any extra predictors if included)
#Aggregate results from all genes into one data.table object using rbindlist from data.table
CompDTUResults1 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, extraPredictors = extraPredictors))
CompDTUmeResults1 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, extraPredictors = extraPredictors))




#Now, examples providing custom specified design matrices under the null and alternative hypotheses

#Load the first gene-level file to get the group (cond) information
key <- CompDTUReg:::loadRData(GeneFiles[1], objNameToGet = "key")
cond <- key$Condition

#Specify null and altertive design matrices
NullDesign2 <- model.matrix(~pred1 + pred2)
AltDesign2 <- model.matrix(~pred1 + pred2 + cond)

#These results are the same as the results above and test for DTU (ie the significance of the Group (cond) variable)
CompDTUResults2 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, customHypTest = TRUE, NullDesign = NullDesign2, AltDesign = AltDesign2))
CompDTUmeResults2 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, customHypTest = TRUE, NullDesign = NullDesign2, AltDesign = AltDesign2))



#Now, we could test for the significance of pred2 keeping cond and pred1 in the model

#Specify null and alternative design matrices
NullDesign3 <- model.matrix(~pred1 + cond)
AltDesign3 <- model.matrix(~pred1 + pred2 + cond)


CompDTUResults3 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, customHypTest = TRUE, NullDesign = NullDesign3, AltDesign = AltDesign3))
CompDTUmeResults3 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, customHypTest = TRUE, NullDesign = NullDesign3, AltDesign = AltDesign3))



#Save Inferential Replicates to R with one file per biological sample
#These files will be saved within the outer
library(CompDTUReg)

#SalmonFilesDir is the directory where the Salmon quantification results are saved
#and where the sample specific inferential replicates will be saved
SalmonFilesDir <- "~/res/SQCCDataReproduceOldResBeforeCommonCode/SalmonReproduceResBeforeCommonCodeBootSamps/"
setwd(SalmonFilesDir)

#Code needs to loop over the biological samples in some form, we provide sample
  #code for a slurm array but a loop could be used as well as long as curr_samp and curr_file_loc get assigned properly
#Array val here needs to match the number of biological samples/replicates in the analysis
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))


#Set the number of observations used in the analysis (called nsamp even if they don't correspond to unique samples)
#Should match the number of rows in key from (1)
nsamp <- 10

QuantFiles <- mixedsort(list.files(pattern = ".sf", recursive = TRUE, full.names = TRUE))

#Names of each element in QuantFiles must be set to "Sample1", "Sample2", etc even if they are not unique biological samples
#This is because this is how the code will expect the columns to be named, just like for the key object
#tximport will name the columns in its created Salmon output object with these names and the code will expect them to be there
#in the format "Sample1", "Sample2", etc
names(QuantFiles) <- paste0("Sample", 1:nsamp)

QuantFiles2 <- QuantFiles[mixedsort(names(QuantFiles))]

array_val2 <- array_val %% nsamp
if(array_val2==0){
  array_val2 <- nsamp
}
curr_samp <- names(QuantFiles2)[array_val2]
curr_file_loc <- QuantFiles2[array_val2]



#Set to true if using Gibbs samples, false if using bootstrap samples
GibbsSamps <- FALSE

#Set value for countsFromAbundnace parameter for use with txImport
#Love (2018) (Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3])
  #recommends "scaledTPM" for DTU analysis
  #See tximport for further options
countsFromAbundance <- "scaledTPM"


dir1 <- "~/res/SQCCDataReproduceOldResBeforeCommonCodeTest/"
load(paste0(dir1,"tx2gene.RData"))


#Files will be saved to SalmonFilesDir/BootSamps for bootstrap samples or /GibbsSamps for Gibbs samples
SaveInfRepDataAsRData(curr_samp = curr_samp, curr_file_loc = curr_file_loc, GibbsSamps = GibbsSamps,
                     countsFromAbundance = countsFromAbundance, direc_to_save_res = SalmonFilesDir)

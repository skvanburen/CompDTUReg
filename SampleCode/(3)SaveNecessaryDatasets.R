#Sleep statement can help performance if there alot of different files running at the same time via an array job
#Sys.sleep(sample(1:100,1))

library(DTUCompReg)


#array_value needs to span the number of parts the data is split up into, specified by nparts above
#An alternative to the array job would be looping over each part number, but it could be time consuming especially
  #with a large number of samples so we recommend the array method if possible
  #if using a loop, ensure curr_part_num gets changed and set properly
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
curr_part_num <- array_val

#Directory where the full inferential replicate datasets and gene level files will be saved
#These files can get quite large, so it is a good idea to specify a temporary
#directory with plenty of extra storage space
save_dir <- "/pine/scr/s/k/skvanbur/SQCCDataReproduceOldResBeforeCommonCodeTest/"
if(!dir.exists(save_dir)){
  dir.create(save_dir)
}


#Directory where the gene level files to be used for the analysis will be saved
direc_to_save <- paste0(save_dir, "GeneLevelFiles/")
if(!dir.exists(direc_to_save)){dir.create(direc_to_save, recursive = T)}

#func_loc <-  "~/code/CompFunctions.R"
#source(func_loc)


#Set countsFromAbundance value to match whatever was used previously
countsFromAbundance <- "scaledTPM"

#Directory previous results are loaded from
#Make sure this matches dir1 from (2)SaveInfRepsAsRData
dir1 <- "~/res/SQCCDataReproduceOldResBeforeCommonCodeTest/"

if(countsFromAbundance=="scaledTPM" | countsFromAbundance=="lengthScaledTPM"){
  load(paste0(dir1, "cntGenecntsScaledTPMFiltered.RData"))
  cntGene <- cntGeneFiltered
}else{
  load(paste0(dir1, "cntGeneFiltered.RData"))
  cntGene <- cntGeneFiltered
}

load(paste0(dir1,"tx2gene.RData"))
load(paste0(dir1, "abDatasetsNoOtherGroupsFiltered.RData"))
fullgenenames <- names(abDatasetsFiltered)


#Specify directory where Salmon quantification files are saved
SalmonFilesDir <- "~/res/SQCCDataReproduceOldResBeforeCommonCode/SalmonReproduceResBeforeCommonCodeBootSamps/"

#Load Salmon quantification object imported into R format using tximport in (1)DataProcessing.R
load(paste0(SalmonFilesDir, "SalmonData.RData"))

def_wd1 <- SalmonFilesDir
setwd(def_wd1)


#Set to "Boot" if using bootstrap samples and "Gibbs" if using Gibbs samples
infReps <- "Boot"
if(infReps=="Boot"){
  GibbsSamps <- FALSE
}else{
  GibbsSamps <- TRUE
}

type <- infReps
dirpiece <- infReps

#Specify number of biological replicates/samples
nsamp <- 10

#Specify the number of parts the dataset should be split into.  The number necessary will depend on the number of biological replicates/samples
  #For example, we use 10 parts for the 10 replicate SEQC data and 100 for the 462 sample E-GEUV-1 data
nparts <- 10


#Will save necessary temporary files that contain all inferential replicates for all biological samples/replicates for a specific set of genes
SaveFullinfRepDat()



#Save the within subject covariance matricies (on the ilr scale)
SaveWithinSubjCovMatrices()



#Save the files that will be directly loaded to run DTUCompReg, with each gene having a separate file
  #If a file already exists for that gene, the function will skip resaving that one to save time
SaveGeneLevelFiles(dir1 = dir1, direc_to_save = direc_to_save)




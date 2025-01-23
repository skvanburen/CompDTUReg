#This file saves the necessary gene-level files for use with the CompDTU regression methods
  #If useInferentialReplicates is set to FALSE the gene-level files will only contain data corresponding to the point estimates (output as Y)
  #and if it is TRUE the gene-level files will also contain data corresponding to the inferential replicates (output as YInfRep)
  #Especially if the number of samples is large it will be very difficult to generate datasets for all inferential replicates without some form of
  #computing cluster available since the total amount of data gets quite large
  #In this case you may be forced to skip use of the inferential replicates and use CompDTU instead of CompDTUme by setting useInferentialReplicates to FALSE

#Set to true if bootstrap/Gibbs samples are available and to be incorporated via the CompDTUme model, FALSE if you only want to use the CompDTU model with no bootstrap or Gibbs samples
useInferentialReplicates <- TRUE

library(CompDTUReg)

#Directory previous results are loaded from
#Make sure this matches def_wd from (2)SaveInfRepsAsRData.R
def_wd <- "/Users/Scott/Documents/Dissertation Data/CompDTURegData/"

#Outer directory where the full inferential replicate datasets (if used) and gene level files will be saved
#These files can get quite large, so it is a good idea to specify a directory with plenty of extra storage space
save_dir <- def_wd

if(!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
}

#Specify number of biological replicates/samples
nsamp <- 10

#Specify the number of parts the full dataset including inferential replicates should be split into.  
  #The number necessary will depend on the number of biological replicates/samples
  #For example, we use 10 parts for the 10 replicate SEQC data and 100 for the 462 sample E-GEUV-1 data
  #It is necessary to split the data up to ensure specific files do not get too large for R to use effectively
  #This becomes especially true as the number of samples increases
  #If these files do become too large to deal with effectively, we recommend using point estimates only via CompDTU by setting useInferentialReplicates to FALSE
nparts <- 10

#Directory where the gene level files to be used for the analysis will be saved
GeneLevelFilesSaveDir <- paste0(save_dir, "GeneLevelFiles/")
if(!dir.exists(GeneLevelFilesSaveDir)){dir.create(GeneLevelFilesSaveDir, recursive = T)}


#Set countsFromAbundance value to match whatever was used previously
countsFromAbundance <- "scaledTPM"


if(countsFromAbundance=="scaledTPM" | countsFromAbundance=="lengthScaledTPM"){
  load(paste0(def_wd, "cntGenecntsScaledTPMFiltered.RData"))
  cntGene <- cntGeneFiltered
}else{
  load(paste0(def_wd, "cntGeneFiltered.RData"))
  cntGene <- cntGeneFiltered
}

load(paste0(def_wd,"tx2gene.RData"))
load(paste0(def_wd, "abDatasetsNoOtherGroupsFiltered.RData"))
filteredgenenames <- names(abDatasetsFiltered)


#Set to "Boot" if using bootstrap samples and "Gibbs" if using Gibbs samples
infReps <- "Boot"
if(infReps=="Boot"){
  GibbsSamps <- FALSE
}else{
  GibbsSamps <- TRUE
}

type <- infReps
dirpiece <- infReps



#######################################
#Code below this line "loops" over nparts and saves the gene-level files necessary for the compositional regression analysis
  #We specify this for an array value approach that can parallelize tasks over a computing cluster for use with inferential replicates
  #array_value would need to span the number of parts the data is split up into, specified by nparts above
  #An alternative to the array job would be looping over each part number in 1:nparts, but a loop could be prohibitively time consuming with inferential replicates
  #and a large number of samples so we recommend a procedure that can run each part separately, such as the array approach mentioned here, if inferential replicates are used
  #If looping over "array_val" instead, ensure curr_part_num gets changed and set properly to range from 1 to nparts and that the SaveGeneLevelFiles function is properly run for each part
  #If useInferentialReplicates is FALSE you will still need to loop over curr_part_num over all code below this line but the code will run much faster and a computing cluster is not likely to be needed
  #so we provide a sample loop below that can be easily used if useInferentialReplicates is FALSE
if(useInferentialReplicates==TRUE){
  array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  curr_part_num <- array_val
  
  #Specify directory where Salmon quantification files are saved
  #should match the same argument from (2)
  #Ensure this SalmonFilesDir ends in a / to ensure code compatibility
  SalmonFilesDir <- paste0(def_wd, "ExampleSalmonQuantifications/")
  
  #Load Salmon quantification object imported into R format using tximport in (1)DataProcessing.R
  load(paste0(SalmonFilesDir, "SalmonData.RData"))
  
  setwd(SalmonFilesDir)
  
  #Will save necessary temporary files that contain all inferential replicates for all biological samples/replicates for a specific set of genes
  SaveFullinfRepDat(SalmonFilesDir = SalmonFilesDir, QuantSalmon = QuantSalmon, abDatasetsFiltered = abDatasetsFiltered, save_dir = save_dir, GibbsSamps = GibbsSamps,
                    filteredgenenames = filteredgenenames, cntGene = cntGene, key = key, nparts = nparts, curr_part_num = curr_part_num)

  print("Saving of Full Inferential Replicate datasets is complete")
  #Save the within subject covariance matricies (on the ilr scale)
  SaveWithinSubjCovMatrices(directory = def_wd, save_dir = save_dir, GibbsSamps = GibbsSamps, curr_part_num = curr_part_num, nsamp = nsamp)
  print("Saving of within subject covariance matrices is complete")
  
  
  #Save the files that will be directly loaded to run CompDTUReg, with each gene having a separate file
  #If a file already exists for that gene, the function will skip resaving that one to save time
  SaveGeneLevelFiles(directory = def_wd, GeneLevelFilesSaveDir = GeneLevelFilesSaveDir, curr_part_num = curr_part_num,
                     useInferentialReplicates = useInferentialReplicates, GibbsSamps = GibbsSamps)
  

}else if(useInferentialReplicates==FALSE){
  for(l in 1:nparts){
    print(paste0("Current part is ", l, " out of ", nparts))
    curr_part_num <- l
    SaveGeneLevelFiles(directory = def_wd, GeneLevelFilesSaveDir = GeneLevelFilesSaveDir, curr_part_num = curr_part_num,
                       useInferentialReplicates = useInferentialReplicates, GibbsSamps = GibbsSamps)
  }
}






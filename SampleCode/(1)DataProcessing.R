#Code to process data quantified from Salmon
library(CompDTUReg)

#def_wd is the top level directory where files will be saved
def_wd <- "/Users/Scott/Documents/Dissertation Data/CompDTURegData/"
if(!dir.exists(def_wd)){
  dir.create(def_wd)
}
setwd(def_wd)

#SalmonFilesDir is the directory where the Salmon quantification results have already been saved
SalmonFilesDir <- paste0(def_wd, "ExampleSalmonQuantifications/")

#func_loc <- "~/code/CompFunctions.R"
#source(func_loc)

#Specify location of annotation to use in maketx2gene
#GENCODE annotations can be downloaded from https://www.gencodegenes.org/human/
#We used the annotations for the reference chromosomes only
txdb_loc <- "~/gencode.v27.annotation.gtf.gz"



#Construct a cluster from parallel package for possible use later.  For no parallelization, use makeCluster(1)
clust <- makeCluster(1)

#Build tx2gene dataframe matching transcripts to genes using GENCODE if it doesn't exist
maketx2gene(txdb_loc = txdb_loc)

#Ensure new tx2gene object can be loaded and load it
load("tx2gene.RData")


#Read in Salmon Files and save results (without inferential replicates for now, as trying to save all of them in one file
  #can easily get prohibitively large)
setwd(SalmonFilesDir)

#Set value for countsFromAbundnace parameter for use with txImport
  #Love (2018) (Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3])
  #recommends "scaledTPM" for DTU analysis
  #See tximport for further options
countsFromAbundance <- "scaledTPM"




#List of Salmon quantification files, which end in .sf
QuantFiles <- mixedsort(list.files(pattern = ".sf", recursive = TRUE, full.names = TRUE))

#Names of each element in QuantFiles must be set to "Sample1", "Sample2", etc even if they are not unique biological samples
  #This is because this is how the code will expect the columns to be named
  #tximport will name the columns in its created Salmon output object with these names and the code will expect them to be there
  #in the format "Sample1", "Sample2", etc
names(QuantFiles) <- paste0("Sample", 1:length(QuantFiles))

#Create key matrix that contains matches samples to conditions and identifiers
#Code later will be expecting key to have columns "Sample", "Condition", and "Identifier"
#"Sample" has a unique identifier for the particular sample/replicate
#"Condition" is a factor variable that corresponds to the different groups/conditions to conduct DTU analysis on
#With "Identifier" containing names as "Sample1", "Sample2", etc even if data isn't corresponding to unique biological samples
  #because this is how future code will expect key to be constructed
key <- matrix(c("SRR950078", "A",
                "SRR950080", "A",
                "SRR950082", "A",
                "SRR950084", "A",
                "SRR950086", "A",
                "SRR950079", "B",
                "SRR950081", "B",
                "SRR950083", "B",
                "SRR950085", "B",
                "SRR950087", "B"), nrow = 10, ncol = 2, byrow = TRUE)
key <- as.data.frame(key, stringsAsFactors = FALSE)
key[3] <- paste0("Sample", 1:nrow(key))
colnames(key) <- c("Sample", "Condition", "Identifier")
key$Condition <- relevel(as.factor(key$Condition), ref = 1)



#Use tximport to load in the results that have been pre-quantified by salmon
  #Drop inferential replicates for now, as they will be read in in file (2)
QuantSalmon <- tximport(QuantFiles, type = "salmon", txOut = TRUE, ignoreTxVersion = FALSE,
                        countsFromAbundance = countsFromAbundance, dropInfReps = TRUE)
fulltransnames <- rownames(QuantSalmon$abundance) #transcript names


#Save the tximport object that contains the results of the salmon quantification as an r quantification
save(QuantSalmon, QuantFiles, key, fulltransnames, countsFromAbundance,
      file = paste0(SalmonFilesDir, "SalmonData.RData"))



setwd(def_wd)

#Files will save to def_wd
  #These will be unfiltered lists with the data with each gene as a separate element
sumToGene(QuantSalmon = QuantSalmon, tx2gene = tx2gene, clust = clust, key = key,
          countsFromAbundance = countsFromAbundance)



##################################################################
#Filter Genes using DRIMSeq's Filtering Approach
#See DRIMSeq Paper for more information
##################################################################
#Load dataframe that contains quantification results for each transcript(rows) and sample (column) along with other information
load("abGene.RData")

#Load list of dataframes that separates data by gene, with element in the list being a data frame of expression for that gene
  #Within this data frame, rows are samples and columns are transcripts
load("abDatasets.RData")



if(countsFromAbundance=="scaledTPM" | countsFromAbundance=="lengthScaledTPM"){
  load("cntGenecntsScaledTPM.RData")
  load("cntDatasetsNoOtherGroupscntsScaledTPM.RData")
  #load("failedgibbssampsCountsScaledTPM.RData")
}else if(countsFromAbundance=="no"){
  load("cntGene.RData")
  load("cntDatasetsNoOtherGroups.RData")
  #load("failedgibbssamps.RData")
}



#The filtering values below can be modified to make the filtering more or less strict
#These are the default values used in
  #Love et al (2018) (Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification [version 3])


#Sample size of smallest condition
n.small <- min(table(key$Condition))
n <- nrow(key)

min_samps_feature_expr <- n.small
min_feature_expr <- 10

min_samps_feature_prop <- n.small
min_feature_prop <- 0.10

min_samps_gene_expr <- n
min_gene_expr <- 10

#Set the list of samples to use to determine which genes/transcripts pass filtering based on
  #Set to a character vector in the form of key$Identifier, for example ("Sample1", "Sample5", etc)
sampstouse <- key$Identifier

DRIMSeqFilter(cntGene = cntGene, key = key, min_samps_feature_expr = min_samps_feature_expr, min_feature_expr = min_feature_expr,
              min_samps_feature_prop = min_samps_feature_prop, min_feature_prop = min_feature_prop, min_samps_gene_expr = min_samps_gene_expr,
              min_gene_expr = min_gene_expr, sampstouse = sampstouse)



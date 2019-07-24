# DTUCompReg

## Introduction

<code>DTUCompReg</code> is an R package that fits the compositional regression models for DTU analysis of RNA-seq data described in (Give bioarxiv link).

## Installation
We recommend installing from Github for the latest version of the code:
```r
install.packages("devtools")
devtools::install_github("skvanburen/DTUCompReg")
library(DTUCompReg)
```

## Model
For details of the model, see (Give bioarxiv link).

## Usage  
For examples of how to use the package, see the files contained within SampleCode within the package.  These files provide a pipeline to run the DTUCompReg method starting with a folder of Salmon quantifications.  We recommend you copy these files and modify these files as needed for your analysis.  The files should be run in order, starting with (1).  Specifically:<br>
 <br>
(1)DataProcessing.R gives example code to import the non-inferential data into R and save initial files needed for further processing, including the tx2gene object that links a transcript to a gene<br>
 <br>
(2)SaveInfRepsAsRData.R gives example code to save the inferential replicate data by sample.  This data must be saved by sample because the data otherwise becomes too large to work with easily in R, especially as the number of samples increases.  If no inferential replicates (ie bootstrap/Gibbs) samples are used, file (2) should be skipped.  <br>
 <br>
(3)SaveNecessaryDatasetsForDTUCompReg.R gives example code to save necessary temporary files that include information about inferential replicates and the GeneLevelFiles that contain all data needed to run the DTU compositional regressions. If no inferential replicates (ie bootstrap/Gibbs) samples are used, useInferentialReplicates (3) should be set to FALSE. <br>
 <br>
(4)RunCompositionalRegressions.R gives example code to run the DTU compositional regression analyses, both with and without incorporating the inferential replicates.  This includes the option to incorporate extra predictors other than condition. <br>
 <br>
 File (2) needs to be run separately for each different sample and (3) needs to be run separately for each part that the data is split up into.  The files need to be split up into parts in this way to make the amount of memory required by each part more manageable. <br>  
   <br>
Our code is written assuming slurm based array jobs, but loops could be used instead by modifying the code (see each file for more details).  We recommend arrays if possible, especially if the number of samples is very high.  Sample array jobs for (2) and (3) respectively are provided as (2)SampleArrayJob.sh and (3)SampleArrayJob.sh.  These files can be run using the following commands: <br>
  <br>
module load r <br>
sbatch --array=1-10 (2)SampleArrayJob.sh <br>
  <br>
module load r <br>
sbatch --array=1-10 (3)SampleArrayJob.sh <br>

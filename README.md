# DTUCompReg

## Introduction

<code>DTUCompReg</code> is an R package that fits the compositional regression models for DTU analysis of RNA-seq data.  This includes the procedure that incorporate inferential replicates into the model and the procedure that does not.

## Installation
We recommend installing from Github for the latest version of the code:
```r
install.packages("devtools")
devtools::install_github("skvanburen/DTUCompReg")
library(DTUCompReg)
```

## Model
For details of the model, see (Give bioarxiv link)

## Usage  
For examples of how to use the package, see the files contained within SampleCode within the package.  These files provide a pipeline to run the DTUCompReg method starting with a folder of Salmon quantifications.  Specifically:

(1)DataProcessing.R shows how to import the non-inferential data into R and save initial files needed for further processing, including the tx2gene object that links a transcript to a gene<br>
(2)SaveInfRepsAsRData.R gives example code to save the inferential replicate data by sample.  This data must be saved by sample because the data otherwise becomes too difficult to work with easily, especially as the number of samples increases.<br>
(3)SaveNecessaryDatasets.R saves necessary temporary files that include information about inferential replicates and the GeneLevelFiles that contain all data needed to run the DTU compositional regressions<br>
(4)RunCompositionalRegressions.R <br>

# CompDTUReg

## Introduction

<code>CompDTUReg</code> is an R package that fits the compositional regression models for DTU analysis of RNA-seq data described in "Differential Transcript Usage Analysis Incorporating Quantification Uncertainty Via Compositional Measurement Error Regression Modeling" (available as a preprint at https://www.biorxiv.org/content/10.1101/2020.05.22.111450v1).  For code to reproduce results from the paper, see the repo located at https://github.com/skvanburen/CompDTUPaperCode. For sample data that can be used with the method, see the repo located at https://github.com/skvanburen/CompDTURegSampleData.

## Installation
We recommend installing from Github for the latest version of the code:
```r
install.packages("devtools")
devtools::install_github("skvanburen/CompDTUReg")
```

## Model
For details of the model, see the preprint at https://www.biorxiv.org/content/10.1101/2020.05.22.111450v1.

## Quick Example
This section will give quick instructions on fitting the `CompDTU` and `CompDTUme` results on data that has already been processed using files (1) through (3) in the SampleCode folder.  For a full walkthrough starting from original data quantified from Salmon, see the "Full Example" section below.
<br/><br/>
This code will use the gene level files in contained in the  "GeneLevelFiles.zip" file from the 'CompDTURegSampleData' repo located at https://github.com/skvanburen/CompDTURegSampleData, which has already been processed. This data contains real transcript-level RNA-Seq abundance data for 10 genes from five replicates of two samples.

In general, data will have to be saved separately for each gene because the size of the data can become too large to load into memory in `R` at once for all genes for all samples if inferential replicates are used and the sample size is large.  Use of this data structure also has the advantage of greatly reducing the amount of memory required to conduct a DTU analysis, and the amount of time required to load small files into `R` is negligible.

<br/><br/>We first load the package, specify the path to the GeneLevelFiles directory, and generate a list of all files saved in the directory.  These files include all necessary information for a given gene but do not need to be loaded into memory by the user.

```r
library(CompDTUReg)
library(data.table)

#Specify the same directory where the "GeneLevelFiles.zip" file was extracted
GeneLevelFilesSaveDir <- "/path/to/GeneLevelFiles/"

#Generate list of all gene level files saved in the directory
GeneFiles <- list.files(GeneLevelFilesSaveDir, full.names = TRUE)
```

The `CompDTU` and `CompDTUme` methods are run separately for each gene, and the easiest way to run the methods is to use `lapply` with the `startCompDTUReg` function.  This function will load in the data for the current gene from the files conatined in `GeneFiles`.  The condition variable is loaded and included automatically, such that the returned `p-value` will be the `p-value` for the significance test of condition.  If inferential replicates are used in the analysis, set `runWithME` to `TRUE` and if inferential replicates are not used, set `runWithME` to `FALSE`.

```r
#Aggregate results from all genes into one data.table object using rbindlist from data.table
CompDTUResults1 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE))
CompDTUmeResults1 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE))
```

The resulting output is given below, and gives the `p-value` as well as the `F` statistic and associated degrees of freedom, which differ for each gene depending on how many transcripts from the gene were incorporated in the model.

```r
> CompDTUResults1
               gene_id pval_CompDTU       FStat NumDF DenomDF
 1: ENSG00000000457.13 3.536641e-01    1.314568     3       6
 2: ENSG00000000460.16 1.718228e-02    7.677140     2       7
 3: ENSG00000000971.15 1.415319e-05   87.143678     1       8
 4: ENSG00000001167.14 1.291726e-04   47.111459     1       8
 5: ENSG00000001460.17 1.470795e-05   86.241708     1       8
 6: ENSG00000001461.16 8.432400e-07  186.831448     2       7
 7: ENSG00000001497.16 4.691125e-11 7196.443058     3       6
 8: ENSG00000001617.11 1.463706e-01    2.610781     3       6
 9: ENSG00000001626.14 2.391880e-03   17.150811     3       6
10: ENSG00000002016.17 9.985113e-06   90.434329     2       7
> CompDTUmeResults1
               gene_id pval_CompDTUme       FStat NumDF DenomDF
 1: ENSG00000000457.13   3.351576e-01    1.384456     3       6
 2: ENSG00000000460.16   9.776200e-02    3.301284     2       7
 3: ENSG00000000971.15   1.096450e-05   93.357461     1       8
 4: ENSG00000001167.14   8.809514e-05   52.555505     1       8
 5: ENSG00000001460.17   9.572183e-06   96.828330     1       8
 6: ENSG00000001461.16   3.779448e-07  235.880349     2       7
 7: ENSG00000001497.16   3.099099e-11 8263.212391     3       6
 8: ENSG00000001617.11   1.131538e-01    3.059782     3       6
 9: ENSG00000001626.14   2.009514e-03   18.312117     3       6
10: ENSG00000002016.17   1.077359e-05   88.416540     2       7

```

Previous results did not incorporate any additional predictors, and we now denomstrate how to control for additional predictors in the model.  First, create two additional predictors:

```r
pred1 <- c(54,23,45,26,78,33,22,44,55,66)
pred2 <- c(5,2,5,2,5,2,5,2,5,2)
```
<br>
Now, load the first gene-level file to get the group (condition) information for each sample and create the null and alternative design matricies.  The rows of the design matrix must be in the same order as key$Identifier.
<br>

```r

key <- CompDTUReg::loadRData(GeneFiles[1], objNameToGet = "key")
cond <- key$Condition

#The rows must be in the same order as key$Identifier, where key is extracted above
NullDesign2 <- model.matrix(~pred1 + pred2)
AltDesign2 <- model.matrix(~pred1 + pred2 + cond)
```


The results below now test for DTU (ie the significance of the Group (cond) variable) while controlling for `pred1` and `pred2`.

```r
CompDTUResults2 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, customHypTest = TRUE, NullDesign = NullDesign2, AltDesign = AltDesign2))
CompDTUmeResults2 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, customHypTest = TRUE, NullDesign = NullDesign2, AltDesign = AltDesign2))
```

The resulting output is given below:


```r
> CompDTUResults2
               gene_id pval_CompDTU       FStat NumDF DenomDF
 1: ENSG00000000457.13 4.462324e-01    1.100050     3       4
 2: ENSG00000000460.16 2.412928e-02    8.589763     2       5
 3: ENSG00000000971.15 7.328083e-05   92.064737     1       6
 4: ENSG00000001167.14 2.837914e-04   56.735082     1       6
 5: ENSG00000001460.17 6.325517e-05   96.954431     1       6
 6: ENSG00000001461.16 3.595979e-06  373.868378     2       5
 7: ENSG00000001497.16 5.486524e-08 7792.988637     3       4
 8: ENSG00000001617.11 1.601098e-01    2.970835     3       4
 9: ENSG00000001626.14 1.319325e-02   14.330870     3       4
10: ENSG00000002016.17 1.561218e-04   80.782661     2       5
> CompDTUmeResults2
               gene_id pval_CompDTUme       FStat NumDF DenomDF
 1: ENSG00000000457.13   3.621268e-01    1.413293     3       4
 2: ENSG00000000460.16   8.251591e-02    4.281487     2       5
 3: ENSG00000000971.15   8.148739e-05   88.682840     1       6
 4: ENSG00000001167.14   2.294471e-04   61.282587     1       6
 5: ENSG00000001460.17   4.293334e-05  111.043312     1       6
 6: ENSG00000001461.16   2.590284e-06  426.642035     2       5
 7: ENSG00000001497.16   3.519355e-08 9730.572883     3       4
 8: ENSG00000001617.11   5.940672e-02    5.915228     3       4
 9: ENSG00000001626.14   1.044313e-02   16.302640     3       4
10: ENSG00000002016.17   1.544114e-04   81.150454     2       5
```


Now, we test for the significance of pred2 controlling for cond and pred1 in the model.  The rows of the design matrix again must be in the same order as key$Identifier, where key is extracted above.
```r
#
NullDesign3 <- model.matrix(~pred1 + cond)
AltDesign3 <- model.matrix(~pred1 + pred2 + cond)

#Run results with the new design matrices
CompDTUResults3 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = FALSE, customHypTest = TRUE, NullDesign = NullDesign3, AltDesign = AltDesign3))
CompDTUmeResults3 <- rbindlist(lapply(GeneFiles, startCompDTUReg, runWithME = TRUE, customHypTest = TRUE, NullDesign = NullDesign3, AltDesign = AltDesign3))
```
Now, the `p-values` correspond to testing for the significance of `pred2` controlling for `cond` and `pred1` in the model.

```r
> CompDTUResults3
               gene_id pval_CompDTU      FStat NumDF DenomDF
 1: ENSG00000000457.13   0.01851216 11.8527488     3       4
 2: ENSG00000000460.16   0.38002707  1.1814101     2       5
 3: ENSG00000000971.15   0.36049530  0.9796941     1       6
 4: ENSG00000001167.14   0.16694083  2.4721324     1       6
 5: ENSG00000001460.17   0.14852338  2.7468248     1       6
 6: ENSG00000001461.16   0.31235118  1.4818303     2       5
 7: ENSG00000001497.16   0.20430431  2.4408919     3       4
 8: ENSG00000001617.11   0.18005943  2.7077321     3       4
 9: ENSG00000001626.14   0.69946907  0.5048401     3       4
10: ENSG00000002016.17   0.57883465  0.6111303     2       5
> CompDTUmeResults3
               gene_id pval_CompDTUme     FStat NumDF DenomDF
 1: ENSG00000000457.13     0.03026367 8.9257441     3       4
 2: ENSG00000000460.16     0.34024339 1.3479031     2       5
 3: ENSG00000000971.15     0.47925429 0.5688665     1       6
 4: ENSG00000001167.14     0.18027793 2.2985731     1       6
 5: ENSG00000001460.17     0.13756737 2.9340220     1       6
 6: ENSG00000001461.16     0.34350545 1.3332448     2       5
 7: ENSG00000001497.16     0.14645040 3.1809141     3       4
 8: ENSG00000001617.11     0.07974494 4.8860872     3       4
 9: ENSG00000001626.14     0.63231157 0.6315158     3       4
10: ENSG00000002016.17     0.56659139 0.6378488     2       5
```

## Full Example
For a full example of how to use the package starting from quantified RNA-Seq data from `Salmon`, see the files contained in the SampleCode folder within the package.    It is recommended you copy these files and modify these files as needed for your analysis.  The files should be run in order, starting with (1).  For example Salmon quantification data for 10 replicates that can be used to test the full pipeline, see the folder 'ExampleSalmonQuantifications' from the 'CompDTURegSampleData' repo located at https://github.com/skvanburen/CompDTURegSampleData.   Specifically:<br>
 <br>
`(1)DataProcessing.R` gives example code to import the non-inferential replicate data into R and save initial files needed for further processing, including the tx2gene object that links a transcript to a gene as well as the datasets after filtering using DRIMSeq's filtering approaches.  Please note that transcript or gene names with special characters (such as :) may cause errors in certain functions and thus should be renamed before using the package.<br>
 <br>
`(2)SaveInfRepsAsRData.R` gives example code to save inferential replicates (ie bootstrap/Gibbs samples) by sample.  This data must be saved separately for each sample because otherwise the data can become too large to work with directly in R, especially as the number of biological samples or inferential replicates becomes very high.  If no inferential replicates  are used, file (2) should be skipped.  <br>
 <br>
`(3)SaveNecessaryDatasetsForCompDTUReg.R` gives example code to save necessary temporary files that include information about inferential replicates and the GeneLevelFiles that contain all data needed to run the DTU compositional regressions. If no inferential replicates (ie bootstrap/Gibbs samples) are used, useInferentialReplicates in file (3) should be set to FALSE. <br>
 <br>
`(4)RunCompositionalRegressions.R` gives example code to run the DTU compositional regression analyses, both with and without incorporating the inferential replicates.  This includes the option to incorporate extra predictors other than condition membership . <br>
 <br>
 File (2) needs to be run separately for each different sample and (3) needs to be run separately for each part that the data is split up into.  The files need to be split up into parts in this way to make the amount of memory required by each part more manageable. <br>  
   <br>
Our sample code is written assuming slurm based array jobs, but loops could be used instead by modifying the code (see each file for more details).  We recommend arrays if possible, especially if the number of samples is very high.  Sample array jobs for (2) and (3) respectively are provided as (2)SampleArrayJob.sh and (3)SampleArrayJob.sh.  These files could be run using the following commands: <br>
  <br>
```
module load r
#Update 10 to be the number of biological samples
sbatch --array=1-10 (2)SampleArrayJob.sh
```
```
module load r
#Update 10 to be the number of parts the data is split up into
sbatch --array=1-10 (3)SampleArrayJob.sh
```

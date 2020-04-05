#This function reads in the major transcript information from abDatasets and adds it into the temporary
#'Add Major Trans Info to Existing abundance and count dataframes
#'
#'
#' \code{addMajorTrans} adds a binary indicator of the major trans for a given gene
#'
#' \code{addMajorTrans} adds a binary indicator variable that takes a value of 1 if
#' that transcript is the major transcript for its gene (and 0 otherwise).  The major transcript
#' is the transcript that has the highest average expression proportion across all samples
#' for the given gene.  The major transcripts are computed in \code{\link{generateData}} and input
#' to this function via the abDatasets file.
#' @inheritParams prepareData
#' @inheritParams generateData
#' @param genestouse A list of genenames that are in use for the compositional analysis.
#' @param abGeneTempF data.frame of abundance information that will have the major transcript info added to it
#' @param cntGeneTempF data.frame of count information that will have the major transcript info added to it. Usually also contains length information.
#' @param abDatasets list of length genestouse containing genewise information per sample.  Also contains the major transcript information as an attribute.  Generally comes from the output of \code{\link{generateData}}.
#'
#' @return a list containing abGene and cntGene dataframes with the major transcript information added
addMajorTrans <- function(genestouse, abGeneTempF, cntGeneTempF, abDatasets, CompMI = FALSE){
  MajorTrans <- data.frame(genestouse, NA)
  colnames(MajorTrans) <- c("gene_id", "MajorTransName")
  MajorTrans$MajorTransName <- lapply(genestouse, function(x) {attr(abDatasets[[x]], "MajorTrans")})

  if(CompMI == TRUE){
    abGene <- merge(abGeneTempF, MajorTrans, by = "gene_id", all = FALSE)
    abGene$MajorTrans <- as.numeric(abGene$MajorTransName==abGene$tx_id)
    rownames(abGene) <- abGene$tx_id

    cntGene <- merge(cntGeneTempF, MajorTrans, by = "gene_id", all = FALSE)
    cntGene$MajorTrans <- as.numeric(cntGene$MajorTransName==cntGene$tx_id)
    rownames(cntGene) <- cntGene$tx_id
  }else{
    abGene <- merge(abGeneTempF, MajorTrans, by = "gene_id", all = TRUE)
    abGene$MajorTrans <- as.numeric(abGene$MajorTransName==abGene$tx_id)
    rownames(abGene) <- abGene$tx_id

    cntGene <- merge(cntGeneTempF, MajorTrans, by = "gene_id", all = TRUE)
    cntGene$MajorTrans <- as.numeric(cntGene$MajorTransName==cntGene$tx_id)
    rownames(cntGene) <- cntGene$tx_id
  }


  abGene <-  abGene[order(abGene$gene_id, abGene$tx_id),]
  cntGene <- cntGene[order(cntGene$gene_id, cntGene$tx_id),]
  return(list(abGene = abGene, cntGene = cntGene))
}


#'Create a tx2gene file from a .gff or .gtf file that maps each transcript to its respective gene
#'
#' \code{maketx2gene} creates a tx2gene file from a .gff or .gtf file that maps each transcript to its respective gene
#'
#' @param txdb_loc is the file location of the .gtf or .gff annotation
#' @param save_loc An optional directory location to save the tx2gene data.frame (will save to current working directory otherwise)
#'
#' @return a data.frame with columns gene_id, tx_id (transcript id) and NTrans (the number of transcripts in the current transcript's gene)
#' @export
maketx2gene <- function(txdb_loc, save_loc = NULL){
  temp1 <- GenomicFeatures::makeTxDbFromGFF(txdb_loc)
  temp2 <- S4Vectors::DataFrame(GenomicFeatures::transcripts(temp1, columns = c("tx_id", "tx_name", "gene_id")))
  tx2genetemp <- data.frame(as.character(temp2$tx_name),as.character(temp2$gene_id), stringsAsFactors = FALSE)
  colnames(tx2genetemp) <- c("tx_id", "gene_id")

  #Get number of transcripts per gene
  numtranspergene <- data.frame(table(tx2genetemp$gene_id), stringsAsFactors = FALSE)
  colnames(numtranspergene) <- c("gene_id", "NTrans")

  tx2gene <- merge(tx2genetemp, numtranspergene, by = "gene_id")


  #Save transcript to gene file for potential use later (in a possibly specified save_loc location)
  if(is.null(save_loc)){
    save(tx2gene, file = "tx2gene.RData")
  }else{
    curr_wd <- getwd()
    setwd(save_loc)
    save(tx2gene, file = "tx2gene.RData")
    setwd(curr_wd)
  }

}






#Outer function called by SumToGene.R that saves data output by SaveSalmonDatatoRData (QuantSalmon) into formats needed to run analyses
#' Summarize the non-inferential rep data from Salmon to gene level (see details)
#'
#' @inheritParams prepareData
#' @param QuantSalmon is the Salmon quantification object output using tximport (see file (1)DataProcessing.R in the package's SampleCode folder for example code)
#' @param clust An optional clust object of class parallel to parallelize within this function.  See \code{\link{makeCluster}} for more information.
#' @param countsFromAbundance character corresponding to the countsFromAbundance parameter used when importing the data with \code{\link{tximport}}.  Possible values are \code{"no"}, \code{"scaledTPM"}, or \code{"lengthScaledTPM"}.
#' @param GenAllGroupCombos is a TRUE/FALSE indicator for generating all possible condition combinations from key$Condition.  Only ever needed for certain power analyses, will almost always be set to FALSE.
#'
#' @return \code{sumToGene} saves initial files from the quantification.  These files include lists of gene-specific expression estimates with and with "OtherGroups",
#' which was a filtering alternative we considered in addition to filters built into \emph{DRIMSeq}.  abDatasets correspond to TPM abundances and cntDatasets correspond to counts that may be scaled relative to TPMs
#' if \code{countsFromAdundance} is either "scaledTPM" or "lengthScaledTPM".
#'
#' abGene and cntGene contain the TPM and (possibly scaled) counts with one row per transcript respectively.  These also contain additional information
#' that may be useful, including total gene expression (TGE) for each biological sample and total expression added up across different genes, mean and total TGE by condition, relative transcript
#' abundance proportions (RTAs), and information about the major transcript for that gene, which is the most highly expressed transcript for that gene across all samples.  See the file (1)DataProcessing.R in the package's SampleCode folder for example code.
#' @export
sumToGene <- function(QuantSalmon, key, tx2gene, clust = NULL, countsFromAbundance, GenAllGroupCombos = FALSE){
  nsamp <- ncol(QuantSalmon$abundance)
  Grouptemp <- key$Condition
  Group <- stats::relevel(as.factor(Grouptemp), ref = 1)

  countsFromAbundance <- QuantSalmon$countsFromAbundance

  abundance <- data.frame(QuantSalmon$abundance)
  #abundance$tx_id <- rownames(abundance)

  counts <- data.frame(QuantSalmon$counts)
  #counts$tx_id <- rownames(counts)

  lengths <- data.frame(QuantSalmon$length)
  #lengths$tx_id <- rownames(lengths)

  ObsData <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                             Group = Group, clust = clust, nsamp = nsamp, key = key,
                             GenAllGroupCombos = GenAllGroupCombos, useExistingMajorTrans = FALSE)

  #Use abundances (ie TPMs) instead of the counts to generate the compositions for compositional analysis
  abGene <- ObsData$abGene
  cntGene <- ObsData$cntGene

  abDatasets <- ObsData$abDatasets
  cntDatasets <- ObsData$cntDatasets

  AllGroupCombinations <- ObsData$AllGroupCombinations

  fullgenenames <- ObsData$fullgenenames

  nonfullrankab <- ObsData$nonfullrankab
  nonfullrankabnames <- ObsData$nonfullrankabnames
  nonfullrankcnt <- ObsData$nonfullrankcnt
  nonfullrankcntnames <- ObsData$nonfullrankcntnames

  abDatasetsCompTime <- ObsData$abDatasetsCompTime
  cntDatasetsCompTime <- ObsData$cntDatasetsCompTime
  abGenecntGeneCompTime <- ObsData$abGenecntGeneCompTime

  ncomb <- nrow(AllGroupCombinations)

  if(countsFromAbundance =="no"){
    cntGenefil <- "cntGene.RData"
    cntDatasetsfil <- "cntDatasets.RData"
    cntDatasetsNoOtherfil <- "cntDatasetsNoOtherGroups.RData"
  }else if(countsFromAbundance =="scaledTPM" | countsFromAbundance=="lengthScaledTPM"){
    cntGenefil <- "cntGenecntsScaledTPM.RData"
    cntDatasetsfil <- "cntDatasetscntsScaledTPM.RData"
    cntDatasetsNoOtherfil <- "cntDatasetsNoOtherGroupscntsScaledTPM.RData"
  }
  #Save abundance (abGene) and counts (cntGene) with nsamp information for use downstream
  save(abGene, nsamp, key, abGenecntGeneCompTime, file = "abGene.RData")
  save(cntGene, nsamp, countsFromAbundance, key, abGenecntGeneCompTime, file = cntGenefil)

  save(abDatasets, nsamp, key, fullgenenames, Group, nonfullrankab, nonfullrankabnames, abDatasetsCompTime, file = "abDatasets.RData")

  save(cntDatasets, nsamp, key, fullgenenames, countsFromAbundance, Group, nonfullrankcnt,
       nonfullrankcntnames, cntDatasetsCompTime, file = cntDatasetsfil)
  if(GenAllGroupCombos==TRUE){
    save(AllGroupCombinations, ncomb, file = "AllGroupCombinations.RData")
  }


  rm(ObsData)
  rm(abGene)
  rm(cntGene)

  #rm(abDatasets)
  rm(cntDatasets)
  rm(AllGroupCombinations)


  #Now, generate the abDatasets without using OtherGroups and save results
  ObsDataNoOtherGroups <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                                          Group = Group, clust = clust, nsamp = nsamp, key = key,
                                          abCompDatasets = abDatasets, useOtherGroups = FALSE, useExistingMajorTrans = TRUE)

  abDatasetsNoOtherGroups <- ObsDataNoOtherGroups$abDatasets
  cntDatasetsNoOtherGroups <- ObsDataNoOtherGroups$cntDatasets

  abDatasetsNoOtherGroupsCompTime <- ObsDataNoOtherGroups$abDatasetsCompTime
  cntDatasetsNoOtherGroupsCompTime <- ObsDataNoOtherGroups$cntDatasetsCompTime

  #Do not need to resave cntGene or abGene when turning other groups off because those do
  #not use other groups at all and won't change


  fullgenenames <- ObsDataNoOtherGroups$fullgenenames

  nonfullrankabNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankab
  nonfullrankabnamesNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankabnames
  nonfullrankcntNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankcnt
  nonfullrankcntnamesNoOtherGroups <- ObsDataNoOtherGroups$nonfullrankcntnames

  save(abDatasetsNoOtherGroups, nsamp, key, fullgenenames, Group, nonfullrankabNoOtherGroups, nonfullrankabnamesNoOtherGroups, abDatasetsNoOtherGroupsCompTime, file = "abDatasetsNoOtherGroups.RData")
  save(cntDatasetsNoOtherGroups, nsamp, key, fullgenenames, countsFromAbundance, Group, nonfullrankcntNoOtherGroups, nonfullrankcntnamesNoOtherGroups, cntDatasetsNoOtherGroupsCompTime, file = cntDatasetsNoOtherfil)
}


sumToGeneHelper <- function(abundance, counts, lengths, tx2gene, Group, clust, nsamp, key = NULL, abCompDatasets = NULL,
                            useExistingOtherGroups = FALSE,  useOtherGroups = TRUE, useExistingMajorTrans = TRUE,
                            useExistingGenes = FALSE, GenAllGroupCombos = FALSE){

  if(sum(colnames(abundance) != paste0("Sample", 1:nsamp)) != 0){
    stop("Columns must be ordered by RepsNames tx_ids must just be row names and not a column")
  }

  if(is.null(clust)){
    clust <- parallel::makeCluster(1)
  }
  #Get data into initial format needed

  ST1 <- proc.time()
  initialData <- prepareData(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                             nsamp = nsamp, key = key)

  abGenecntGeneCompTimeP1 <- proc.time() - ST1
  abGeneTempF <- initialData$abGeneTempF
  cntGeneTempF <- initialData$cntGeneTempF


  ##################################################################################################################
  #Generate list of data frames to be used in compositional analysis
  ##################################################################################################################

  #Cant use if the gene only has 1 trans (since no isoform switching/differential splicing to test for otherwise)
  #or if gene expression across all samples is 0 for all samples
  #So remove genes with only 1 total transcript or ones with no expression in any trans/sample combination
  CompAbGene <- subset(abGeneTempF, (abGeneTempF$NTrans!=1 & abGeneTempF$SumTGE!=0))
  CompCntGene <- subset(cntGeneTempF, (cntGeneTempF$NTrans!=1 & cntGeneTempF$SumTGE!=0))

  #Gene names, make sure to sort this so the right transcript go with the right genes down stream
  #Genenames should be the same using abundance or count data, confirm this with the line below (should be 0)
  if(sum(sort(unique(CompAbGene$gene_id))!=sort(unique(CompCntGene$gene_id))) !=0){
    stop("Something is wrong the the genenames, they should be the same between count and TPM data but are not")
  }

  if(useExistingGenes==TRUE){
    fullgenenames <- names(abCompDatasets)
  }else{
    fullgenenames <- sort(unique(CompAbGene$gene_id))
  }

  genestouse <- fullgenenames

  if(useExistingOtherGroups==TRUE | useExistingMajorTrans==TRUE){
    abD <- abCompDatasets
  }else{
    abD <- NULL
  }

  ST2 <- proc.time()
  abDatasets <- parallel::parLapply(clust, genestouse, generateData, dat = CompAbGene,
                          nsamp = length(Group), abundance = TRUE, abData = CompAbGene, abCompDatasets = abD,
                          useExistingOtherGroups = useExistingOtherGroups, useOtherGroups = useOtherGroups,
                          useExistingMajorTrans = useExistingMajorTrans)
  # library(plyr)
  #
  # abDatasets <- laply(genestouse, generateData, dat = CompAbGene,
  #                         nsamp = length(Group), abundance = TRUE, abData = CompAbGene, abCompDatasets = abD,
  #                         useExistingOtherGroups = useExistingOtherGroups, useOtherGroups = useOtherGroups,
  #                         useExistingMajorTrans = useExistingMajorTrans, .progress = "text", .inform = TRUE)

  names(abDatasets) <- genestouse
  abDatasetsCompTime <- proc.time() - ST2

  #Use the other groups from abDatasets created just above
  ST3 <- proc.time()
  cntDatasets <- parallel::parLapply(clust, genestouse, generateData, dat = CompCntGene,
                           nsamp = length(Group), abundance = FALSE, abData = CompAbGene,
                           abCompDatasets = abDatasets, useExistingOtherGroups = TRUE,
                           useOtherGroups = useOtherGroups, useExistingMajorTrans = TRUE)

  # cntDatasets <- laply(genestouse, generateData, dat = CompCntGene,
  #                          nsamp = length(Group), abundance = FALSE, abData = CompAbGene,
  #                          abCompDatasets = abDatasets, useExistingOtherGroups = TRUE,
  #                          useOtherGroups = useOtherGroups, useExistingMajorTrans = TRUE,
  #                          .progress = "text", .inform = TRUE)
  names(cntDatasets) <- genestouse
  cntDatasetsCompTime <- proc.time() - ST3
  #Use this laply loop to help with debugging (will say which gene the code failed at)




  #Add the major transcript to the abundance and count dataframes abGene/cntGene
  #Do this before filtering out genes to get abDatasets and cntDatasets since you could always
  #filter those out later if needed
  #Major trans is the trans witin a gene that has highest average TPM measurement across al samples, even for the count data
  #This is to ensure that the major trans for a gene is the same between the TPM and count measurements
  #This becomes relevant for the power analyses, when the counts of the major transcripts aer modified
  ST4 <- proc.time()
  Temp1 <- addMajorTrans(genestouse = genestouse, abGeneTempF = abGeneTempF, cntGeneTempF = cntGeneTempF, abDatasets = abDatasets)
  abGene <- Temp1$abGene
  cntGene <- Temp1$cntGene
  abGenecntGeneCompTimeP2 <- proc.time() - ST4

  abGenecntGeneCompTime <- abGenecntGeneCompTimeP1 + abGenecntGeneCompTimeP2

  #Generate all possible group arrangements to be able to easily rerun code on each group arrangement
  #This is only possible if there aren't too many samples and is only working for the 10 SQCCData samples for now
  #So by default this is turned of
  if(GenAllGroupCombos==TRUE){
    nsamp <- length(Group)
    numCond1 <- sum(Group==levels(Group)[1])

    #Choose elements corresponding to level 1
    Combs <- gtools::combinations(nsamp, numCond1)
    ncomb <- nrow(Combs)
    #Construct all complete Group combinations, each row is a possible arrangement of group
    AllGroupCombinations <- matrix(NA, nrow = ncomb, ncol = nsamp)
    for (i in 1:ncomb){
      for (j in Combs[i,]){
        AllGroupCombinations[i,j] <- "A"
      }
      for(k in 1:nsamp)
        if(is.na(AllGroupCombinations[i,k])){
          AllGroupCombinations[i,k] <- "B"
        }
    }
  }else{
    AllGroupCombinations <- NULL
  }

  #Keep track of how many genes within abDatasets/cntDatasets are not full rank/ which genes specifically are not full rank
  nonfullrankab <- 0
  nonfullrankabnames <- c()

  nonfullrankcnt <- 0
  nonfullrankcntnames <- c()

  for (i in 1:length(genestouse)){
    if(is.null(attr(abDatasets[[i]], "FullRank"))){
      next
    }
    if(attr(abDatasets[[i]], "FullRank")==FALSE){
      nonfullrankab <- nonfullrankab + 1
      nonfullrankabnames <- c(nonfullrankabnames, names(abDatasets)[i])
    }

    #If there is an error, I think it is at this line
    if(attr(cntDatasets[[i]]$Counts, "FullRank")==FALSE){
      nonfullrankcnt <- nonfullrankcnt + 1
      nonfullrankcntnames <- c(nonfullrankcntnames, names(cntDatasets)[i])
    }
  }


  #Quick check to make sure the majortrans are the same between abundances or counts
  ndiff <- 0
  diff_i <- c()
  for (i in 1:length(genestouse)){
    if(is.null(attr(abDatasets[[i]], "MajorTrans"))){
      next
    }
    if(attr(abDatasets[[i]], "MajorTrans")!=attr(cntDatasets[[i]]$Counts, "MajorTrans")){
      ndiff <- ndiff + 1
      diff_i <- c(diff_i, i)

    }
  }
  if(ndiff!=0){
    stop("Something is wrong with the major trans, it should be the same between TPMs and counts but is not.")
  }

  return(list(abGene = abGene, cntGene = cntGene, abDatasets = abDatasets,
              cntDatasets = cntDatasets, AllGroupCombinations = AllGroupCombinations,
              fullgenenames = fullgenenames, nonfullrankab = nonfullrankab,
              nonfullrankabnames = nonfullrankabnames, nonfullrankcnt = nonfullrankcnt,
              nonfullrankcntnames = nonfullrankcntnames, abDatasetsCompTime = abDatasetsCompTime,
              cntDatasetsCompTime = cntDatasetsCompTime, abGenecntGeneCompTime = abGenecntGeneCompTime))

}



#Prepare raw count/TPM data for use with generateData and later functions
#abundance, counts, and lengths should be data.frames with rows being transcript level values and
#columns corresponding to each sample expression (names should only be Sample1, Sample2, etc not Sample1TPM, etc)
#Samps argument gives sample names
#Only give columns tx_id, Sample1Cnt(or TPM) (or just Sample1, Sample2), etc - nothing else including gene_id
#If key is not specified, the sum of the total gene expression (sumTGE) cannot be calculated for each condition (neither can meanTGE)

#' Prepare raw adbundance/count/length data for later analysis
#'
#' \code{prepareData} generates abundance and count data to be used later, notably in \code{\link{generateData}}
#'
#' @inheritParams generateData
#' @param abundance is a dataframe with nsamp+1 columns, with names Sample1, Sample2, etc and a column for tx_id (that often comes from the rownames).  Rows are transcript level quantification estimates.  Column names should not include "TPM".
#' @param counts is a dataframe with nsamp+1 columns, with names Sample1, Sample2, etc and a column for tx_id (that often comes from the rownames).  Rows are transcript level quantification estimates. Column names should not include "Cnt".
#' @param lengths is a dataframe with nsamp+1 columns, with names Sample1, Sample2, etc and a column for tx_id (that often comes from the rownames).  Rows are transcript level effective length information. Column names should not include "Length".
#' @param tx2gene is a dataframe that matches transcripts to genes. Can be created by \code{\link{maketx2gene}}.
#' @param nsamp is the number of biological samples/replicates used in the analysis
#' @param key is a data.frame with columns "Sample" (corresponding to the unique biological identifier for the analysis), "Condition" (giving the condition/treatment effect variables for the data),
#'  and "Identifier", which should be named "Sample1", "Sample2", ... up to the number of rows of key.  This "Identifier" needs to be created like this even if
#'  the observations don't correspond to unique biological samples.
#' @param samps is an optional vector containing the sample names.  Need to specify this if sample names are not just paste0("Sample", 1:nsamp) without any missing.
#'
#' @return list of length 2 with the first element being the abundance data (abGeneTempF) and the second being the count data (cntGeneTempF) for use with \code{\link{generateData}}
#' @export
prepareData <- function(abundance, counts, lengths, tx2gene, nsamp, key = NULL, infReps = "none", samps = NULL) {
  abundance <- data.frame(abundance)

  if(is.null(abundance$tx_id)){
    abundance$tx_id <- rownames(abundance)
  }

  counts <- data.frame(counts)
  if(is.null(counts$tx_id)){
    counts$tx_id <- rownames(counts)
  }


  lengths <- data.frame(lengths)
  if(is.null(lengths$tx_id)){
    lengths$tx_id <- rownames(lengths)
  }

  if(is.null(abundance$NTrans)){
    abGeneTemp1 <- merge(tx2gene, abundance, by = "tx_id")
  }else{
    abGeneTemp1 <- abundance
  }


  if(is.null(counts$NTrans)){
    cntGeneTemp1 <- merge(tx2gene, counts, by = "tx_id")
  }else{
    cntGeneTemp1 <- counts
  }

  abGeneTemp2 <- as.data.frame(abGeneTemp1[order(abGeneTemp1$gene_id, abGeneTemp1$tx_id),])
  cntGeneTemp2 <- as.data.frame(cntGeneTemp1[order(cntGeneTemp1$gene_id, cntGeneTemp1$tx_id),])

  #t3 <- subset(abGeneTemp2, abGeneTemp2$tx_id=="ENST00000161559.10")
  if(infReps == "none"){
    rownames(abGeneTemp2) <- abGeneTemp2$tx_id
    rownames(cntGeneTemp2) <- cntGeneTemp2$tx_id
  }



  #column numbers of the TPM/count info for all samples
  if(is.null(samps)==FALSE){
    abcolnums <- which(colnames(abGeneTemp2) %in% samps)
    cntcolnums <- which(colnames(cntGeneTemp2) %in% samps)

    if(length(abcolnums)==0){
      abcolnums <- which(colnames(abGeneTemp2) %in% paste0(samps, "TPM"))
    }

    if(length(cntcolnums)==0){
      cntcolnums <- which(colnames(cntGeneTemp2) %in% paste0(samps, "Cnt"))
    }

    colnames(abGeneTemp2)[abcolnums] <- paste0(samps, "TPM")
    colnames(cntGeneTemp2)[cntcolnums] <- paste0(samps, "Cnt")
  }else{
    abcolnums <- which(colnames(abGeneTemp2) %in% paste0("Sample", 1:nsamp))
    cntcolnums <- which(colnames(cntGeneTemp2) %in% paste0("Sample", 1:nsamp))

    if(length(abcolnums)==0){
      abcolnums <- which(colnames(abGeneTemp2) %in% paste0("Sample", 1:nsamp, "TPM"))
    }

    if(length(cntcolnums)==0){
      cntcolnums <- which(colnames(cntGeneTemp2) %in% paste0("Sample", 1:nsamp, "Cnt"))
    }

    colnames(abGeneTemp2)[abcolnums] <- paste0("Sample", 1:nsamp, "TPM")
    colnames(cntGeneTemp2)[cntcolnums] <- paste0("Sample", 1:nsamp, "Cnt")
  }


  #Generate total sums for each sample/gene combo by
  #aggregating over each gene (ie total expression for that gene in that sample)
  abgenetotals <- stats::aggregate(abGeneTemp2[,abcolnums], by = list(gene_id = abGeneTemp2$gene_id), FUN = "sum", drop = FALSE)
  cntgenetotals <- stats::aggregate(cntGeneTemp2[,cntcolnums], by = list(gene_id = cntGeneTemp2$gene_id), FUN = "sum", drop = FALSE)

  #TGE short for total gene expression, ie total expression of that gene for that sample across all transcripts

  if(is.null(samps)==FALSE){
    colnames(abgenetotals) <- c("gene_id", paste0(samps, "TGE"))
    colnames(cntgenetotals) <- c("gene_id", paste0(samps, "TGE"))
  }else{
    colnames(abgenetotals) <- c("gene_id", paste0("Sample", 1:nsamp, "TGE"))
    colnames(cntgenetotals) <- c("gene_id", paste0("Sample", 1:nsamp, "TGE"))
  }


  #Calculate total gene expression across all samples because a gene will have to be dropped from analysis if its expression
  # is 0 for all samples - ie, will later filter out genes with SumTGE=0
  abgenetotals$SumTGE <- rowSums(abgenetotals[,-1], na.rm = TRUE)
  cntgenetotals$SumTGE <- rowSums(cntgenetotals[,-1], na.rm = TRUE)

  #Also calculate MeanTGE level (which is juse mean of SumTGE above)
  abgenetotals$MeanTGE <- abgenetotals$SumTGE/nsamp
  cntgenetotals$MeanTGE <- cntgenetotals$SumTGE/nsamp

  #Also add in information on condition specific TGE
  #If key is not specified, the sum of the total gene expression (sumTGE) cannot be calculated for each condition only (neither can meanTGE)
  if(!is.null(key)){
    ncond <- length(unique(key$Condition))
    Group <- stats::relevel(as.factor(key$Condition), ref = 1)
    for(i in 1:ncond){
      assign(paste0("SampsCond", i), subset(key$Identifier, Group==levels(key$Condition)[i]))
      Samps <- get(paste0("SampsCond", i))
      cols <- paste0(Samps, "TGE")
      condSumAb <- rowSums(abgenetotals[,cols])
      condSumCnt <- rowSums(cntgenetotals[,cols])
      val1 <- paste0("SumTGECond", i)
      abgenetotals[,val1] <- condSumAb
      cntgenetotals[,val1] <- condSumCnt

      val2 <- paste0("MeanTGECond", i)
      abgenetotals[,val2] <- condSumAb/length(get(paste0("SampsCond", i)))
      cntgenetotals[,val2] <- condSumCnt/length(get(paste0("SampsCond", i)))

      cols2 <- paste0(Samps, "TPM")
      tempp <- abGeneTemp2[,cols2]
      val3 <- paste0("MeanTPMCond", i)
      abGeneTemp2[,val3] <- rowMeans(tempp)

      cols3 <- paste0(Samps, "Cnt")
      tempp <- cntGeneTemp2[,cols3]
      val4 <- paste0("MeanCntCond", i)
      cntGeneTemp2[,val4] <- rowMeans(tempp)

    }
  }


  #Merge gene/sample totals back in
  abGeneTemp3 <- merge(abGeneTemp2, abgenetotals, by = "gene_id")
  cntGeneTemp3 <- merge(cntGeneTemp2, cntgenetotals, by = "gene_id")


  #Create relative abundances (proportions) for each sample/gene combo
  #For use in the compositional regression analysis and possibly in later filtering
  if(is.null(samps)==FALSE){
    for (i in 1:length(samps)){
      abGeneTemp3[paste0(samps[i], "RTA")] <-  abGeneTemp3[paste0(samps[i], "TPM")]/abGeneTemp3[paste0(samps[i], "TGE")]
      cntGeneTemp3[paste0(samps[i], "RTA")] <-  cntGeneTemp3[paste0(samps[i], "Cnt")]/cntGeneTemp3[paste0(samps[i], "TGE")]
    }
  }else{
    for (i in 1:nsamp) {
      abGeneTemp3[paste0("Sample", i, "RTA")] <-  abGeneTemp3[paste0("Sample", i, "TPM")]/abGeneTemp3[paste0("Sample", i, "TGE")]
      cntGeneTemp3[paste0("Sample", i, "RTA")] <-  cntGeneTemp3[paste0("Sample", i, "Cnt")]/cntGeneTemp3[paste0("Sample", i, "TGE")]
    }
  }


  #extract which columns are the RTA columns
  if(is.null(samps)==FALSE){
    abRTAcols <- which(colnames(abGeneTemp3) %in% paste0(samps, "RTA"))
    cntRTAcols <- which(colnames(cntGeneTemp3) %in% paste0(samps, "RTA"))
  }else{
    abRTAcols <- which(colnames(abGeneTemp3) %in% paste0("Sample", 1:nsamp, "RTA"))
    cntRTAcols <- which(colnames(cntGeneTemp3) %in% paste0("Sample", 1:nsamp, "RTA"))
  }


  abGeneTemp3$meanRTA <- rowMeans(abGeneTemp3[,abRTAcols], na.rm = TRUE)
  cntGeneTemp3$meanRTA <- rowMeans(cntGeneTemp3[,cntRTAcols], na.rm = TRUE)


  #Merge in length information into count data file for possible use later
  #Keep in mind that the length is the effective length which varies per sample
  len <- lengths
  if(!("tx_id" %in% colnames(len))){
    if(is.null(samps)==FALSE){
      colnames(len) <- paste0(samps, "Len")
    }else{
      colnames(len) <- paste0("Sample", 1:nsamp, "Len")
    }

    len$tx_id <- rownames(len)
  }else{
    pos <- which(colnames(len)=="tx_id")
    if(is.null(samps)==FALSE){
      colnames(len)[-pos] <- paste0(samps, "Len")
    }else{
      colnames(len)[-pos] <- paste0("Sample", 1:nsamp, "Len")
    }

  }

  #len$tx_id <- sapply(strsplit(as.character(rownames(len)), "\\."), "[[", 1)
  cntGeneTemp4 <- merge(cntGeneTemp3, len, by = "tx_id")


  #Order/sort dataframe by gene then by trans within a gene
  abGeneTempF <- abGeneTemp3[order(abGeneTemp3$gene_id, abGeneTemp3$tx_id),]
  cntGeneTempF <- cntGeneTemp4[order(cntGeneTemp4$gene_id, cntGeneTemp4$tx_id),]

  if(infReps=="none"){
    rownames(abGeneTempF) <- abGeneTempF$tx_id
    rownames(cntGeneTempF) <- cntGeneTempF$tx_id
  }

  return(list(abGeneTempF = abGeneTempF, cntGeneTempF = cntGeneTempF))
}



#Generate the datasets with all information that will be needed to do downstream analyses
#Need to run this on abundance data first then use those results for count data
#Generally used with an apply loop and not called directly by user
#Agruments:
#x: genename
#dat: the observed counts/TPM data
#n: number of samples
#abundance; (T/F) is dat adundance (TPM) data or not
#abData: results from dat being abundance (abundance=T) (only non-NULL when running on count data)
#abCompDatasets: list of gene-wise dataframes output when dat is abundance (abundance=T)
#useExistingOtherGroups: (T/F) Should the "Other" groups from the observed abundance be used or new ones computed
#useOtherGroups: Should other groups be created (default Yes)
#infReps: One of "none", "Gibbs", and "Boot", corresponding to none, Gibbs samples, and Boot Samples respectively

#' Generate a list of dataframes with TPM data by Gene
#'
#' \code{generateData} generates a list of data by gene that will be needed for downstream DTU compositional analysis
#'
#' \code{generateData} exports a gene-wise list of all data that will be needed for downstream compositional analysis. This can be done using DRIMSeq's filters or with an approach we considered based on "OtherGroups".
#' This includes combining and transcripts that have <5\% RTA across all samples into an ``Other'' category to ensure proper computation can be done downstream.
#' Note that if there is exactly one transcript that has <5\% RTA it is dropped since there are no other transcript with a low RTA to combine it with and we would not want to combine it with a transcript with high RTA.
#' This also computes the MajorTranscript for each gene, which is the transcript with the highest expression level across all samples and stores it as an attribute.
#' The MajorTranscript is always computed based on the abundance data.  For this reason, need to run this on abundance data first then use those results for count data.
#'
#' @import data.table
#' @param x is a genename of interest.  Genes that have <2 transcripts or have total expression across all samples of 0 are filtered out before calling this function since they can never by used for any kind of DTU analysis
#' @param dat is the observed count or TPM data.  Usually this is the output from \code{\link{prepareData}} that has additionally been filtered by excluding genes with 1 transcript of those that have a total expression level of 0.
#' @param nsamp is the number of biological samples/replicates
#' @param abundance is TRUE/FALSE and indicates whether the data is abundance (TPM) or not.  abundance=F means the length information will also be output for each gene.
#' @param abData is the abundance data as would be used in dat.  This argument is only needed if abundance=F and the results are being run on the count data.  This is used to generate the ``offset'', which is not currently used.
#' @param abCompDatasets is the list of dataframes output by running \code{generateData} on the abundance data.  This argument is only non-Null when abundance=F (ie when running on count data) and is only needed to ensure the calculated ``Other'' Transcripts and major transcripts for each gene are the same when count data is input as when TPM data is input.  Other Groups are currently not used.
#' @param useExistingOtherGroups is a TRUE/FALSE indicator.  If true it will use the ExistingOtherGroups from the abCompDatasets, regardless of the RTA for that particular dataset.  Useful to keep the other groups the same to be able to compare results easier.  This is used for the power analysis, when the other transcript groups should be the same regardless of the current data.
#' @param useOtherGroups is a TRUE/FALSE indicator of whetner other groups should be used or not. Default is FALSE
#' @param useExistingMajorTrans is a TRUE/FALSE indicator of whether to load MajorTrans information from the existing input file.  Useful for the power analyses from the paper or when generating the files corresponding to Bootstrap samples.
#' @param infReps is a character variable indicating what kind of inferential replicates (if any) are to be analyzed by the current function call.  Values to be used should be "none", "Boot", and "Gibbs".  Default is "none".
#' @param ninfreps is the number of inferential replicates being used by the current call. Default is NA, corresponding to none.
#' @param samps is an optional vector containing the sample names.  Need to specify this if sample names are not just paste0("Sample", 1:nsamp) without any missing.
#' @param CompMI is a TRUE/FALSE corresponding to whether datasets for the multiple imputation based analysis are being used.  This will add columns for transcripts that may be missing in the inferential replicates that were't missing in the non-inferential replicate data.  Default is FALSE.
#' @return A list with one element per gene containing TPM or count information for each transcript and the other transcript group.
#' Each element is a dataframe with one row per sample and one column per transcript that is not combined into other (if the OtherGroups are used).
#'  If the data is TPM level each element of the list is the TPM values for that gene, broken down by transcipt and "Other".
#'  If the data is at the count level each element of the list has two elements corresponding to Counts and Lengths.
#'  A list of transcripts that make up the Other category can be viewed in the attribute "OtherTrans", as can a full list of transcripts for that gene ("FullTrans") and a list of transcripts that did not contribute to the Other category "NotOtherTrans".
generateData <- function(x, dat, nsamp, abundance, abData, abCompDatasets = NULL, useExistingOtherGroups, useOtherGroups = FALSE, useExistingMajorTrans = TRUE, infReps = "none", ninfreps = NA, samps = NULL, CompMI = FALSE){
  #for (i in 1:length(fullgenenames)) {
  #x <- genestouse[1]
  temp <- subset(dat, dat$gene_id==x)

  #print(paste0("Currently Running generateData for gene", x))

  #If no rows in temp, no data exists for that gene
  #This can occur especially when using the Gibbs samples if no gibbs samples
  #exist for an entire gene that has "regular" observations for some reason
  if(nrow(temp)==0){
    return(NULL)
  }
  #temp <- subset(dat, dat$gene_id=="ENSG00000002919")

  if(infReps=="Gibbs" | infReps=="GibbsThin16" | infReps=="GibbsThin100"){
    rownames(temp) <- paste0(temp$tx_id, "Gibbs", temp$infRepNum)
    tnames <- unique(temp$tx_id)[unique(temp$tx_id) %in% attributes(abCompDatasets[[x]])$FullTrans]
  }else if(infReps=="Boot"){
    rownames(temp) <- paste0(temp$tx_id, "Boot", temp$infRepNum)
    tnames <- unique(temp$tx_id)[unique(temp$tx_id) %in% attributes(abCompDatasets[[x]])$FullTrans]
  }else if(infReps=="none"){
    tnames <- temp$tx_id
  }

  if(abundance==TRUE){
    if(is.null(samps)==FALSE){
      subcols <- which(colnames(temp) %in% paste0(samps, "TPM"))
    }else{
      subcols <- which(colnames(temp) %in% paste0("Sample", 1:nsamp, "TPM"))
    }

  }else{
    if(is.null(samps)==FALSE){
      subcols <- which(colnames(temp) %in% paste0(samps, "Cnt"))
    }else{
      subcols <- which(colnames(temp) %in% paste0("Sample", 1:nsamp, "Cnt"))
    }
  }

  GDinf <- function(x, data, subcols, infReps){
    data2 <- subset(data, data$infRepNum==x)
    data3 <- data2[,subcols, with = FALSE]
    #data3 <- data.table:::`[.data.table`(data2, , subcols, with = F)

    data4 <- data.frame(data3)
    data5 <- t(data4)
    if(infReps=="Gibbs"| infReps=="GibbsThin16" | infReps=="GibbsThin100"){
      rownames(data5) <- paste0(rownames(data5), "Gibbs", x)
    }else if(infReps=="Boot"){
      rownames(data5) <- paste0(rownames(data5), "Boot", x)
    }
    colnames(data5) <- data2$tx_id
    data6 <- data.frame(data5)
    data6$id <- rownames(data6)
    return(data6)
  }

  if(infReps=="none"){

    sub <- temp[,subcols]
    rownames(sub) <- tnames
    d1 <- t(sub)

  }else{
    sub <- temp[,subcols, with = FALSE]
    #sub <- data.table:::`[.data.table`(temp, , subcols, with = F)
    rownames(sub) <- rownames(temp)
    d1t <- data.frame(data.table::rbindlist(apply(as.matrix(1:ninfreps), 1, GDinf, data = temp,
                                      subcols = subcols, infReps = infReps), use.names = TRUE))
    rownames(d1t) <- d1t$id
    d1t$id <- NULL

    #Use the same trans for Gibbs/Boot replicates datasets as is used in the observed datasets with no other groups
    #Otherwise, the number of columns may not match and the analysis will fail
    #Here abCompDatasets is abDatasetsNoOtherGroups
    d1 <- as.matrix(d1t[,colnames(d1t) %in% attributes(abCompDatasets[[x]])$FullTrans], nrow = nsamp)
    colnames(d1) <- colnames(d1t)[colnames(d1t) %in% attributes(abCompDatasets[[x]])$FullTrans]
    rownames(d1) <- rownames(d1t)
  }



  #If d1 is NULL or 0, gene has no valid data and its dataframe is set to NULL
  if(is.null(ncol(d1))){
    d3 <- NULL
    return(d3)
  }

  if(ncol(d1)==0){
    d3 <- NULL
    return(d3)
  }

  #Ignore transcripts with total count across all samples of 0
  #Just drop these transcripts now so they aren't included in with the Other category if it is used
  #These couldn't ever have a "switching" event so we don't really need to include them
  #Use the existing cols if you want to use the existing other groups or if this is constructing
  #the abDatasets for the Gibbs replicates because the columns/number of columns has to match
  #for the analysis to work right
  if(useExistingOtherGroups==TRUE){
    cols <- which(colnames(d1) %in% attributes(abCompDatasets[[x]])$FullTrans)
    #cols <- 1:ncol(d1)
    #names(cols) <- colnames(d1)
  }else{
    cols <- colSums(d1)>0
  }


  #Need to create the return object in this weird way because of the default r behavior or making a
  # matrix with only one column into a numeric vector
  # Need to include drop = FALSE option so dimensions are not dropped in the case of 1 column
  #As dropping the dimensions would break the R code downstream
  d2 <- data.frame(as.matrix(d1[,cols], nrow = nsamp))
  colnames(d2) <- tnames[cols]

  #Choose which transcripts contribute to the "Other" category based on TPM to ensure the Other category includes the same
  # transcripts for counts and abundance data - otherwise the makeup of this Other category would differ between
  # counts and TPM category
  # So need to run the results on the adundance data first and load in those results into the count results
  if(useOtherGroups==TRUE & useExistingOtherGroups==FALSE){
    #Combine all transcripts with RTA across all samples less than 5% into an "Other" category if there is more than 1
    #If there is exactly 1 trans with RTA < 5%, just drop it for now
    fullcolnames <- colnames(d2)

    #Col #s that have RTA < 5% after summing across all samples- not the same as RTA quantity, which is calculted per sample
    rarecols <- which(colSums(d2)/sum(d2) < 0.05)
    rarecolnames <- colnames(d2)[rarecols]

    #Need to generate like this and not just -rarecols because there could be no rarecols
    # and code would break in that case
    notrarecolnames <- colnames(d2)[!(colnames(d2) %in% rarecolnames)]
  }else if(useOtherGroups==TRUE & useExistingOtherGroups==TRUE){
    fullcolnames <- attr(abCompDatasets[[x]], "FullTrans")
    rarecolnames <- attr(abCompDatasets[[x]], "OtherTrans")
    #rarecols <- which(colnames(d2) %in% attr(abCompDatasets[[x]], "OtherTrans"))
    rarecols <- attr(abCompDatasets[[x]], "OtherTrans")
    notrarecolnames <- attr(abCompDatasets[[x]], "NotOtherTrans")
  }else if(useOtherGroups==FALSE){
    fullcolnames <- colnames(d2)
    rarecolnames <- NULL
    rarecols <- NULL
    notrarecolnames <- NULL
  }



  if(useOtherGroups==TRUE & useExistingOtherGroups==FALSE){
    #Combine all transcripts with RTA across all samples less than 5% into an "Other" category if there is more than 1
    #If there is exactly 1 trans with RTA < 5%, just drop it for now
    fullcolnames <- colnames(d2)

    #Col #s that have RTA < 5% after summing across all samples- not the same as RTA quantity, which is calculted per sample
    rarecols <- which(colSums(d2)/sum(d2) < 0.05)
    rarecolnames <- colnames(d2)[rarecols]

    #Need to generate like this and not just -rarecols because there could be no rarecols
    # and code would break in that case
    notrarecolnames <- colnames(d2)[!(colnames(d2) %in% rarecolnames)]
  }else if(useOtherGroups==TRUE & useExistingOtherGroups==TRUE){
    fullcolnames <- attr(abCompDatasets[[x]], "FullTrans")
    rarecolnames <- attr(abCompDatasets[[x]], "OtherTrans")
    #rarecols <- which(colnames(d2) %in% attr(abCompDatasets[[x]], "OtherTrans"))
    rarecols <- attr(abCompDatasets[[x]], "OtherTrans")
    notrarecolnames <- attr(abCompDatasets[[x]], "NotOtherTrans")
  }else if(useOtherGroups==FALSE & is.null(abCompDatasets)){
    fullcolnames <- colnames(d2)
    rarecolnames <- NULL
    rarecols <- NULL
    notrarecolnames <- colnames(d2)
  }else if(useOtherGroups==FALSE & !is.null(abCompDatasets)){
    fullcolnames <- attr(abCompDatasets[[x]], "FullTrans")
    rarecolnames <- NULL
    rarecols <- NULL
    notrarecolnames <- NULL
  }

  #If length of tnames is not the same as length(union(rarecolnames, notrarecolnames)),
  #there is a transcript that is not a member of the other groups
  #that is also not in the gibbs replicates.  This will cause a problem when trying
  #to calculate the covariances on the ilr scale because the number of columns
  #won't match
  #So, in this case set all values of this transcript to 0 to be able to get the dataset
  #To generate and have the proper dimensions to match the regular counts/TPMs
  if((infReps!="none" & length(tnames)!=length(fullcolnames)) | CompMI==TRUE){
    misstrans <- fullcolnames[!(fullcolnames %in% tnames)]
    if(length(misstrans)!=0){
      for(l in 1:length(misstrans)){
        if(is.na(misstrans[l])){next}
        d2[,misstrans[l]] <- 0
      }
    }
  }

  # If d2 has no columns now, move on to the next gene
  if(ncol(d2)==0){
    return(NULL)
  }


  #If there is only 1 rare column, for now drop that trans instead of combining it with anything else
  #d3 is the count or abundance data frame that will be output for each gene

  #If no trans with RTA < 5%, just use d2 dataset from above
  if(length(rarecols)==0){
    d3 <- d2
  }

  #If exactly 1 trans with a trans with RTA < 5%, just drop that trans for now to avoid needing to combine with a higher one
  #Could maybe reevaluate this or perhaps add an option for how this could be handled?
  if(length(rarecols)==1){
    d3 <- d2[notrarecolnames]
  }

  #If multiple trans with RTA < 5%, combine them into the "Other" category
  if(length(rarecols) > 1){
    d2$Other <- as.numeric(rowSums(d2[,rarecols]))
    d3 <- d2[c(notrarecolnames, "Other")]
  }

  #As of now, don't do anything in this file if all trans are "rare" (ie all have RTA <5%)
  #just filter that out before doing any analysis on it- this should be very rare
  #This is a problem since isn't really possible to pick a transcript to choose to not be included in the other category
  #And having only an Other category means we can't do isoform switching analysis (and it looses interpretation also)

  #If there are no cols left, all transcripts are "rare" and the gene can't be used accurately as of now
  #So, just return it as null and go on to the next gene
  if(ncol(d3)==0){
    d3 <- NULL
    return(d3)
  }

  #Use the MajorTrans from the observed TPM values to be the major trans both for the counts
  #and for any data based on the Gibbs replicates
  #If there are multiple "major" transcripts (ie the meanRTA for two transcripts across all samples
  #are exactly tied for the max value) just assign the first transcript of those to be the major transcript
  if(useExistingMajorTrans==FALSE){
    #Now, return column number and name of the major transcript (transcript with highest
    # average TPM across) but don't allow the "Other" category to be the major trans
    if(length(notrarecolnames)==1){
      majtrans <- notrarecolnames
    }else{
      colm <- colMeans(d3[,notrarecolnames], na.rm = TRUE)
      majtrans <- names(which(colm==max(colm)))
    }

    if(length(majtrans) >1){
      majtrans <- majtrans[1]
    }

    attr(d3, "MajorTrans") <-   majtrans
  }else{
    #Use the major trans from TPM for the count data so they are the same
    attr(d3, "MajorTrans") <- attr(abCompDatasets[[x]], "MajorTrans")
  }

  if(length(attr(d3, "MajorTrans")) >1){stop("Something is wrong with the major transcript calculation, there are multiple of them when there should only be 1")}


  #Remember that a transcripts membership in the Other group is determined from observed TPMs always, not the counts
  #or from any specific Gibbs replicate that may be currently being used
  attr(d3, "OtherTrans") <- rarecolnames
  attr(d3, "NotOtherTrans") <- notrarecolnames
  attr(d3, "FullTrans") <- fullcolnames


  #calculating rank for the Gibbs/Boot datasets causes computational issues, and is never needed based on the way the analysis is done, so don't bother calculating it
  if(infReps=="none"){
    attr(d3, "FullRank") <- (Matrix::rankMatrix(d3) == ncol(d3))
  }

  ##################################################################################################################
  #Now ,create length data files for use with count data in case they would ever be needed
  #Don't include these with abundance data (since theyve already been normalized)
  #And, don't include if using Gibbs/Boot Samps since that will just serve to take up space and these
  #can be loaded from the regular cntDatasets file
  ##################################################################################################################
  if(abundance==FALSE & infReps=="none"){
    #First, get the lengths

    subcols2 <- which(colnames(temp) %in% paste0("Sample", 1:nsamp, "Len"))
    sub2 <- temp[,subcols2]
    rownames(sub2) <- tnames
    ln <- data.frame(t(sub2))

    #Again, drop transcripts that are always 0- these are the ones in cols from above
    #Need to create the return object in this weird way because of the default r behavior or making a
    # matrix with only one column into a numeric vector
    # This would break code downstream that relies on dim(), etc so need to create data this way
    ln2 <- data.frame(as.matrix(ln[,fullcolnames], nrow = nsamp))

    #As before, if there are no trans with RTA < 5% just use ln result
    if(length(rarecols)==0){
      ln3 <- ln2
      colnames(ln3) <- notrarecolnames
    }

    #With exactly one trans with <5% RTA, drop for now- could change later
    if(length(rarecols)==1){
      ln3 <- data.frame(as.matrix(ln2[,notrarecolnames], nrow = nsamp))
      colnames(ln3) <- notrarecolnames

    }

    #For >=2 RTA, create other category, with the length for that category just
    #being the sum of the offsets for the trans that make up the other category
    if(length(rarecols) > 1){
      ln2$Other <- rowSums(ln2[,rarecols])
      ln3 <- ln2[,c(notrarecolnames, "Other")]
    }

    rownames(ln3) <- paste0("Sample", 1:nsamp, "Len")
    if(useOtherGroups==FALSE){
      colnames(ln3) <- colnames(d3)
    }
    attr(ln3, "OtherTrans") <- rarecolnames
    attr(ln3, "NotOtherTrans") <- notrarecolnames
    attr(ln3, "FullTrans") <- fullcolnames

    ###############################################################################################
    #This code below was creating an offset, which as of now is never used so is being skipped
    #Offset as of now is the log difference between
    #the count value and the the TPM value (ie offset = log(TPM/Count))- check on this
    #To get the offset for the "other" category, for now sum TPM for the categories used in other and
    #sum counts for the other categories and do Offset =log(sum(TPM)/sum(Count)), where
    #the sums are over all transcripts that make up the other category
    ###############################################################################################

    # tempoff <- subset(abData, abData$gene_id==x)
    # #tempoff <- subset(abData, abData$gene_id=="ENSG00000002919")
    #
    # tnames <- tempoff$tx_id
    #
    # subcols3 <- which(colnames(tempoff) %in% paste0("Sample", 1:nsamp, "TPM"))
    # suboff <- tempoff[,subcols3]
    # rownames(suboff) <- tnames
    #
    # off1 <- t(suboff)
    #
    # #Drop columns (transcripts) that always have 0 expression, as is done above with creating the count object
    # off2 <- as.matrix(off1[,fullcolnames], nrow = nsamp)
    # colnames(off2) <- fullcolnames
    # rownames(off2) <- paste0("Sample", 1:nsamp)
    #
    # denom <- as.matrix(d3[,notrarecolnames], nrow = nsamp)
    # colnames(denom) <- notrarecolnames
    # rownames(denom) <- paste0("Sample", 1:nsamp)
    #
    # #Offset for "regular" (not rare transcripts) is TPM value divided by count value
    # OffsetNotOther <- data.frame(as.matrix(off2[,notrarecolnames]/denom, nrow = nsamp))
    # colnames(OffsetNotOther) <- notrarecolnames
    # rownames(OffsetNotOther) <- paste0("Sample", 1:nsamp)
    #
    #
    # #Offset for "Other" transcripts is the sum of TPM for the other trans divided
    # #by the sum of counts for other trans
    # #Note that if there is exactly 1 rare transcript, it is just dropped above instead of trying to combine
    # #it with another one.  So, only create this "Other" group offset if there are 2 or more raretrans
    # if(length(rarecolnames) > 1){
    #   AbOther <- data.frame(as.matrix(rowSums(as.matrix(off2[,rarecols], nrow = nsamp)), nrow = nsamp, ncol = 1))
    #   colnames(AbOther) <- "Other"
    #   rownames(AbOther) <- paste0("Sample", 1:nsamp, "Off")
    #   OffsetOther <- AbOther/(as.matrix(d3$Other, nrow = nsamp))
    #   #OffsetFull <- log(cbind(OffsetNotOther, OffsetOther))
    #   OffsetFull <- cbind(OffsetNotOther, OffsetOther)
    # }else{
    #   #OffsetFull <- log(OffsetNotOther)
    #   OffsetFull <- OffsetNotOther
    # }
    #
    # if(ncol(d3) != ncol(OffsetFull)){
    #   stop(print(paste("ncol of Counts not equal to ncol of Offsets, somethings wrong with gene", x)))
    # }
    #
    # #Return a list with 1st element being the counts, second element being the lengths, third being Offset
    # ret <- list(d3, ln3, OffsetFull)
    # names(ret) <- c("Counts", "Lengths", "Offsets")


    ret <- list(d3, ln3)
    names(ret) <- c("Counts", "Lengths")


    #dataloop[[i]] <- ret
    return(ret)
  }else{
    return(d3)
  }

} #end generateData






#Filter functions using procedure from \emph{DRIMSeq}
#' Filter data using filtering procedure built into \emph{DRIMSeq} via the \code{\link{dmFilter}} function.  Will automatically save
#' the filtered versions of the various datasets described in \code{\link{sumToGene}}
#' @param abGene is the data.frame of abundances (TPMs) for each sample saved by \code{\link{sumToGene}}
#' @param cntGene is the data.frame of counts and lengths for each sample saved by \code{\link{sumToGene}}
#' @inheritParams prepareData
#' @inheritParams sumToGene
#' @param min_samps_feature_expr From \code{\link{dmFilter}} documentation: Minimal number of samples where features (transcripts) should be expressed
#' @param min_feature_expr From \code{\link{dmFilter}} documentation: Minimal feature (transcript) expression.
#' @param min_samps_feature_prop From \code{\link{dmFilter}} documentation: Minimal number of samples where features (transcripts) should be expressed.
#' @param min_feature_prop From \code{\link{dmFilter}} documentation: Minimal proportion for feature (transcript) expression. This value should be between 0 and 1.
#' @param min_samps_gene_expr From \code{\link{dmFilter}} documentation: Minimal number of samples where genes should be expressed.
#' @param min_gene_expr From \code{\link{dmFilter}} documentation: Minimal gene expression.
#' @param sampstouse is a vector of sample names (in the form of "Sample1", "Sample2", etc) to be used in the analysis.
#' This argument should be used if you only want to run a subset of all sample ID's from key$Identifier.
#' @param failedinfRepsamps is an optional parameter that gives names of samples (in the form of "Sample1", "Sample2", etc) that had the infRep sampler fail.
#' This should not be needed, as newer versions of Salmon don't seem to have this issue but is left for backward compatability.
#'
#' @details This function internally calls \code{\link{dmFilter}}.  See the documentation for that function for more information, 
#' including a discussion of setting all filtering parameters to zero to only remove features with zero expression across all samples and genes with only one non-zero feature (since DTU analysis cannot be performed if a gene has only
#' one transcript.  See also the file (1)DataProcessing.R in the package's SampleCode folder for example code.
#'
#' @return This function will save versions of abGene, cntGene, abDatasets, and cntDatasets containing information for only those genes and transcripts that pass filtering with the given input parameters.
#' For more information on the output datasets see \code{\link{sumToGene}}.
#'
#' @export DRIMSeqFilter
DRIMSeqFilter <- function(abGene, cntGene, key, min_samps_feature_expr, min_feature_expr, min_samps_feature_prop,
                          min_feature_prop, min_samps_gene_expr, min_gene_expr, tx2gene, countsFromAbundance, sampstouse = NULL, failedinfRepsamps = NULL){

  temp1 <- PreDRIMSeq(cntGene = cntGene, failedinfRepsamps = failedinfRepsamps, key = key)
  cnts <- temp1$cnts
  samp <- temp1$samp

  if(is.null(sampstouse)){
    samp2 <- samp
  }else{
    samp2 <- subset(samp, samp$sample_id %in% sampstouse)
    samp2$group <- stats::relevel(factor(samp2$group), ref = 1)
  }

  nsamp <- nrow(samp2)
  Group <- samp2$group

  DRIMData <- DRIMSeq::dmDSdata(counts = cnts, samples = samp2)

  DRIMData2 <- DRIMSeq::dmFilter(DRIMData,
                        min_samps_feature_expr=min_samps_feature_expr, min_feature_expr=min_feature_expr,
                        min_samps_feature_prop=min_samps_feature_prop, min_feature_prop=min_feature_prop,
                        min_samps_gene_expr=min_samps_gene_expr, min_gene_expr=min_gene_expr)


  DRIMSeqFilteredData <- DRIMSeq::counts(DRIMData2)
  transcriptstouse <- DRIMSeqFilteredData$feature_id
  fullgenenames_filtered <- unique(DRIMSeqFilteredData$gene_id)

  counts <- cntGene[cntGene$tx_id %in% transcriptstouse, c("tx_id", paste0(key$Identifier, "Cnt"))]
  rownames(counts) <- counts$tx_id
  counts$tx_id <- NULL
  colnames(counts) <- key$Identifier
  abundance <- abGene[abGene$tx_id %in% transcriptstouse, c("tx_id", paste0(key$Identifier, "TPM"))]
  rownames(abundance) <- abundance$tx_id
  abundance$tx_id <- NULL
  colnames(abundance) <- key$Identifier
  lengths <- cntGene[cntGene$tx_id %in% transcriptstouse, c("tx_id", paste0(key$Identifier, "Len"))]
  rownames(lengths) <- lengths$tx_id
  lengths$tx_id <- NULL
  colnames(lengths) <- key$Identifier

  ST1 <- proc.time()
  initialData <- prepareData(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene,
                             nsamp = nsamp, key = key)

  abGenecntGeneCompTimeP1 <- proc.time() - ST1
  abGeneTempF <- initialData$abGeneTempF
  cntGeneTempF <- initialData$cntGeneTempF

  #library(parallel)
  #clust <- makeCluster(1)
  #FilteredDat <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene, Group = Group, clust = NULL, nsamp = length(Group), key = key, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, useExistingMajorTrans = FALSE)

  #abDatasets are only input to extract existing MajorTrans information from
  FilteredDat <- sumToGeneHelper(abundance = abundance, counts = counts, lengths = lengths, tx2gene = tx2gene, Group = Group, clust = NULL, nsamp = length(Group),
                                 key = key, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, useExistingMajorTrans = FALSE, abCompDatasets = NULL)


  abDatasetsFiltered <- FilteredDat$abDatasets
  cntDatasetsFiltered <- FilteredDat$cntDatasets

  abDatasetsFilteredCompTime <- FilteredDat$abDatasetsCompTime
  cntDatasetsFilteredCompTime <- FilteredDat$cntDatasetsCompTime

  abGeneFiltered <- FilteredDat$abGene
  cntGeneFiltered <- FilteredDat$cntGene

  abGeneFilteredCompTime <- FilteredDat$abGeneCompTime
  cntGeneFilteredCompTime <- FilteredDat$cntGeneCompTime

  NTransFilabGene <- data.frame(table(abGeneFiltered$gene_id))
  colnames(NTransFilabGene) <- c("gene_id", "NTransFiltered")

  NTransFilcntGene <- data.frame(table(cntGeneFiltered$gene_id))
  colnames(NTransFilcntGene) <- c("gene_id", "NTransFiltered")

  abGeneFiltered <- merge(abGeneFiltered, NTransFilabGene, by = "gene_id")
  cntGeneFiltered <- merge(cntGeneFiltered, NTransFilcntGene, by = "gene_id")

  rownames(abGeneFiltered) <- abGeneFiltered$tx_id
  rownames(cntGeneFiltered) <- cntGeneFiltered$tx_id

  if(countsFromAbundance=="scaledTPM" | countsFromAbundance=="lengthScaledTPM"){
    save(abGeneFiltered, nsamp, key, file = "abGeneFiltered.RData")
    save(cntGeneFiltered, nsamp, countsFromAbundance, key, file = "cntGenecntsScaledTPMFiltered.RData")

    save(abDatasetsFiltered, nsamp, key, fullgenenames_filtered, Group, abDatasetsFilteredCompTime, file = "abDatasetsNoOtherGroupsFiltered.RData")
    save(cntDatasetsFiltered, nsamp, key, fullgenenames_filtered, countsFromAbundance, Group, cntDatasetsFilteredCompTime, file = "cntDatasetscntsScaledTPMNoOtherGroupsFiltered.RData")
  }else if(countsFromAbundance=="no"){
    save(abGeneFiltered, nsamp, key, abGeneFilteredCompTime, file = "abGeneFiltered.RData")
    save(cntGeneFiltered, nsamp, countsFromAbundance, key, cntGeneFilteredCompTime, file = "cntGeneFiltered.RData")

    save(abDatasetsFiltered, nsamp, key, fullgenenames_filtered, Group, abDatasetsFilteredCompTime, file = "abDatasetsNoOtherGroupsFiltered.RData")
    save(cntDatasetsFiltered, nsamp, key, fullgenenames_filtered, countsFromAbundance, Group, cntDatasetsFilteredCompTime, file = "cntDatasetsNoOtherGroupsFiltered.RData")
  }
}



#Only specify key_full if you are using a random subset of the samples for the analysis
#This is needed to extract the total number of samples, whcih is needed to be able to get all of the
#Possible count columns of interest which will later be matched to the samples used for that specific analysis
PreDRIMSeq <- function(cntGene, key, failedinfRepsamps = NULL){
  if(is.null(failedinfRepsamps)==TRUE){
    failedinfRepsamps <- integer(0)
  }

  # if(is.null(key_full)==FALSE){
  #   nsamp <- nrow(key_full)
  # }
  #Only include genes that have >1 trans, since only these can be used for differential transcript analysis
  cntGene2 <- cntGene[(cntGene$NTrans > 1 & cntGene$SumTGE!=0),]

  #Only using columns corresponding to counts that had gibbs replicates that worked for comparison with Compositional method
  #sampscols <- paste0("Sample", 1:nsamp, "Cnt")
  if(length(failedinfRepsamps)==0){
    sampstouse <- key$Identifier
    sampstousecols <- paste0(sampstouse, "Cnt")
    cnts <- cntGene2[,c("tx_id", "gene_id", sampstousecols)]

    Group <- key$Condition[key$Identifier %in% sampstouse]

  }else{
    failedsamps <- paste0("Sample", failedinfRepsamps)
    failedsampscols <- paste0(failedsamps, "Cnt")

    sampstouse <- key$Identifier[!(key$Identifier %in% failedsamps)]
    sampstousecols <- paste0(sampstouse, "Cnt")
    #sampstousecols <- sampscols[!(sampscols %in% failedsampscols)]
    cntGene3 <- cntGene2[,!(colnames(cntGene2) %in% failedsampscols)]

    #cntcolnums <- which(colnames(cntGene3) %in% sampstousecols)
    cnts <- cntGene3[,c("tx_id", "gene_id", sampstousecols)]

    Group <- key$Condition[key$Identifier %in% sampstouse]

  }

  colnames(cnts)[colnames(cnts) == "tx_id"] <- "feature_id"
  colnames(cnts)[colnames(cnts) %in% sampstousecols] <- sampstouse
  rownames(cnts) <- cnts$feature_id

  grp <- stats::relevel(factor(Group), ref = 1)
  samp <- data.frame(sample_id = sampstouse, group = grp)

  return(list(cnts = cnts, samp = samp))
}



#This function runs for one biological sample at a Time, and thus needs to have the analysis repeated for each biological sample
#Set GibbsSamps to TRUE if the infReps were Gibbs Samples draws and FALSE if they're boostrap sample draws
#' Save the infRep data for all genes across each unique biological sample/replicate.
#' @inheritParams prepareData
#' @inheritParams sumToGene
#' @param curr_samp Gives the current sample number in the form of "Sample1", "Sample2", etc.  Used in conjunction with key$Identifier, so needs to be "Sample1", "Sample2", etc
#' even if the observation does not correspond to a unique biological sample
#' @param curr_file_loc Gives the location of the .sf file for the Salmon quantification of the current Sample specified by the \code{curr_samp} parameter
#' @param GibbsSamps is TRUE if the inferential replicates are Gibbs samples and FALSE if the replicates are bootstrap samples
#' @param direc_to_save_res is an optional directory to save the sample specific dataset containing bootstrap or Gibbs replicates.  If unspecifying it will be saved to the current working directory.
#'
#' @details The count values will be scaled by the TPM values if \code{countsFromAbundance} is \code{scaledTPM} or \code{lengthScaledTPM} See the file (2)SaveInfRepsAsRData.R in the package's SampleCode folder for example code.
#'
#' @return This function saves a file that contains all bootstrap or Gibbs samples for the given biological sample input by the \code{curr_samp} parameter.  Data is saved in
#' the subdirectory BootSamps or GibbsSamps off of the current working directory unless the parameter \code{direc_to_save_res} is specified.
#' @export SaveInfRepDataAsRData
SaveInfRepDataAsRData <- function(curr_samp, tx2gene, curr_file_loc, GibbsSamps = FALSE, countsFromAbundance = "no", direc_to_save_res = NULL){

  startTime <- proc.time()
  QuantSalmon <- tryCatch(tximport::tximport(curr_file_loc, type = "salmon", txOut = TRUE, countsFromAbundance = countsFromAbundance,
                                   ignoreTxVersion = FALSE, dropInfReps = F), error=function(e){})
  if(is.null(QuantSalmon)){
    stop("The infRep sampler seems to have failed for this sample, so no Gibbs output will be able to be saved")
  }


  t1 <- data.frame(QuantSalmon$infReps)

  #This could be nboot or ngibbs depending on whether gibbs samps or boot samps are used- just call it ninfreps here regardless
  ninfreps <- ncol(t1)
  colnames(t1) <- 1:ninfreps
  rownames(t1) <- rownames(QuantSalmon$abundance)

  t1$tx_id <- rownames(t1)

  #Add in gene and transcript information
  t2 <- merge(t1, tx2gene, by = "tx_id")
  t3 <- t2[order(t2$gene_id, t2$tx_id),]
  rownames(t3) <- t3$tx_id

  attr(t3, "Identifier") <- curr_samp

  #Want to move the tx_id column to be the last column- this is the easiest way to do this since the tx info is in the row names
  txcol <- which(colnames(t3) %in% "tx_id")
  t4 <- t3[,-txcol]
  t4$tx_id <- rownames(t3)

  #Extract length information and order it by gene/transcript
  l1 <- data.frame(QuantSalmon$length)
  colnames(l1) <- paste0(curr_samp, "Len")
  l1$tx_id <- rownames(l1)
  l2 <- merge(l1, tx2gene, by = "tx_id")

  l3 <- l2[order(l2$gene_id, l2$tx_id),]
  rownames(l3) <- l3$tx_id


  #Stack all gibbs/boot samples for this biological sample on top of each other
  outp <- data.table::rbindlist(apply(as.matrix(1:ninfreps), 1, SaveInfRepDataAsRDataHelper, t4 = t4, l3 = l3, curr_samp = curr_samp))
  #outp <- rbindlist(apply(as.matrix(1:1), 1, SaveInfRepDataAsRDataHelper, t4 = t4, l3 = l3, curr_samp = curr_samp))
  outp2 <- merge(outp, tx2gene, by = "tx_id")
  SaveInfRepsAsRCompTime <- proc.time() - startTime

  if(GibbsSamps==TRUE){
    if(is.null(direc_to_save_res)){
      save_dir <- paste0(getwd(), "GibbsSamps/")
    }else{
      save_dir <- paste0(direc_to_save_res, "GibbsSamps/")
    }

    attr(outp2, "save_dir") <- save_dir
    attr(outp2, "ngibbs") <- ninfreps
    attr(outp2, "ninfreps") <- ninfreps
    attr(outp2, "Sample") <- curr_samp
    nam <- paste0("GibbsSamps", curr_samp)
    assign(nam, outp2)

    if(!dir.exists(save_dir)){
      dir.create(save_dir)
    }
    save(list = c(nam, "SaveInfRepsAsRCompTime", countsFromAbundance), file = paste0(save_dir, "GibbsSamps", curr_samp, ".RData" ))
  }else{
    if(is.null(direc_to_save_res)){
      save_dir <- paste0(getwd(), "BootSamps/")
    }else{
      save_dir <- paste0(direc_to_save_res, "BootSamps/")
    }
    attr(outp2, "save_dir") <- save_dir
    attr(outp2, "nboot") <- ninfreps
    attr(outp2, "ninfreps") <- ninfreps
    attr(outp2, "Sample") <- curr_samp
    nam <- paste0("BootSamps", curr_samp)
    assign(nam, outp2)

    if(!dir.exists(save_dir)){
      dir.create(save_dir)
    }
    save(list = c(nam, "SaveInfRepsAsRCompTime", "countsFromAbundance"), file = paste0(save_dir, "BootSamps", curr_samp, ".RData" ))
  }
}

SaveInfRepDataAsRDataHelper <- function(x, t4, l3, curr_samp){
  cntcol <- as.character(x)
  colst5 <- c(cntcol, "gene_id", "NTrans", "tx_id")
  t5 <- t4[,colst5]

  colnames(t5)[colnames(t5)==cntcol] <- paste0(curr_samp, "Cnt")

  temp1 <- cntsToTPM(cnts = t5, nsamp = 1, len = l3, samps = curr_samp)
  temp1$infRepNum <- x

  t6 <- t5[,c("tx_id", paste0(curr_samp, "Cnt"))]
  temp2 <- merge(temp1, t6, by = "tx_id")
  return(temp2)

}



#Inputs
#cnts dataframe of counts with columns with names Sample1Cnt, Sample2Cnt, ... and a column for tx_id
#Additionally, if this does nt have length columns (Sample1Len, Sample2Len, etc) specify len argument
#nsamp: number of biological samples
#len: An optional argument with length information.  dataframe with (Sample1Len, Sample2Len, etc) and a tx_id column at the end
#samps: An optional vector containing the sample names.  Need this if sample names are not just paste0("Sample", 1:nsamp) without any missing.

#' Compute TPMs from counts and effective lengths
#'
#' \code{cntsToTPM} computes TPMs from counts and lengths.  Rownames need
#'
#' @param cnts is a dataframe of counts that has thecolum nnames Sample1Cnt, Sample2Cnt, ... and a column for tx_id.  Needs to also contain the (effective) length information Sample1Len, Sample2Len, ... unless the argument len is also specified
#' @param nsamp is the number of samples/biological replicates
#' @param len is an optional dataframe of length information.  len has (nsamp+1) columns with columnnames Sample1Len, Sample2Len, ... anda column for tx_id. Not needed if length information is contained in cnts.
#' @param samps is an optional vector containing the sample names.  Need to specify this if sample names are not just paste0("Sample", 1:nsamp) without any missing.
#' @return data.frame of TPM information, with column names Sample1TPM, Sample2TPM, etc.
cntsToTPM <- function(cnts, nsamp, len = NULL, samps = NULL){
  #print("head of samps within cntsToTPM is ")
  #print(head(samps))
  if(is.null(len)==TRUE){
    cnts2 <- cnts
  }else{
    #Confirm that the transcripts are in the same order in the new TPM and abGene
    if(sum(rownames(len)!=rownames(cnts)) != 0){
      stop("Transcript row names are not aligned, dataframes are not in the same order")
    }
    len2 <- data.frame(len)
    if(is.null(len2$tx_id)){
      len2$tx_id <- rownames(len2)
    }
    if(is.null(cnts$tx_id)){
      cnts$tx_id <- rownames(cnts)
    }
    cnts2 <- merge(cnts, len2, by = "tx_id")
    rownames(cnts2) <- cnts2$tx_id
    #lengthCols <- which(colnames(len) %in% paste0("Sample", 1:nsamp, "Len"))
    #t3 <- len[,lengthCols]
  }
  if(is.null(samps)==TRUE){
    countCols <- which(colnames(cnts2) %in% paste0("Sample", 1:nsamp, "Cnt"))
    lengthCols <- which(colnames(cnts2) %in% paste0("Sample", 1:nsamp, "Len"))
  }else{
    countCols <- which(colnames(cnts2) %in% paste0(samps, "Cnt"))
    lengthCols <- which(colnames(cnts2) %in% paste0(samps, "Len"))
  }

  #print("head of countCols in cntsToTPM is ")
  #print(head(countCols))

  #print("head of lengthCols in cntsToTPM is ")
  #print(head(lengthCols))

  t2 <- cnts2[,countCols, drop = FALSE]
  t3 <- cnts2[,lengthCols, drop = FALSE]

  #print("head of colnames of t2 (countscols) in cntsToTPM is ")
  #print(head(colnames(t2)))
  #print("head of colnames of t3 (lengthcols) in cntsToTPM is ")
  #print(head(colnames(t3)))

  z <- t2/t3 # Will do element wise division of counts and (effective) lengths as desired

  #print(paste0("length of samps within cntstoTPM is ", length(samps)))
  #print(paste0("nsamp within cntstoTPM is ", nsamp))
  if(is.null(samps)==TRUE){
    TPM <- as.data.frame(1e6 * apply(as.matrix(1:nsamp), 1, function(x) {z[,x]/sum(z[,x])}))
  }else{
    TPM <- as.data.frame(1e6 * apply(as.matrix(1:length(samps)), 1, function(x) {z[,x]/sum(z[,x])}))
  }

  #print("head of colnames of TPM before modification within cntsToTPM is ")
  #print(head(colnames(TPM)))

  rownames(TPM) <- rownames(t2)
  if(is.null(samps)==TRUE){
    colnames(TPM) <- paste0("Sample", 1:nsamp, "TPM")
  }else{
    colnames(TPM) <- paste0(samps, "TPM")
  }

  #print("head of colnames of TPM after modification within cntsToTPM is ")
  #print(head(colnames(TPM)))
  TPM$tx_id <- rownames(TPM)
  return(TPM)
}



#' Save the datasets that contain data across all inferential replicates and samples for a particular subset of genes
#' @inheritParams SaveInfRepDataAsRData
#' @inheritParams DRIMSeqFilter
#' @inheritParams prepareData
#' @inheritParams sumToGene
#' @param SalmonFilesDir is the directory the Salmon quantification results are saved in
#' @param save_dir is the outer directory to save the full inferential replicate datasets in.  Datasets can get quite large with a large number of samples or small number of parts so choose a directory with plenty of free space.
#' Specified directory should be the same in \code{\link{SaveFullinfRepDat}}, \code{\link{SaveWithinSubjCovMatrices}}, and \code{\link{SaveGeneLevelFiles}}.
#' @param filteredgenenames is a character vector of all genenames that pass filtering
#' @param abDatasetsFiltered is a list of dataframes that contains filtered TPM measurements for genes/transcripts that pass filtering.  Is output by \code{\link{DRIMSeqFilter}}
#' @param nparts is the total number of parts to split the \code{filteredgenenames} list into when saving the datasets.  See details.
#' @param curr_part_num is the current part number that is being run.  The genelist specified in \code{filteredgenenames} is split into \code{nparts} equally sized chunks. See the example in (3)SaveNecessaryDatasetsForCompDTUReg.R within the package's SampleCode folder
#'
#' @details This function is used to save the necessary inferential replicate datasets.  These files can get quite large with a large number of genes and/or a large number or biological samples or inferential replicates.  The parameter \code{nparts} controls how many chunks the full data is split into based
#' on splitting gene name list \code{filteredgenenames} into \code{nparts} chunks.  Increasing \code{nparts} can help mitigate issues with the resulting files becoming too large.
#' For example, for a ten sample analysis with 100 bootstrap samples we set nparts to be 10 to result in files that are around 150 MB each.  The files saved by this function will be used by \code{\link{SaveWithinSubjCovMatrices}} and \code{\link{SaveGeneLevelFiles}}.
#'  See the file (3)SaveNecessaryDatasetsForCompDTUReg.R in the package's SampleCode folder for example code.
#'
#' @return The function will save three separate .RData files that contains all bootstrap or Gibbs replicates for all samples for the genes that exist for the current part of the data (specified by \code{curr_part_num}).  Specifically, the abDataset cntDataset files containing all bootstrap/Gibbs samples
#' for the current part will be saved in the sub directories "infRepsabDatasets/" and "infRepscntDatasets/" respectively.  See the documentation for the function \code{\link{sumToGene}} for more information on the abDataset and cntDataset files.  Additionally, a data.frame that contains TPM and count information for all genes
#' in the current part for all samples is saved in the subdirectory "infRepsFullinfRepDat/".
#' @export SaveFullinfRepDat
SaveFullinfRepDat <- function(SalmonFilesDir, QuantSalmon, abDatasetsFiltered, save_dir, GibbsSamps, filteredgenenames, cntGene, key, nparts, curr_part_num){

  startTime <- proc.time()

  #Subdirectories within save_dir where the various files will be saved
  sd1 <- "infRepsabDatasets/"
  sd2 <- "infRepscntDatasets/"
  sd3 <- "infRepsFullinfRepDat/"

  if(GibbsSamps==TRUE){
    load_dir <- paste0(SalmonFilesDir, "GibbsSamps")
    infReps <- "Gibbs"
    type <- "Gibbs"
  }else{
    load_dir <- paste0(SalmonFilesDir, "BootSamps")
    infReps <- "Boot"
    type <- "Boot"
  }

  nsamp <- nrow(key)

  #Now, need to generate dataframe that has the following columns: tx_id, Sample1TPM, Sample2TPM, ..., GibbsRep#
  #To do this, need to load in the GibbsSamps dataframe for each sample, select only a current subset of genes
  #(since using all of them at once would take up too much memory) and merge them all together
  #and then run generateData as done below

  genestouse <- filteredgenenames

  #Function to split character vector into n (approximately) equally sized chunks
  #Taken from https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  #Originally by mathheadinclouds and edited by Dis Shishkov
  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  genestouse_split <- chunk2(genestouse, nparts)

  #Gene list for all genes included in cntGene, even if some of them only have 1 transcript or are filtered out, etc
  #If passing cntGeneFiltered as the argument cntGene, this will be equivalent to only using the genes to be used in the compositional analysis
  #Generally this will be the case
  all_genes_annotation_split <- chunk2(unique(sort(cntGene$gene_id)), nparts)

  #Gene list only for those genes that will be used in the compositional analysis
  curr_genes <- genestouse_split[[curr_part_num]]


  curr_genes_all_genes_annotation <- all_genes_annotation_split[[curr_part_num]]


  infRepFiles <- gtools::mixedsort(list.files(load_dir, pattern = ".RData", full.names = TRUE))

  #Load the first sample to get the ngibbs/nboot value and assign that value to ninfreps
  d <- loadRData(infRepFiles[1])
  if(GibbsSamps==TRUE){
    ninfreps <- attr(d, "ngibbs")
  }else{
    ninfreps <- attr(d, "nboot")
  }


  lengths <- data.frame(QuantSalmon$length)
  lengths$tx_id <- rownames(lengths)

  InfRepsKey <- data.frame(rep(key$Sample, ninfreps), rep(key$Condition, ninfreps), rep(key$Identifier, ninfreps))
  colnames(InfRepsKey) <- c("Sample", "Condition", "Identifier")



  #x is which of the samples from infRepFiles is being used
  f1 <- function(x, genes, infRepFiles){
    print(paste0("Currently loading file ", x))
    loadRData <- function(fileName){
      #loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }
    d <- loadRData(x)
    #Apply over each of the genestouse_split subsets
    d2 <- subset(d, d$gene_id %in% genes)
    #d2 <- lapply(genestouse_split, function(x, d){return(subset(d, d$gene_id %in% x))}, d=d)
    return(d2)
  }

  #Combine the data for the current subset of genes to all biological samples
  #Note that this subset comes from all genes from cntGene, not just the genes that will be used for the comp analysis
  #This could be useful if you ever wanted a dataset that combined inferential replicates for genes that even have
  #only one transcript
  AllDat <- lapply(infRepFiles, f1, genes = curr_genes_all_genes_annotation, infRepFiles = infRepFiles)

  #Create data file, and call it datt
  datt <- AllDat[[1]]

  #i here indexes biological samples, not the different parts the array_val indexes
  for (i in 2:length(AllDat)){
    print(paste0("Currently processing biological sample ", i, " out of ", length(AllDat)))
    curr <- AllDat[[i]]

    #Check and make sure the files are in the same order such that the trans/Gibbs rep combos match up properly
    if(sum(paste0(datt$tx_id, "GibbsRep", datt$GibbsRep)!=paste0(curr$tx_id, "GibbsRep", curr$GibbsRep))!=0){
      stop("Dataframes are not in the right order and the code will not work properly")
    }
    col <- colnames(curr)[grep("TPM", colnames(curr))]
    datt$newcol <-  curr[,col, with = F]
    #datt$newcol <- data.table:::`[.data.table`(curr, , col, with = F)
    colnames(datt)[colnames(datt)=="newcol"] <- col

    col2 <- colnames(curr)[grep("Cnt", colnames(curr))]
    datt$newcol2 <-  curr[,col2, with = F]
    #datt$newcol2 <-  data.table:::`[.data.table`(curr, , col2, with = F)
    colnames(datt)[colnames(datt)=="newcol2"] <- col2

    rm(curr)
    gc()

  }

  #datt2 <- datt[order(datt$gene_id, datt$tx_id, datt$GibbsRep),]
  if(!dir.exists(paste0(save_dir, sd3))){
    dir.create(paste0(save_dir, sd3))
  }

  assign(paste0("CompTimedattPart", curr_part_num), proc.time() - startTime)
  nam <- c("datt", paste0("CompTimedattPart", curr_part_num))
  print(paste0("Directory where files will be saved is ", paste0(save_dir, sd3, "FullinfRepDatPart")))
  save(list = nam, file = paste0(save_dir, sd3, "FullinfRepDatPart", curr_part_num, ".RData"))


  rm(datt)
  rm(AllDat)
  gc()


  #stop("Code below generated abDatasets/cntDatasets for all Gibbs replicates, only run above this line for now")

  #Combine the data for the current subset of genes to all biological samples
  #Note that this subset is only genes that are used in abDatasets, unlike the same lines of code above
  AllDat2 <- lapply(infRepFiles, f1, genes = curr_genes, infRepFiles = infRepFiles)

  #Create data file, datt
  datt2 <- AllDat2[[1]]

  #i here indexes biological samples, not the different parts the array_val indexes
  for (i in 2:length(AllDat2)){
    print(paste0("Currently processing biological sample ", i, " out of ", length(AllDat2)))
    curr <- AllDat2[[i]]

    #Check and make sure the files are in the same order such that the trans/Gibbs rep combos match up properly
    if(sum(paste0(datt2$tx_id, "GibbsRep", datt2$GibbsRep)!=paste0(curr$tx_id, "GibbsRep", curr$GibbsRep))!=0){
      stop("Dataframes are not in the right order and the code will not work properly")
    }
    col <- colnames(curr)[grep("TPM", colnames(curr))]
    datt2$newcol <-  curr[,col, with = F]
    #datt2$newcol <- data.table:::`[.data.table`(curr, , col, with = F)
    colnames(datt2)[colnames(datt2)=="newcol"] <- col

    col2 <- colnames(curr)[grep("Cnt", colnames(curr))]
    datt2$newcol2 <-  curr[,col2, with = F]
    #datt2$newcol2 <- data.table:::`[.data.table`(curr, , col2, with = F)
    colnames(datt2)[colnames(datt2)=="newcol2"] <- col2

    rm(curr)
    gc()

  }

    abGNO <- lapply(curr_genes, generateData, dat = datt2, nsamp = nsamp, abundance = TRUE, abCompDatasets = abDatasetsFiltered, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, infReps = infReps, ninfreps = ninfreps)
    names(abGNO) <- curr_genes

    assign(paste0("abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num), abGNO)
    obj2 <- c(paste0("abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num))

    rm(abGNO)
    gc()

    if(!dir.exists(paste0(save_dir, sd1))){
      dir.create(paste0(save_dir, sd1))
    }

    save(InfRepsKey, ninfreps, nsamp, list = obj2, file = paste0(save_dir, sd1, "abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num, ".RData"))


    rm(list = paste0("abDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num))
    gc()


    cntGNO <- lapply(curr_genes, generateData, dat = datt2, nsamp = nsamp, abundance = FALSE, abCompDatasets = abDatasetsFiltered, useOtherGroups = FALSE, useExistingOtherGroups = FALSE, infReps = infReps, ninfreps = ninfreps)
    names(cntGNO) <- curr_genes

    assign(paste0("cntDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num), cntGNO)
    obj2 <- c(paste0("cntDatasets", type, "NoOtherGroupsFilteredPart", curr_part_num))

    rm(cntGNO)
    gc()

    if(!dir.exists(paste0(save_dir, sd2))){
      dir.create(paste0(save_dir, sd2))
    }

    #save(abDatasetsGibbsRepsNoOtherGroups, InfRepsKey, ninfreps, nsamp, file = "abDatasetsGibbsRepsNoOtherGroups.RData")
    save(InfRepsKey, ninfreps, nsamp, list = obj2, file = paste0(save_dir, sd2, "cntDatasets", type,  "NoOtherGroupsFilteredPart", curr_part_num, ".RData"))
  


}


#Neat little function that will load an .RData object and save its contents to whatever as the results of the output
#For example can do d <- loadRData("file.RData") and whatever is loaded is saved as d
#Modified from stackexchange user ricardo at this link https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName, objNameToGet = NULL){
  #loads an RData file, and returns it
  load(fileName)
  #print(ls()[ls() != "fileName"])
  if(is.null(objNameToGet)){
    rm(objNameToGet)
    #print(ls()[ls() != "fileName"])
    return(get(ls()[ls() != "fileName"]))
  }else{
    return(get(objNameToGet))
  }
  
}


#' Save the ilr-transformed within subject covariance matricies needed for a \emph{CompDTUme} analysis
#' @inheritParams SaveInfRepDataAsRData
#' @inheritParams SaveFullinfRepDat
#' @inheritParams prepareData
#' @inheritParams SaveGeneLevelFiles
#' @inheritParams CorrectLowExpression
#' @param CLE is TRUE if the \code{\link{CorrectLowExpression}} procedure is to be used to correct low expression and FALSE if not.  Should be the same value in both \code{\link{SaveWithinSubjCovMatrices}} and \code{\link{SaveGeneLevelFiles}}.
#'
#' @details This function is used to save the necessary within subject covariance matrices. See the file (3)SaveNecessaryDatasetsForCompDTUReg.R in the package's SampleCode folder for example code.
#'
#' @return This function will save the sample-specific mean and covariance values across all inferential replicates on the \code{\link{ilr}} scale for all genes within the current part number \code{curr_part_num}.  These files will be saved within the subdirectory
#' "ilrMeansCovs/" in the \code{save_dir} folder and will be used by \code{\link{SaveGeneLevelFiles}} if present.
#' @export SaveWithinSubjCovMatrices
SaveWithinSubjCovMatrices <- function(directory, save_dir, GibbsSamps, curr_part_num, nsamp, CLE = TRUE, CLEParam = 0.05){
  #Directory where the infRepabDatasets are saved by SaveFullinfRepDat
  abDatasetsSubDir <- "infRepsabDatasets/"

  #Directory to save the within subject covariance matrices
  #ilrMeansCovsSubDir <- "ilrMeansCovs/"
  ilrMeansCovsSubDir <- "WithinSubjectCovarianceMatrices/"

  if(GibbsSamps==TRUE){
    dirpiece <- "Gibbs"
  }else{
    dirpiece <- "Boot"
  }


    load(paste0(save_dir, abDatasetsSubDir, "abDatasets", dirpiece, "Part", curr_part_num, ".RData"))
    load(paste0(save_dir, abDatasetsSubDir, "abDatasets", dirpiece, "NoOtherGroupsPart", curr_part_num, ".RData"))

    abDatasetsToUse1 <- get(paste0("abDatasets", dirpiece, "Part", curr_part_num))
    genestouse1 <- names(abDatasetsToUse1)

    abDatasetsToUse2 <- get(paste0("abDatasets", dirpiece, "NoOtherGroupsPart", curr_part_num))
    genestouse2 <- names(abDatasetsToUse2)

    load(paste0(save_dir, abDatasetsSubDir, "abDatasets", dirpiece, "NoOtherGroupsFilteredPart", curr_part_num, ".RData"))

    abDatasetsToUse3 <- get(paste0("abDatasets", dirpiece, "NoOtherGroupsFilteredPart", curr_part_num))
    genestouse3 <- names(abDatasetsToUse3)

    
    load(paste0(directory, "cntGenecntsScaledTPMFiltered.RData"))
    #genen <- names(abDatasetsToUse2)
    filterabDatasets <- function(x, cntGeneFiltered, abDatasetsToUse2){
      curr_gene <- x
      txs <- cntGeneFiltered$tx_id[cntGeneFiltered$gene_id==curr_gene]
      curr <- abDatasetsToUse2[[curr_gene]]
      curr2 <- curr[,colnames(curr) %in% txs, drop = FALSE]
      if(is.null(ncol(curr2))){
        return(NULL)
      }

      #if(ncol(curr2)!=length(txs)){stop("ree")}

      if(ncol(curr2)==0){
        curr2 <- NULL
      }
      return(curr2)
    }

    #abDatasetsToUse3 <- lapply(genen, filterabDatasets, cntGeneFiltered = cntGeneFiltered, abDatasetsToUse2 = abDatasetsToUse2)
    #abDatasetsToUse3 <- laply(genen, filterabDatasets, cntGeneFiltered = cntGeneFiltered, abDatasetsToUse2 = abDatasetsToUse2, .inform=T, .progress = "text")

    #names(abDatasetsToUse3) <- genen

    genestouse3 <- names(abDatasetsToUse3)

    d1 <- proc.time()
    ilrMeansCovsNoOtherGroupsFiltered <- lapply(genestouse3, calcIlrMeansCovs, dat = abDatasetsToUse3, nsamp = nsamp, CLE = CLE, CLEParam = CLEParam)
    proc.time() - d1
    #ilrMeansCovsNoOtherGroupsFiltered <- laply(genestouse3, calcIlrMeansCovs, dat = abDatasetsToUse3, nsamp = nsamp, .inform=T, .progress = T)
    names(ilrMeansCovsNoOtherGroupsFiltered) <- genestouse3

    assign(paste0("ilrMeansCovsNoOtherGroupsFilteredPart", curr_part_num), ilrMeansCovsNoOtherGroupsFiltered)
    #assign(paste0("WithinSubjectCovarianceMatricesFilteredPart", curr_part_num), ilrMeansCovsNoOtherGroupsFiltered)

    if(!dir.exists(paste0(save_dir, ilrMeansCovsSubDir))){
      dir.create(paste0(save_dir, ilrMeansCovsSubDir))
    }
    save(list = paste0("ilrMeansCovsNoOtherGroupsFilteredPart", curr_part_num),
         file = paste0(save_dir, ilrMeansCovsSubDir,"ilrMeansCovsNoOtherGroupsFilteredPart", curr_part_num, ".RData"))
    # save(list = paste0("WithinSubjectCovarianceMatricesFilteredPart", curr_part_num),
    #      file = paste0(save_dir, ilrMeansCovsSubDir,"WithinSubjectCovarianceMatricesFilteredPart", curr_part_num, ".RData"))

  
}





calcIlrMeansCovs <- function(x, dat, nsamp, CLE = TRUE, CLEParam = 0.05){
  test1 <- dat[[x]]
  if(is.null(test1)==TRUE){
    return(NULL)
  }
  if(ncol(test1)==1){
    return(NULL)
  }
  test2 <- test1[gtools::mixedorder(rownames(test1)),]

  SampNames <- plyr::laply(rownames(test2), function(x){strsplit(x, "TPM")[[1]][1]})
  UniqueSampNames <- gtools::mixedsort(unique(SampNames))

  UniqueSampNumbers <- as.numeric(plyr::laply(UniqueSampNames, function(x){strsplit(x, "Sample")[[1]][2]}))
  #codetest <- apply(test2, 1, ilr)
  #Calculate the ilr values and their means and covs separately for each biological sample
  Means <- list()
  Covs <- list()


    vals <- UniqueSampNumbers

  for (i in 1:length(vals)){
    #minrange <- (i-1)*100 + 1
    #maxrange <- i * 100

    curr_samp <- UniqueSampNames[i]
    #curr_samp <- SampNames[i]

    #test3 <- test2[minrange:maxrange,]
    #test3 <- test2[UniqueSampNames==curr_samp,]

    #Dimension of test3 is nsamp*ncol of the current abDataset file

    #Only select observations corresponding to the current sample
    if(length(SampNames)!=nrow(test2)){
      stop("Subsetting of test2 within calcilrMeansCovs has the incorrect number of dimensions, check it and try again")
    }
    test3 <- test2[SampNames==curr_samp,]

    test3SampNames <- plyr::laply(rownames(test3), function(x){strsplit(x, "TPM")[[1]][1]})
    if(length(unique(test3SampNames))!=1){
      stop("Something is wrong in calcIlrMeansCovs, observations from more than 1 sample are being used at the same time")
    }
    #If no rows, there are no Gibbs samples for this sample and we can move to the next sample
    if(nrow(test3)==0){
      Means[[i]] <- "No Gibbs/Boot Samples Available"
      Covs[[i]] <- "No Gibbs/Boot Samples Available"
      next
    }

    #Note that we can just apply the ilr transform on the TPM (and not on the relative proportion abundances) because the ilr values will be the same
    # in either case since the ilr of a normalized vector (or closed, which is normalilzed summing to 1) is the same as the ilr of the same vector
    # non-normalized.  Can confirm that test4T below is the same as test4 (except for small computation related error)

    # test3T <- ccomp(test3, total=1)
    # test4T <- apply(test3T, 1, ilr)

    #Dim of test 4 will be (ncol of current abDataset file - 1) * nsamp
    #This is not the final form, will be transposed by test5

    if(CLE==TRUE){
      test3_2 <- CorrectLowExpression(test3)
      test4 <- apply(test3_2, 1, function(x){compositions::ilr(x)})
    }else{
      test4 <- apply(test3, 1, function(x){compositions::ilr(x)})
    }


    #Dim of test 5 will be nsamp * (ncol of current abDataset file - 1), as intended
    if(is.null(ncol(test4))){
      test5 <- data.frame(as.matrix(test4, ncol = 1))
    }else{
      test5 <- data.frame(t(test4))
    }
    #Columns don't have a practical interpretation here, so change the names to null
    #See Analyzing Compositional Data with R book, p 45 for an explanation of this
    colnames(test5) <- NULL

    ilrMeans <- as.data.frame(t(as.matrix(colMeans(test5))))
    colnames(ilrMeans) <- NULL
    rownames(ilrMeans) <- NULL

    ilrCov <- stats::cov(test5)

    Means[[i]] <- ilrMeans
    #names(Means)[i] <- paste0("Sample", i)
    Covs[[i]] <- ilrCov
    #names(Covs)[i] <- paste0("Sample", i)
  }
  #names(Means) <- paste0("Sample", 1:nsamp)
  names(Means) <- UniqueSampNames
  #names(Covs) <- paste0("Sample", 1:nsamp)
  names(Covs) <- UniqueSampNames
  #Means2 <- rbindlist(Means)
  #rownames(Means2) <- paste0("Sample", 1:nsamp)
  return(list(ilrMeans = Means, ilrCovs = Covs))
}


#' Save a file for each gene containing all information necessary to run CompDTUReg
#' @inheritParams SaveInfRepDataAsRData
#' @inheritParams SaveWithinSubjCovMatrices
#' @inheritParams CorrectLowExpression
#' @param directory is the directory where previous datasets are saved by \code{\link{sumToGene}} and \code{\link{DRIMSeqFilter}}
#' @param GeneLevelFilesSaveDir is the directory the gene level files will be saved to
#' @param curr_part_num if the current part number to save results for.  See \code{\link{SaveFullinfRepDat}} for more details.
#' @param useInferentialReplicates is set to TRUE if inferential replicates are to be used in the analysis and FALSE otherwise.   If useInferentialReplicates if set to FALSE the gene-level files will only contain data corresponding to the point estimates (output as Y)
#' and if it is TRUE the gene-level files will also contain data corresponding to the inferential replicates (output as YInfRep)  If FALSE the argument \code{GibbsSamps} is ignored.
#'
#' @return The function saves gene-level .RData files within the subdirectory \code{GeneLevelFilesSaveDir} of \code{save_dir} that contain all necessary
#' information to run a \code{\link{CompDTUReg}} analysis, including inferential replicates (if any).  The name of each file is the name of the current gene.  Specifically, these files contain: \cr
#'  (1) key, which contains the sample names (in the Sample column), condition information (in the Condition column), and an identifier in the form of Sample1, Sample2, Sample3, etc \cr
#'  (2) genename, the name of the current gene (that will be the name of the outputted file as well) \cr
#'  (3) Group, which is a factor of the condition information specified in the Condition column of key \cr
#'  (4) samps, the identifiers of the sample (Sample1, Sample2, Sample3, etc) contained within key$Sample \cr
#'  (5) Y, the ilr-transformed matrix of response values (ie the ilr coordinates) corresponding to the point estimates for all samples for the current gene for use with CompDTU. \cr
#'  (6) YInfRep, the ilr-transformed matrix  of response values (ie the ilr coordinates) corresponding to all inferential replicates for all samples for the current gene for use with CompDTUme.\cr
#'  Will be in the order Sample1InfRep1, Sample2InfRep1, Sample3InfRep1, etc and will only be output if useInferentialReplicates is TRUE. \cr
#'  (7) mean.withinhat The overall within-sample covariance matrix (that is the mean of each sample's witin-subject covariance matrix) used to correct the total covariance in the CompDTUme method.  Only output if useInferentialReplicates is TRUE. \cr
#'
#' @details See the file (3)SaveNecessaryDatasetsForCompDTUReg.R in the package's SampleCode folder for example code.
#' @export SaveGeneLevelFiles
SaveGeneLevelFiles <- function(directory, GeneLevelFilesSaveDir, curr_part_num, useInferentialReplicates = TRUE, GibbsSamps,
                               CLE = TRUE, CLEParam = 0.05){
    setwd(directory)


    load("abGeneFiltered.RData")
    load("cntGenecntsScaledTPMFiltered.RData")
    load("abDatasetsNoOtherGroupsFiltered.RData")
    load("tx2gene.RData")

    abGeneFiltered <- abGeneFiltered
    cntGeneFiltered <- cntGeneFiltered
    abDatasetsFiltered <- abDatasetsFiltered

    abGene <- abGeneFiltered
    cntGene <- cntGeneFiltered
    abDatasets <- abDatasetsFiltered

  key <- key

  Group <- key$Condition
  sampstouse <- key$Identifier

  genestouse <- names(abDatasets)

  if(GibbsSamps==TRUE){
    dirpiece <- "Gibbs"
  }else{
    dirpiece <- "Boot"
  }

  if(useInferentialReplicates==TRUE){
    abDatasetsSubDir <- "infRepsabDatasets/"
    #ilrMeansCovsSubDir <- "ilrMeansCovs/"
    ilrMeansCovsSubDir <- "WithinSubjectCovarianceMatrices/"

      ilrMeansCovs <- loadRData(paste0(directory, ilrMeansCovsSubDir, "ilrMeansCovsNoOtherGroupsFilteredPart", curr_part_num, ".RData"))
      newAbDatasetsinfRepsFinal <- loadRData(paste0(directory, abDatasetsSubDir, "abDatasets", dirpiece, "NoOtherGroupsFilteredPart", curr_part_num,".RData"))

    newAbDatasetsinfRepsFinalSub <- newAbDatasetsinfRepsFinal[names(newAbDatasetsinfRepsFinal) %in% genestouse]
    ilrMeansCovsSub <- ilrMeansCovs[names(ilrMeansCovs) %in% genestouse]

    numg <- max(length(ilrMeansCovsSub), length(newAbDatasetsinfRepsFinalSub))
    gnames <- names(ilrMeansCovsSub)
  }else{
    gnames <- genestouse
    numg <- length(genestouse)
  }


  if(!dir.exists(GeneLevelFilesSaveDir)){
    dir.create(GeneLevelFilesSaveDir)
  }



  for(j in 1:numg){
    curr_gene <- gnames[j]

    if(j %% 100==0){
      print(paste0("Currently saving files for gene number ", j, " out of ", numg))
    }

    if(file.exists(file = paste0(GeneLevelFilesSaveDir, curr_gene, ".RData"))){next}

    if(!(curr_gene %in% genestouse)){
      next
    }

    abDatasetsToUse <- abDatasets[curr_gene]
    if(useInferentialReplicates==TRUE){
      ilrMeansCovs <- ilrMeansCovsSub[curr_gene]
      newAbDatasetsInfRepsFinal <- newAbDatasetsinfRepsFinalSub[curr_gene]
    }


    if(is.null(abDatasetsToUse[[1]])){
      next
    }

    if(ncol(abDatasetsToUse[[1]])==1){
      next
    }

    if(CLE==TRUE){
      Y <- unclass(compositions::ilr(compositions::ccomp(CorrectLowExpression(abDatasetsToUse[[1]], CLEParam), total = 1)))
    }else{
      Y <- unclass(compositions::ilr(compositions::ccomp(abDatasetsToUse[[1]], total = 1)))
    }

    #Remove "TPM" from rownames of Y if it is still present, as these no longer correspond to TPMs
    rownames(Y) <- strsplit(rownames(Y), "TPM")

    #Files by default contain both means and covariances of inferential reps (on ilr scale)
    #For this purpose, return only the covariances
    if(useInferentialReplicates==TRUE){
      if(is.null(ilrMeansCovs)==FALSE){
        ReturnilrCovsOnly <- function(x, ilrMeansCovs){
          temp1 <- ilrMeansCovs[[x]]$ilrCovs
          return(temp1)
        }

        ilrCovsToUse <- ReturnilrCovsOnly(x = curr_gene, ilrMeansCovs = ilrMeansCovs)
        # if(!is.null(ilrCovsToUse)){
        #   names(ilrCovsToUse) <- curr_gene
        # }

      }else{
        ilrCovsToUse <- NULL
      }

      ilrCovsOnly <- TRUE

      GibbsCovsToUse <- ilrCovsToUse

      if(is.null(newAbDatasetsInfRepsFinal)){
        YInfRep <- NULL
      }else if(ncol(newAbDatasetsInfRepsFinal[[1]])==1){
        YInfRep <- NULL
      }else{

        if(CLE == TRUE){
          YInfRep <- unclass(compositions::ilr(compositions::ccomp(CorrectLowExpression(newAbDatasetsInfRepsFinal[[1]], CLEParam), total = 1)))
        }else{
          YInfRep <- unclass(compositions::ilr(compositions::ccomp(newAbDatasetsInfRepsFinal[[1]], total = 1)))
        }

        rownames(YInfRep) <- lapply(strsplit(rownames(YInfRep), "TPM"), FUN = function(x){paste0(x[1], x[2])})
      }

      if(is.null(GibbsCovsToUse)){
        mean.withinhat <- NULL
      }else{
        nsamp2 <- length(GibbsCovsToUse)
        mean.withinhat <- tryCatch(Reduce("+", GibbsCovsToUse)/nsamp2, error = function(x){})
      }
    }




    Group <- key$Condition
    samps <- key$Identifier
    genename <- curr_gene

    if(useInferentialReplicates==TRUE){
      save(key, Group, Y, samps, mean.withinhat, YInfRep, genename, file = paste0(GeneLevelFilesSaveDir, curr_gene, ".RData"))
      rm(abDatasetsToUse)
      rm(GibbsCovsToUse)
      rm(ilrMeansCovs)
      rm(newAbDatasetsInfRepsFinal)
      rm(Y)
      rm(mean.withinhat)
      rm(YInfRep)
      rm(genename)
      gc()
    }else{
      save(key, Group, Y, samps, genename, file = paste0(GeneLevelFilesSaveDir, curr_gene, ".RData"))
      rm(abDatasetsToUse)
      rm(Y)
      rm(genename)
      gc()
    }

  }
}


#' Correct sample/gene combinations that have expression values of 0 or close to 0 to stabilize results
#' @param y is the data for the current gene/sample combination
#' @param CLEParam is the parameter (betwen 0 and 1) that controls the correction threshold (see details of \code{\link{CorrectLowExpression}} for more information)
#' @details The CLEParam parameter a works as follows: any TPM value that is less than `100*a' percent of the total gene-level
#' expression for the sample is replaced by `100*a' percent of this expression.  Mathematically, let \eqn{T_{ij}} be the TPM value for
#' transcript $j=1,..., D$ for sample $i = 1,..., n$ within a given gene with $D$ transcripts.  Any \deqn{T_{ij} < a * (T_{i1}+...+T_{iD})} will be replaced by \deqn{a * (T_{i1}+...+T_{iD})}.
#' This procedure results in relative transcript abundance fractions (RTAFs)
#' being zero only when every \eqn{T_{ij}} is equal to zero.  The value \eqn{a} could be increased or decreased to result in more or less modification to
#' the observed TPM values.  As $a$ increases, the proportions are driven closer to each other, with each converging towards \eqn{(1/D)} as \eqn{a} converges to 1
#' (and each equal to \eqn{(1/D)} for \eqn{a > 1}).  We find \deqn{a=0.05} is a good compromise that is large enough to stabilize the ilr coordinates sufficiently while
#' additionally not over-modifying the observed data, and use this value by default for all \emph{CompDTU} and \emph{CompDTUme} results. We recommend keeping this value at the default of 0.05.
#' @export CorrectLowExpression
CorrectLowExpression <- function(y, CLEParam = 0.05){
  if(is.null(y)){
    return(NULL)
  }

  if(ncol(y)==1){
    y[y==0] <- CLEParam
    return(y)
  }else{
    curr_dat <- y

    lowExp_corrected_dat <- t(apply(curr_dat, 1, CorrectLowExpressionHelper, CLEParam = CLEParam))
    return(lowExp_corrected_dat)
  }


}

CorrectLowExpressionHelper <- function(x, CLEParam = 0.05){
  #print(x)
  # if(nrow(x)!=1){
  #   stop("Number of rows is not 1")
  # }
  curr_rowSum <- sum(x)
  if(is.na(curr_rowSum)){
    return(x)
  }else if(curr_rowSum==0){
    return(x)
  }else{
    curr_dat2 <- x
    curr_dat2[curr_dat2 < CLEParam*curr_rowSum] <- CLEParam*curr_rowSum
    return(curr_dat2)
  }
}


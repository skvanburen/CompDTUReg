
# Starts the compositional DTU regression with or with bootstrap samples and return the results
#'
#'
#' \code{startCompDTUReg} start the compositional DTU regression model corresponding to a specific gene-level file.
#'
#' \code{startCompDTUReg} runs the compositional DTU regression models.  This is a helper function to run \code{\link{CompDTUReg}} from gene level files originating
#' from the file (4)RunCompositionalRegressions.R in the sample code for the package.
#' @inheritParams CompDTUReg
#' @param x The file path to a specific gene-level file
#'
#' @details This function loads results from \code{\link{SaveGeneLevelFiles}} thatcontain all input arguments aside from \code{extraPredictors} (if any),
#' and input arguments to \code{\link{CompDTUReg}} are automatically set by \code{\link{startCompDTUReg}}.
#' See the file (4)RunCompositionalRegressions.R in the package's SampleCode folder for example code.
#'
#' @return a data.frame containing the gene_id being used, p-value from the CompDTU or CompDTUme significance test for condition, and various information on
#' the current dataset being used.
#'
#' @export startCompDTUReg
startCompDTUReg <- function(x, runWithME, extraPredictors = NULL, customHypTest = FALSE, NullDesign = NULL, AltDesign = NULL){

  if(customHypTest==TRUE & (is.null(NullDesign) | is.null(AltDesign))){
    stop("To conduct a custom hypothesis test, both NullDesign and AltDesign must be specified")
  }
  
  if(customHypTest==FALSE & (!is.null(NullDesign) | !is.null(AltDesign))){
    stop("You have input a custom design matrix but customHypTest is FALSE.  Set customHypTest to TRUE and ensure both NullDesign and AltDesign are specified to run a custom hypothesis test.")
  }
  
  
  if(!is.null(NullDesign) & !is.null(AltDesign) ){
    if(ncol(NullDesign)>=ncol(AltDesign)){
      stop("The Alternative design matrix does not have more predictors than the null design matrix.  Check the specification of the two design matrices.")
    }
  }
  
  load(x)
  genename <- genename
  Group <- Group
  if(runWithME==TRUE){
    if(is.null(YInfRep)){
      stop("You have specified runWithME to be TRUE to run the CompDTUme model but the inferential replicate data appears to not exist within the gene-level file.  Did you generate the gene-level files using SaveGeneLevelFiles with useInferentialReplicates set to FALSE by mistake?")
    }
    mean.withinhat <- mean.withinhat
    
    if(customHypTest==TRUE){
      if(!all.equal(rownames(NullDesign), rownames(YInfRep)) | !all.equal(rownames(AltDesign), rownames(YInfRep))){
        stop("The rownames of the inputted design matrices do not match the rownames of the response variable, YInfRep.  Verify that the rows of the custom design matrices are in the correct order.")
      }
    }
    res <- CompDTUReg(genename = genename, Y = NULL, Group = Group, runWithME = TRUE, YInfRep = YInfRep, mean.withinhat = mean.withinhat,
                      extraPredictors = extraPredictors, customHypTest = customHypTest, NullDesign = NullDesign, AltDesign = AltDesign)
  }else{
    
    if(customHypTest==TRUE){
      if(!all.equal(rownames(NullDesign), rownames(Y)) | !all.equal(rownames(AltDesign), rownames(Y))){
        stop("The rownames of the inputted design matrices do not match the rownames of the response variable, Y.  Verify that the rows of the custom design matrices are in the correct order.")
      }
    }  
      
    res <- CompDTUReg(genename = genename, Y = Y, Group = Group, runWithME = FALSE, YInfRep = NULL, extraPredictors = extraPredictors, 
                      customHypTest = customHypTest, NullDesign = NullDesign, AltDesign = AltDesign)
  }
  return(res)
}

#' Implement the CompDTU and CompDTUme regression model
#'
#'
#' \code{CompDTUReg} runs the CompDTU and CompDTUme regression models for a specific gene.
#'
#'
#' @param genename The name of the current gene to run the method on.
#' @param Y corresponds is the ilr transformed matrix of response values for the non-inferential replicate data.  Matrix
#' should have a number of rows corresponding to the number of samples and a number of columns corresponding to the number of ilr coordinates
#' (which one less than the number transcripts remaining after filtering).  Rows need to be in the same order as either Group or the design matrices specified by AltDesign and NullDesign.  Set to NULL if you are running the model with measurement error, as it is not used.
#' @param Group A vector of condition assignments corresponding to the samples.  Needs to be ordered the same as Y or YInfRep.  Set to NULL if custom design matrices will be specified by the NullDesign and AltDesign arguments.
#' @param runWithME is TRUE/FALSE indicating whether the model should be run with measurement error or not (corresponding to CompDTU
#' and CompDTUme models respectively).  If runWithME is TRUE ensure YInfRep is non-NULL and if runWithME is FALSE ensure Y is non-NULL.
#' @param YInfRep corresponds is the ilr transformed matrix of response values for the inferential replicate data.  Matrix
#' should have a number of rows corresponding to the number of samples times the number of replicates and a number of columns corresponding to the number of ilr coordinates
#' (which is one less than the number transcripts remaining after filtering).  Rows must be ordered as Sample1InfRep1, Sample2InfRep1, Sample3Infrep1, etc. NOT Sample1InfRep1, Sample1InfRep2, Sample1InfRepe,  etc.  Set to NULL if you are running the model without measurement error, as it is not used.
#' @param mean.withinhat is the mean of the sample-specific covariance matrices of the inferential replicates (calculated on the \code{\link{ilr}} scale).
#' @param extraPredictors is an optional matrix of additional predictor values.  This should have one row per sample and one column per predictor.  The column names of the matrix will be taken as the names of the predictor.
#' The condition variable should not be included in this matrix because it is included automatically via the \code{Group} parameter.  Number of rows needs to be the number of samples, even if CompDTUme is to be run (the matrix will be replicated automatically as needed). Set to NULL if custom design matrices will be specified by the NullDesign and AltDesign arguments.
#' @param customHypTest should be set to TRUE if custom design matricies will be specified using the NullDesign and AltDesign arguments.  If TRUE the Group and extraPredictors arguments must be set to NULL.  Default is FALSE.
#' @param AltDesign specifies the design matrix corresponding to the alternative hypothesis if a custom hypothesis test is specified via customHypTest being set to TRUE.  Number of rows needs to be the number of samples, even if CompDTUme is to be run (the matrix will be replicated automatically as needed).
#' @param NullDesign specifies the design matrix corresponding to the null hypothesis if a custom hypothesis test is specified via customHypTest being set to TRUE.  Number of rows needs to be the number of samples, even if CompDTUme is to be run (the matrix will be replicated automatically as needed).
#'
#' @details  This function is run separately for each gene and the easiest way to run it will be to follow the pipeline given in the SampleCode folder of the package.
#' In particular if gene-level files created from \code{\link{SaveGeneLevelFiles}} are used the easiest way to use this is to use the helper function \code{\link{startCompDTUReg}}, which loads the gene-level results.
#' These gene-level files contain all input arguments aside from \code{extraPredictors} (if any) or custom specified null and alternative hypothesis design matrices, and input arguments are automatically set by \code{\link{startCompDTUReg}}.
#' See the file (4)RunCompositionalRegressions.R in the package's SampleCode folder for example code.
#'
#' @return a data.frame containing the gene_id being used, p-value from the CompDTU or CompDTUme significance test for condition, and various information on
#' the current dataset being used.
#'
#' @export CompDTUReg
CompDTUReg <- function(genename, Y = NULL, Group = NULL, runWithME = TRUE, YInfRep = NULL, mean.withinhat = NULL,
                       extraPredictors = NULL, customHypTest = FALSE, NullDesign = NULL, AltDesign = NULL){
  
  #If a custom hypothesis test is to be specified, set Group and extraPredictors to be NULL
  if(customHypTest==TRUE){
    Group <- NULL
    extraPredictors <- NULL
  }
  
  if(is.null(Group) & ((is.null(NullDesign) | is.null(AltDesign)))){
      stop("If Group is not specified design matrices must be specified under the null and alternative hypothesis using the NullDesign and AltDesign arguments")
  }
  
  if(customHypTest==TRUE & (is.null(NullDesign) | is.null(AltDesign))){
    stop("To conduct a custom hypothesis test, both NullDesign and AltDesign must be specified")
  }
  
  if(customHypTest==FALSE & (!is.null(NullDesign) | !is.null(AltDesign))){
    stop("You have input a custom design matrix but customHypTest is FALSE.  Set customHypTest to TRUE and ensure both NullDesign and AltDesign are specified to run a custom hypothesis test.")
  }
  
  if(!is.null(NullDesign) & !is.null(AltDesign) ){
    if(ncol(NullDesign)>=ncol(AltDesign)){
      stop("The Alternative design matrix does not have more predictors than the null.  Check the specification of the two design matrices.")
    }
  }
  
  if(runWithME==TRUE & is.null(YInfRep)){
    stop("You are attempting to run the CompDTUme model but YInfRep is NULL")
  }
  
  if(runWithME==FALSE & is.null(Y)){
    stop("You are attempting to run the CompDTU model but Y is NULL")
  }
  
  if(runWithME==TRUE & !is.null(Y)){
    stop("Set Y to be NULL if you are attempting to run the CompDTUme model")
  }
  
  if(runWithME==FALSE & !is.null(YInfRep)){
    stop("Set YInfRep to be NULL if you are attempting to run the CompDTU model")
  }
  

  if(customHypTest==FALSE){
    
    if(length(unique(Group))!=length(levels(Group))){
      warning("Check the number of levels in the Group factor, there appear to be more possible levels than are actually used by the data and this will result in incorrect inference.")
    }
    
    Group2 <- Group
    nsamp <- length(Group2)
    ncond <- length(unique(Group2))
    
  }else{
    nsamp <- nrow(AltDesign)
  }
  




  if(runWithME==TRUE){
    
      if(customHypTest==FALSE){
        ninfreps <- nrow(YInfRep)/length(Group)
        Group2InfReps <- rep(Group2, ninfreps)
      }else{
        ninfreps <- nrow(YInfRep)/nrow(AltDesign)
      }


      if(!is.null(AltDesign)){
        #XAlt <- AltDesign
        XAlt <- do.call(rbind, replicate(ninfreps, AltDesign, simplify=FALSE))
        
        if(nrow(XAlt)!=nrow(YInfRep)){
          stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you input a proper design matrix for AltDesign?")
        }
        
      }else if(!is.null(extraPredictors)){
      #Repeat the extraPredictors the necessary number of times
        #Stack the new predictor matrix the necessary number of times
        ExtraPredMatrixInfReps <- do.call(rbind, replicate(ninfreps, extraPredictors, simplify=FALSE))
        XAlt <- stats::model.matrix(~Group2InfReps + ExtraPredMatrixInfReps)

        #Reassign column names to match colunm names of extraPredictors
        colnames(XAlt)[ncol(XAlt):(ncol(XAlt) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
        
        if(nrow(XAlt)!=nrow(YInfRep)){
          stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you a proper matrix for extraPredictors and that the input Group variable has the correct number of samples?")
        }
      }else{
        XAlt <- stats::model.matrix(~Group2InfReps)
        
        if(nrow(XAlt)!=nrow(YInfRep)){
          stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure the input Group variable has the correct number of samples?")
        }
      }
      



      ns <- nrow(YInfRep)
      XAltT <- t(XAlt)
      bhatalt <- solve(crossprod(XAlt)) %*% (XAltT %*% YInfRep)
      pie1 <- YInfRep - (XAlt%*%bhatalt)

      SigmaTildeAltNewModeling <- (crossprod(pie1))/ns

      #resAlt <- fastLmSVB(Y = YInfRep, X = XAlt)
      #SigmaTildeAltNewModeling <- resAlt$EstCov

      if(!is.null(NullDesign)){
        XNull <- do.call(rbind, replicate(ninfreps, NullDesign, simplify=FALSE))
        
        if(nrow(XNull)!=nrow(YInfRep)){
          stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you input a proper design matrix for NullDesign?")
        }
      }else if(!is.null(extraPredictors)){
        ExtraPredMatrixInfReps <- do.call(rbind, replicate(ninfreps, extraPredictors, simplify=FALSE))
        XNull <- stats::model.matrix(~ExtraPredMatrixInfReps)

        #Reassign column names to match colunm names of extraPredictors
        colnames(XNull)[ncol(XNull):(ncol(XNull) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
        
        if(nrow(XNull)!=nrow(YInfRep)){
          stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you a proper matrix for extraPredictors and that the input Group variable has the correct number of samples?")
        }
      }else{
        XNull <- stats::model.matrix(~1, data = Group2InfReps)
        
        if(nrow(XNull)!=nrow(YInfRep)){
          stop("The number of rows of the design matrix under the null should match the number of samples and does not.  Are you sure the input Group variable has the correct number of samples?")
        }
      }

      XNullT <- t(XNull)
      bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% YInfRep)
      pie2 <- YInfRep - (XNull%*%bhatnull)

      SigmaTildeNullNewModeling <- (crossprod(pie2))/ns

      gene_id <- genename
      ret <- data.frame(gene_id, pval_CompDTUme = NA, FStat = NA, NumDF = NA, DenomDF = NA, stringsAsFactors = F)

      if(!is.null(mean.withinhat)){
        UpdatedCovAlt <- SigmaTildeAltNewModeling - mean.withinhat
        UpdatedCovNull <- SigmaTildeNullNewModeling - mean.withinhat

        #df_residual needs to be based on nsamp not nreps or nreps*nsamp
        
        #If there is a negative variance term, the pvalue for CompDTUme will be undefined
        statement1 <- sum(diag(UpdatedCovAlt)<=0) !=0
        statement2 <- sum(diag(UpdatedCovNull)<=0) !=0
        
        if((statement1==TRUE | statement2==TRUE)){
          ret$pval_CompDTUme <- NA
        }else{
          if(customHypTest==TRUE){
            qq <- ncol(XAlt) - ncol(XNull)
          }else{
            qq <- ncond - 1
          }
          
          CompDTUmeRes <- tryCatch(calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
                                                                  lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
                                                   error = function(x){})

          if(is.null(CompDTUmeRes)==FALSE){
            ret$pval_CompDTUme <- CompDTUmeRes$pval_pillai
            ret$FStat <- CompDTUmeRes$FStat
            ret$NumDF <- CompDTUmeRes$NumDF
            ret$DenomDF <- CompDTUmeRes$DenomDF
          }


        }

      }
      return(ret)
    

  }else if(runWithME==FALSE){
    if(!is.null(AltDesign)){
      XAlt <- AltDesign
      
      if(nrow(XAlt)!=nrow(Y)){
        stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you input a proper design matrix for AltDesign?")
      }
    }else if(!is.null(extraPredictors)){
      XAlt <- stats::model.matrix(~Group2 + extraPredictors)

      #Reassign column names to match colunm names of extraPredictors
      colnames(XAlt)[ncol(XAlt):(ncol(XAlt) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
      
      if(nrow(XAlt)!=nrow(Y)){
        stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you a proper matrix for extraPredictors and that the input Group variable has the correct number of samples?")
      }
    }else{
      XAlt <- stats::model.matrix(~Group2)
      
      if(nrow(XAlt)!=nrow(Y)){
        stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure the input Group variable has the correct number of samples?")
      }
    }
    ns <- nrow(Y)
    XAltT <- t(XAlt)
    bhatalt <- solve(crossprod(XAlt)) %*% (XAltT  %*% Y)
    pie1 <- Y - XAlt%*%bhatalt
    SigmaTildeAlt <- (crossprod(pie1))/ns

    if(!is.null(NullDesign)){
      XNull <- NullDesign
      
      if(nrow(XNull)!=nrow(Y)){
        stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you input a proper design matrix for NullDesign?")
      }
    }else if(!is.null(extraPredictors)){
      XNull <- stats::model.matrix(~1 + extraPredictors, data = Group2)

      #Reassign column names to match column names of extraPredictors
      colnames(XNull)[ncol(XNull):(ncol(XNull) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
      
      if(nrow(XNull)!=nrow(Y)){
        stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure you a proper matrix for extraPredictors and that the input Group variable has the correct number of samples?")
      }

    }else{
      XNull <- stats::model.matrix(~1, data = Group2)
      
      if(nrow(XAlt)!=nrow(Y)){
        stop("The number of rows of the design matrix should match the number of samples and does not.  Are you sure the input Group variable has the correct number of samples?")
      }
    }

    XNullT <- t(XNull)
    bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% Y)
    pie2 <- Y - XNull%*%bhatnull
    #SigmaTildeNull <- (t(pie1)%*%pie1)/ns
    SigmaTildeNull <- (crossprod(pie2))/ns



    gene_id <- genename
    ret <- data.frame(gene_id, pval_CompDTU = NA,  FStat = NA, NumDF = NA, DenomDF = NA, stringsAsFactors = F)

    if(customHypTest==TRUE){
      qq <- ncol(XAlt) - ncol(XNull)
    }else{
      qq <- ncond - 1
    }
    
    CompDTURes <- tryCatch(calcPillaiPval(SigmaTildeNull = SigmaTildeNull, SigmaTildeAlt = SigmaTildeAlt,
                                           lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
                            error = function(x){})

    if(is.null(CompDTURes)==FALSE){
      ret$pval_CompDTU <- CompDTURes$pval_pillai
      ret$FStat <- CompDTURes$FStat
      ret$NumDF <- CompDTURes$NumDF
      ret$DenomDF <- CompDTURes$DenomDF
    }


  }
  return(ret)

  #End compReg
}





#' Calculate the Pillai pvalue for a categorical predictor (usually condition)
#'
#'
#' \code{calcPillaiPval} calculate the Pillai pvalue for a categorical predictor
#' @param SigmaTildeNull is the null covariance to be used
#' @param SigmaTildeAlt is the alternative covariance to be used
#' @param lm_model_fit is a lm object from the regression model if interest.  If specified it can be used to
#' extract the degrees of freedom of the residual (df_residual) otherwise this must be specified
#' @param q is degrees of freedom of the categorical predictor being tested in the current hypothesis test (usually condition).
#' For example, with 2 condition levels q is 1, with 3 levels q is 2, etc
#' @param nsamp is the number of samples used in the analysis. Note that this is the number of unique biological samples and is thus not
#' a function of the number of inferential replicates used in an analysis.
#' @param df_residual is the degrees of freedom of the residual. Can be extracted from the lm_model_fit object otherwise has to
#' be specified.  Equal to the number of samples used minus the number of total coefficients fit by the model.
#'
#' @return The p-value from the Pillai significance test.
calcPillaiPval <- function(SigmaTildeNull, SigmaTildeAlt, lm_model_fit = NULL, q, nsamp, df_residual = NA){
  if(is.null(lm_model_fit) & is.na(df_residual)){
    stop("df. residual must be specified to calcPillaiPval or else an lm object must specified in lm_model_fit to extract df.residual from")
  }

  Etilde <- nsamp * SigmaTildeAlt
  Htilde <- nsamp * (SigmaTildeNull - SigmaTildeAlt)

  vals <- diag((Htilde %*%solve(Etilde + Htilde)))
  pill_stat <- sum(vals)


  #See the Multivariate ANOVA Testing pdf document (from the SAS help file) for the necessary formulas
  #This proved to be the easiest way to calculate the statistic and are confirmed to match R's anova.mlm function
  #v <- nsamp - ncol(stats::model.matrix(lm_model_fit))

  if(!is.na(df_residual)){
    v <- df_residual
  }else{
    v <- lm_model_fit$df.residual
  }
  #v is the error/residual df- also extract from the r anova fit


  #p is the number of eigenvales (ie rank of Etilde + Htilde)
  p <- length(vals)

  s <- min(p,q)

  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (v - p - 1)

  #Formulas come from the SAS help file Multivariate ANOVA Testing and are confirmed to match R's anova.mlm function
  piece1 <- 2*n + s + 1
  piece2 <- 2*m + s + 1
  fstat_pillai <- (piece1/piece2) * (pill_stat/(s - pill_stat))

  if(fstat_pillai < 0){
    return(list(FStat = NA, NumDF = NA, DenomDF = NA, pval_pillai = NA))
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-stats::pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  
  
  return(list(FStat = fstat_pillai, NumDF = numdf_pillai, DenomDF = denomdf_pillai, pval_pillai = pval_pillai))
}

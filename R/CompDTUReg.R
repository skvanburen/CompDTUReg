
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
startCompDTUReg <- function(x, runWithME, extraPredictors = NULL){

  load(x)
  genename <- genename
  Y <- Y
  Group <- Group
  YInfRep <- YInfRep
  mean.withinhat <- mean.withinhat
  if(runWithME==TRUE){
    res <- CompDTUReg(genename = genename, Y = Y, Group = Group, runWithME = TRUE, YInfRep = YInfRep, mean.withinhat = mean.withinhat,
                      extraPredictors = extraPredictors)
  }else{
    res <- CompDTUReg(genename = genename, Y = Y, Group = Group, runWithME = FALSE, extraPredictors = extraPredictors)
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
#' (which one less than the number transcripts remaining after fltering).  Set to NULL if you are running the model with measurement error, as it is not used.
#' @param Group A vector of condition assignments corresponding to the samples
#' @param runWithME is TRUE/FALSE indicating whether the model should be run with measurement error or not (corresponding to CompDTU
#' and CompDTUme models respectively).  If runWithME is TRUE ensure YInfRep is non-NULL and if runWithME is FALSE ensure Y is non-NULL.
#' @param YInfRep corresponds is the ilr transformed matrix of response values for the inferential replicate data.  Matrix
#' should have a number of rows corresponding to the number of samples times the number of replicates and a number of columns corresponding to the number of ilr coordinates
#' (which one less than the number transcripts remaining after fltering).  Set to NULL if you are running the model without measurement error, as it is not used.
#' @param mean.withinhat is the mean of the sample-specific covariance matrices of the inferential replicates (calculated on the \code{\link{ilr}} scale).
#' @param  extraPredictors is an optional matrix of additional predictor values.  This should have one row per sample and one column per predictor.  The column names of the matrix will be taken as the names of the predictor.
#' The condition variable should not be included in this matrix because it is included automatically via the \code{Group} parameter.
#'
#' @details  This function is run separately for each gene and the easiest way to run it will be to follow the pipeline given in the SampleCode folder of the package.
#' In particular the easiest way to call this is to use the helper function \code{\link{startCompDTUReg}}, which loads results from \code{\link{SaveGeneLevelFiles}}.
#' These results contain all input arguments aside from \code{extraPredictors} (if any), and input arguments are automatically set by \code{\link{startCompDTUReg}}.
#' See the file (4)RunCompositionalRegressions.R in the package's SampleCode folder for example code.
#'
#' @return a data.frame containing the gene_id being used, p-value from the CompDTU or CompDTUme significance test for condition, and various information on
#' the current dataset being used.
#'
#' @export CompDTUReg
CompDTUReg <- function(genename, Y = NULL, Group, runWithME = TRUE, YInfRep = NULL, mean.withinhat = NULL,
                       extraPredictors = NULL){
  if(length(unique(Group))!=length(levels(Group))){
    stop("Check the number of levels in the Group factor, there appear to be more levels than are used and this will result in incorrect inference.")
  }

  Group2 <- Group
  nsamp <- length(Group2)
  if(runWithME==TRUE){
    ninfreps <- nrow(YInfRep)/length(Group)
    # if(is.null(ninfreps)){
    #   stop("ninfreps used must be specified when running the CompME model")
    # }
    if(is.null(YInfRep)){
      gene_id <- genename
      ret <- data.frame(gene_id, pval_CompDTUme = NA,
                        MENullCovNegVar = NA, MEAltCovNegVar = NA,
                        numYs = NA,
                        nsamp = NA, ncond = NA, stringsAsFactors = F)
      ret$numYs <- ncol(Y)
      ret$nsamp <- nsamp
      ncond <- length(unique(Group2))
      ret$ncond <- ncond

      return(ret)
    }else{
      Group2InfReps <- rep(Group2, ninfreps)

      #Repeat the extraPredictors the necessary number of times
      if(!is.null(extraPredictors)){
        #Stack the new predictor matrix the necessary number of times
        ExtraPredMatrixInfReps <- do.call(rbind, replicate(ninfreps, extraPredictors, simplify=FALSE))
        XAlt <- stats::model.matrix(~Group2InfReps + ExtraPredMatrixInfReps)

        #Reassign column names to match colunm names of extraPredictors
        colnames(XAlt)[ncol(XAlt):(ncol(XAlt) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
        XAltT <- t(XAlt)
      }else{
        XAlt <- stats::model.matrix(~Group2InfReps)
        XAltT <- t(XAlt)
      }


      ns <- nrow(YInfRep)
      bhatalt <- solve(crossprod(XAlt)) %*% (XAltT %*% YInfRep)
      pie1 <- YInfRep - (XAlt%*%bhatalt)

      SigmaTildeAltNewModeling <- (crossprod(pie1))/ns

      #resAlt <- fastLmSVB(Y = YInfRep, X = XAlt)
      #SigmaTildeAltNewModeling <- resAlt$EstCov

      if(!is.null(extraPredictors)){
        ExtraPredMatrixInfReps <- do.call(rbind, replicate(ninfreps, extraPredictors, simplify=FALSE))
        XNull <- stats::model.matrix(~ExtraPredMatrixInfReps)

        #Reassign column names to match colunm names of extraPredictors
        colnames(XNull)[ncol(XNull):(ncol(XNull) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
        XNullT <- t(XNull)
      }else{
        XNull <- stats::model.matrix(~1, data = Group2InfReps)
        XNullT <- t(XNull)
      }

      bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% YInfRep)
      pie2 <- YInfRep - (XNull%*%bhatnull)

      SigmaTildeNullNewModeling <- (crossprod(pie2))/ns

      gene_id <- genename
      ret <- data.frame(gene_id, pval_CompDTUme = NA,
                        MENullCovNegVar = NA, MEAltCovNegVar = NA,
                        numYs = NA,
                        nsamp = NA, ncond = NA, stringsAsFactors = F)
      ret$numYs <- ncol(Y)
      ret$nsamp <- nsamp
      ncond <- length(unique(Group2))
      ret$ncond <- ncond

      if(!is.null(mean.withinhat)){
        UpdatedCovAlt <- SigmaTildeAltNewModeling - mean.withinhat
        UpdatedCovNull <- SigmaTildeNullNewModeling - mean.withinhat

        statement1 <- sum(diag(UpdatedCovAlt)<=0) !=0
        if(statement1==TRUE){
          #print(paste0("UpdatedCovAlt for the new modeling approach has negative variance terms for gene ", gene_id))
          ret$MEAltCovNegVar <- TRUE
        }else{
          ret$MEAltCovNegVar <- FALSE
        }
        statement2 <- sum(diag(UpdatedCovNull)<=0) !=0
        if(statement2==TRUE){
          #print(paste0("UpdatedCovNull for the new modeling approach has negative variance terms for gene ", gene_id))
          ret$MENullCovNegVar <- TRUE
        }else{
          ret$MENullCovNegVar <- FALSE
        }

        #df_residual needs to be based on nsamp not nreps or nreps*nsamp (which is equal to nrow(YGibbs))
        if((statement1==TRUE | statement2==TRUE)){
          ret$pval_CompDTUme <- NA
        }else{
          qq <- ncond - 1
          pval_CompDTUme <- tryCatch(calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
                                                                  lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
                                                   error = function(x){})

          if(is.null(pval_CompDTUme)==FALSE){
            ret$pval_CompDTUme <- pval_CompDTUme
          }


        }

      }
      return(ret)
    }

  }else{
    if(!is.null(extraPredictors)){
      XAlt <- stats::model.matrix(~Group2 + extraPredictors)

      #Reassign column names to match colunm names of extraPredictors
      colnames(XAlt)[ncol(XAlt):(ncol(XAlt) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
      XAltT <- t(XAlt)
    }else{
      XAlt <- stats::model.matrix(~Group2)
      XAltT <- t(XAlt)
    }
    ns <- nrow(Y)
    bhatalt <- solve(crossprod(XAlt)) %*% (XAltT  %*% Y)
    pie1 <- Y - XAlt%*%bhatalt
    SigmaTildeAlt <- (crossprod(pie1))/ns

    if(!is.null(extraPredictors)){
      XNull <- stats::model.matrix(~1 + extraPredictors, data = Group2)

      #Reassign column names to match column names of extraPredictors
      colnames(XNull)[ncol(XNull):(ncol(XNull) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)

      XNullT <- t(XNull)
    }else{
      XNull <- stats::model.matrix(~1, data = Group2)
      XNullT <- t(XNull)
    }

    #bhatnull <- solve(XNullT %*% XNull) %*% XNullT %*% Y
    bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% Y)
    pie2 <- Y - XNull%*%bhatnull
    #SigmaTildeNull <- (t(pie1)%*%pie1)/ns
    SigmaTildeNull <- (crossprod(pie2))/ns



    gene_id <- genename
    ret <- data.frame(gene_id, pval_CompDTU = NA,
                      numYs = NA,
                      nsamp = NA, ncond = NA, stringsAsFactors = F)
    ret$numYs <- ncol(Y)
    ret$nsamp <- nsamp
    ncond <- length(unique(Group2))
    ret$ncond <- ncond

    qq <- ncond - 1
    pval_CompDTU <- tryCatch(calcPillaiPval(SigmaTildeNull = SigmaTildeNull, SigmaTildeAlt = SigmaTildeAlt,
                                           lm_model_fit = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
                            error = function(x){})

    if(is.null(pval_CompDTU)==FALSE){
      ret$pval_CompDTU <- pval_CompDTU
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
    return(NA)
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-stats::pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  return(pval_pillai)
}

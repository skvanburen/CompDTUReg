startDTUCompReg <- function(x, runWithME, extraPredictors = NULL){

  load(x)
  if(runWithME==TRUE){
    res <- DTUCompReg(genename = genename, Y = Y, Group = Group, runWithME = TRUE, YInfRep = YInfRep, mean.withinhat = mean.withinhat,
                      extraPredictors = extraPredictors)
  }else{
    res <- DTUCompReg(genename = genename, Y = Y, Group = Group, runWithME = FALSE, extraPredictors = extraPredictors)
  }
  return(res)
}


DTUCompReg <- function(genename, Y = NULL, Group, runWithME = TRUE, YInfRep = NULL, mean.withinhat = NULL,
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
      ret <- data.frame(gene_id, pval_pillai_ME = NA,
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
        XAlt <- model.matrix(~Group2InfReps + ExtraPredMatrixInfReps)

        #Reassign column names to match colunm names of extraPredictors
        colnames(XAlt)[ncol(XAlt):(ncol(XAlt) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
        XAltT <- t(XAlt)
      }else{
        XAlt <- model.matrix(~Group2InfReps)
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
        XNull <- model.matrix(~ExtraPredMatrixInfReps)

        #Reassign column names to match colunm names of extraPredictors
        colnames(XNull)[ncol(XNull):(ncol(XNull) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
        XNullT <- t(XNull)
      }else{
        XNull <- model.matrix(~1, data = Group2InfReps)
        XNullT <- t(XNull)
      }

      bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% YInfRep)
      pie2 <- YInfRep - (XNull%*%bhatnull)

      SigmaTildeNullNewModeling <- (crossprod(pie2))/ns

      gene_id <- genename
      ret <- data.frame(gene_id, pval_pillai_ME = NA,
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
          ret$pval_pillai_ME <- NA
        }else{
          qq <- ncond - 1
          pval_pillai_ME <- tryCatch(calcPillaiPval(SigmaTildeNull = UpdatedCovNull, SigmaTildeAlt = UpdatedCovAlt,
                                                                  res1 = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
                                                   error = function(x){})

          if(is.null(pval_pillai_ME)==FALSE){
            ret$pval_pillai_ME <- pval_pillai_ME
          }


        }

      }
      return(ret)
    }

  }else{
    if(!is.null(extraPredictors)){
      XAlt <- model.matrix(~Group2 + extraPredictors)

      #Reassign column names to match colunm names of extraPredictors
      colnames(XAlt)[ncol(XAlt):(ncol(XAlt) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)
      XAltT <- t(XAlt)
    }else{
      XAlt <- model.matrix(~Group2)
      XAltT <- t(XAlt)
    }
    ns <- nrow(Y)
    bhatalt <- solve(crossprod(XAlt)) %*% (XAltT  %*% Y)
    pie1 <- Y - XAlt%*%bhatalt
    SigmaTildeAlt <- (crossprod(pie1))/ns

    if(!is.null(extraPredictors)){
      XNull <- model.matrix(~1 + extraPredictors, data = Group2)

      #Reassign column names to match colunm names of extraPredictors
      colnames(XNull)[ncol(XNull):(ncol(XNull) - (ncol(extraPredictors) - 1))] <- colnames(extraPredictors)

      XNullT <- t(XNull)
    }else{
      XNull <- model.matrix(~1, data = Group2)
      XNullT <- t(XNull)
    }

    #bhatnull <- solve(XNullT %*% XNull) %*% XNullT %*% Y
    bhatnull <- solve(crossprod(XNull)) %*% (XNullT %*% Y)
    pie2 <- Y - XNull%*%bhatnull
    #SigmaTildeNull <- (t(pie1)%*%pie1)/ns
    SigmaTildeNull <- (crossprod(pie2))/ns



    gene_id <- genename
    ret <- data.frame(gene_id, pval_pillai = NA,
                      numYs = NA,
                      nsamp = NA, ncond = NA, stringsAsFactors = F)
    ret$numYs <- ncol(Y)
    ret$nsamp <- nsamp
    ncond <- length(unique(Group2))
    ret$ncond <- ncond

    qq <- ncond - 1
    pval_pillai <- tryCatch(calcPillaiPval(SigmaTildeNull = SigmaTildeNull, SigmaTildeAlt = SigmaTildeAlt,
                                           res1 = NULL, q = qq, nsamp = nsamp, df_residual = nsamp - ncol(XAlt)),
                            error = function(x){})

    if(is.null(pval_pillai)==FALSE){
      ret$pval_pillai <- pval_pillai
    }


  }
  return(ret)

  #End compReg
}





#q is df of condition variable
#code has been verified compared to base R ANOVA results
calcPillaiPval <- function(SigmaTildeNull, SigmaTildeAlt, res1 = NULL, q, nsamp, df_residual = NA){
  if(is.null(res1) & is.na(df_residual)){
    stop("df. residual must be specified to calcPillaiPval or else an lm object must specified in res1 to extract df.residual from")
  }

  Etilde <- nsamp * SigmaTildeAlt
  Htilde <- nsamp * (SigmaTildeNull - SigmaTildeAlt)

  vals <- diag((Htilde %*%solve(Etilde + Htilde)))
  pill_stat <- sum(vals)


  #See the Multivariate ANOVA Testing pdf document (from the SAS help file) for the necessary formulas
  #This proved to be the easiest way to calculate the statistic
  #v <- nsamp - ncol(model.matrix(res1))

  if(!is.na(df_residual)){
    v <- df_residual
  }else{
    v <- res1$df.residual
  }
  #v is the error/residual df- also extract from the r anova fit


  #p is the number of eigenvales (ie rank of Etilde + Htilde)
  p <- length(vals)

  s <- min(p,q)

  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (v - p - 1)

  #Formulas come from the SAS help file Multivariate ANOVA Testing
  piece1 <- 2*n + s + 1
  piece2 <- 2*m + s + 1
  fstat_pillai <- (piece1/piece2) * (pill_stat/(s - pill_stat))

  if(fstat_pillai < 0){
    return(NA)
  }
  numdf_pillai <- s * piece2
  denomdf_pillai <- s * piece1
  pval_pillai <- 1-pf(fstat_pillai, df1 = numdf_pillai, df2 = denomdf_pillai)
  return(pval_pillai)
}

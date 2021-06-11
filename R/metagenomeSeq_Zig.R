# To cite the original statistical method and normalization method implemented in metagenomeSeq use
#
# Paulson JN, Stine OC, Bravo HC, Pop M (2013). “Differential abundance analysis for microbial
#   marker-gene surveys.” _Nat Meth_, *advance online publication*. doi: 10.1038/nmeth.2658 (URL:
#   https://doi.org/10.1038/nmeth.2658), <URL: http://www.nature.com/nmeth/journal/vaop/ncurrent/abs/nmeth.2658.html>.
#
# To cite the metagenomeSeq software/vignette guide use
#
# Paulson JN, Olson ND, Braccia DJ, Wagner J, Talukder H, Pop M, Bravo HC (2013). _metagenomeSeq:
#   Statistical analysis for sparse high-throughput sequncing._. Bioconductor package, <URL:
#   http://www.cbcb.umd.edu/software/metagenomeSeq>.
#
# To cite time series analysis/function fitTimeSeries use
#
# Paulson* JN, Talukder* H, Bravo HC (2017). “Longitudinal differential abundance analysis of
#   marker-gene surveys using smoothing splines.” _biorxiv_. doi: 10.1101/099457 (URL:
#   https://doi.org/10.1101/099457), <URL: https://www.biorxiv.org/content/10.1101/099457v1>.


#### metagenomeSeq Functions ####
#-------------------------------#

#' Computes the weighted fold-change estimates and t-statistics.
#'
#' Wrapper to actually run the Expectation-maximization algorithm and estimate
#' $f_count$ fits.  Maximum-likelihood estimates are approximated using the EM
#' algorithm where we treat mixture membership $delta_ij = 1$ if $y_ij$ is
#' generated from the zero point mass as latent indicator variables. The
#' density is defined as $f_zig(y_ij = pi_j(S_j)*f_0(y_ij) +(1-pi_j (S_j)) *
#' f_count(y_ij; mu_i, sigma_i^2)$. The log-likelihood in this extended model
#' is: $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij) log (1-pi_j (s_j))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#'
#' @param abun Abundance table with samples as columns and features as rows
#' @param mod The model for the count distribution.
#' @param zeroMod The zero model, the model to account for the change in the
#' number of OTUs observed as a linear effect of the depth of coverage.
#' @param useCSSoffset Boolean, whether to include the default scaling
#' parameters in the model or not.
#' @param control The settings for fitZig.
#' @param useMixedModel Estimate the correlation between duplicate
#' features or replicates using duplicateCorrelation.
#' @param ... Additional parameters for duplicateCorrelation.
#'
#' @return A list of objects including:
#' \itemize{
#' 	\item{call - the call made to fitZig}
#' 	\item{fit  - 'MLArrayLM' Limma object of the weighted fit}
#' 	\item{countResiduals - standardized residuals of the fit}
#' 	\item{z - matrix of the posterior probabilities}
#' 	\item{eb - output of eBayes, moderated t-statistics, moderated F-statistics, etc}
#' 	\item{taxa - vector of the taxa names}
#' 	\item{counts - the original count matrix input}
#' 	\item{zeroMod - the zero model matrix}
#' 	\item{zeroCoef - the zero model fitted results}
#' 	\item{stillActive - convergence}
#' 	\item{stillActiveNLL - nll at convergence}
#' 	\item{dupcor - correlation of duplicates}
#' }
#' @examples
#'
#' # This is a simple demonstration
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' k = which(rowSums(MRcounts(lungTrim)>0)<30)
#' lungTrim = cumNorm(lungTrim)
#' lungTrim = lungTrim[-k,]
#' smokingStatus = pData(lungTrim)$SmokingStatus
#' mod = model.matrix(~smokingStatus)
#' # The maxit is not meant to be 1 - this is for demonstration/speed
#' settings = zigControl(maxit=1,verbose=FALSE)
#' fit = fitZig(obj = lungTrim,mod=mod,control=settings)
#'
mvfitZig <- function(abun,
                     mod,
                     zeroMod=NULL,
                     useCSSoffset=TRUE,
                     control=mvzigControl(),
                     useMixedModel=FALSE,
                   ...) {
  nf <- unname(unlist(mvcalcNormFactors(abun)))

  libsize <- colSums(abun)

  if(any(is.na(nf))) stop("At least one NA normalization factors")
  # if(any(is.na(libSize(obj)))) stop("Calculate the library size first!")

  y <- abun
  nc <- ncol(y) #nsamples
  nr <- nrow(y) #nfeatures

  # Normalization step
  Nmatrix <- log2(y + 1)

  # Initializing the model matrix
  if (useCSSoffset == TRUE){
    if (any(is.na(nf))) {
      stop("Calculate the normalization factors first!")
    }
    mmCount <- cbind(mod, log2(nf/1000 + 1))
    colnames(mmCount)[ncol(mmCount)] <- "scalingFactor"
  } else {
    mmCount <- mod
  }

  if (is.null(zeroMod)) {
    if (any(is.na(libsize))) {
      stop("Calculate the library size first!")
    }

    mmZero <- model.matrix(~1+log(libsize))
  } else {
    mmZero <- zeroMod
  }

  # dat <- mv.do_fitZig(Nmatrix, mmCount, mmZero, control=control, useMixedModel=useMixedModel, ...)
  dat <- mv.do_fitZig(Nmatrix, mmCount, mmZero, control=control, useMixedModel=useMixedModel)

  # assayData(obj)[["z"]] <- dat$z
  # assayData(obj)[["zUsed"]] <- dat$zUsed
  dat$zUsed <- NULL

  dat <- c(dat, list(call=match.call(),taxa=rownames(abun),counts=y))

  # old way of outputting results with list
  # dat <- c(dat, list(call=match.call(),taxa=rownames(obj),counts=y))

  # new output with defined results class
  dat <- new("fitZigResults", fit=dat$fit, countResiduals=dat$countResiduals,
             z=dat$z, zUsed=dat$zUsed, eb=dat$eb, zeroMod=dat$zeroMod, stillActive=dat$stillActive,
             stillActiveNLL=dat$stillActiveNLL, zeroCoef=dat$zeroCoef, dupcor=dat$dupcor, call = dat$call,
             taxa = rownames(abun), counts = dat$counts)
  dat
}

mv.do_fitZig <- function(y,
                       count_model_matrix,
                       zero_model_matrix,
                       control=zigControl(),
                       useMixedModel=FALSE,
                       ...) {
  # Initialization
  tol <- control$tol
  maxit <- control$maxit
  verbose <- control$verbose
  dfMethod <- control$dfMethod
  pvalMethod <- control$pvalMethod

  nr <- nrow(y)
  nc <- ncol(y)

  zeroIndices <- (y == 0)
  z <- matrix(0, nrow=nr, ncol=nc)
  z[zeroIndices] <- 0.5
  zUsed <- z

  curIt <- 0
  nllOld <- rep(Inf, nr)
  nll <- rep(Inf, nr)
  nllUSED <- nll
  stillActive <- rep(TRUE, nr)
  stillActiveNLL <- rep(1, nr)
  dupcor <- NULL

  modRank <- ncol(count_model_matrix)
  # E-M Algorithm
  while (any(stillActive) && (curIt < maxit)) {

    # M-step for count density (each feature independently)
    if(curIt == 0){
      fit <- mvdoCountMStep(z, y, count_model_matrix, stillActive, dfMethod=dfMethod)
    } else {
      fit <- mvdoCountMStep(z, y, count_model_matrix, stillActive, fit2=fit, dfMethod=dfMethod)
    }

    # M-step for zero density (all features together)
    zeroCoef <- mvdoZeroMStep(z, zeroIndices, zero_model_matrix)

    # E-step
    z <- mvdoEStep(fit$residuals, zeroCoef$residuals, zeroIndices)
    zzdata <- mvgetZ(z, zUsed, stillActive, nll, nllUSED);
    zUsed <- zzdata$zUsed;

    # NLL
    nll <- mvgetNegativeLogLikelihoods(z, fit$residuals, zeroCoef$residuals)
    eps <- mvgetEpsilon(nll, nllOld)
    active <- mvisItStillActive(eps, tol,stillActive,stillActiveNLL,nll)
    stillActive <- active$stillActive;
    stillActiveNLL <- active$stillActiveNLL;
    if (verbose == TRUE){
      cat(sprintf("it=%2d, nll=%0.2f, log10(eps+1)=%0.2f, stillActive=%d\n", curIt, mean(nll,na.rm=TRUE), log10(max(eps,na.rm=TRUE)+1), sum(stillActive)))
    }
    nllOld <- nll
    curIt <- curIt + 1

    if (sum(rowSums((1-z) > 0) <= modRank, na.rm=TRUE) > 0) {
      k <- which(rowSums((1-z) > 0) <= modRank)
      stillActive[k] <- FALSE;
      stillActiveNLL[k] <- nll[k]
    }
  }

  if (useMixedModel == TRUE) {
    dupcor <- duplicateCorrelation(y, count_model_matrix, weights=(1-z), ...)
    fit$fit <- limma::lmFit(y, count_model_matrix, weights=(1-z), correlation=dupcor$consensus, ...)
    countCoef <- fit$fit$coefficients
    countMu <- tcrossprod(countCoef, count_model_matrix)
    fit$residuals <- sweep((y-countMu), 1, fit$fit$sigma, "/")
  }

  eb <- limma::eBayes(fit$fit)
  dat <- list(fit=fit$fit, countResiduals=fit$residuals,
              z=z, zUsed=zUsed, eb=eb, zeroMod=zero_model_matrix, stillActive=stillActive,
              stillActiveNLL=stillActiveNLL, zeroCoef=zeroCoef, dupcor=dupcor)
  dat
}

#' Settings for the fitZig function
#'
#' @param tol The tolerance for the difference in negative log likelihood estimates for a feature to remain active.
#' @param maxit The maximum number of iterations for the expectation-maximization algorithm.
#' @param verbose Whether to display iterative step summary statistics or not.
#' @param dfMethod Either 'default' or 'modified' (by responsibilities).
#' @param pvalMethod Either 'default' or 'bootstrap'.
#' @return The value for the tolerance, maximum no. of iterations, and the verbose warning.
#'
#' @aliases settings2
#' @examples
#' control =  zigControl(tol=1e-10,maxit=10,verbose=FALSE)
#'
mvzigControl <- function(tol=1e-4,maxit=10,verbose=TRUE,dfMethod="modified",pvalMethod="default"){
  # to do: add stop if not
  DFMETHODS <- c("default", "modified")
  PMETHODS  <- c("default", "bootstrap")
  dfMethod  <- DFMETHODS[pmatch(dfMethod, DFMETHODS)]
  pvalMethod<- PMETHODS[pmatch(pvalMethod,PMETHODS)]

  stopifnot(dfMethod%in%DFMETHODS)
  stopifnot(pvalMethod%in%PMETHODS)

  set <-list(tol=tol,maxit=maxit,verbose=verbose,dfMethod=dfMethod,pvalMethod=pvalMethod);
  return(set)
}

#' Compute the Maximization step calculation for features still active.
#'
#' Maximization step is solved by weighted least squares.  The function also
#' computes counts residuals.
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The density is defined as
#' $f_zig(y_ij = pi_j(S_j)*f_0(y_ij) +(1-pi_j (S_j)) *
#' f_count(y_ij;mu_i,sigma_i^2)$. The log-likelihood in this extended model is
#' $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (s_j))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a
#' count comes from a spike distribution at 0).
#' @param y Matrix (m x n) of count observations.
#' @param mmCount Model matrix for the count distribution.
#' @param stillActive Boolean vector of size M, indicating whether a feature
#' converged or not.
#' @param fit2 Previous fit of the count model.
#' @param dfMethod Either 'default' or 'modified' (by responsibilities)
#' @return Update matrix (m x n) of estimate responsibilities (probabilities
#' that a count comes from a spike distribution at 0).
mvdoCountMStep <- function(z, y, mmCount, stillActive,fit2=NULL,dfMethod="modified") {
  if (is.null(fit2)){
    fit=limma::lmFit(y[stillActive,],mmCount,weights = (1-z[stillActive,]))
    if(dfMethod=="modified"){
      df = rowSums(1-z[stillActive,,drop=FALSE]) - ncol(mmCount)
      fit$df[stillActive] = df
      fit$df.residual[stillActive] = df
    }
    countCoef = fit$coefficients
    countMu=tcrossprod(countCoef, mmCount)
    residuals=sweep((y[stillActive,,drop=FALSE]-countMu),1,fit$sigma,"/")
    dat = list(fit = fit, residuals = residuals)
    return(dat)
  } else {

    residuals = fit2$residuals
    fit2 = fit2$fit

    fit=limma::lmFit(y[stillActive,,drop=FALSE],mmCount,weights = (1-z[stillActive,,drop=FALSE]))

    fit2$coefficients[stillActive,] = fit$coefficients
    fit2$stdev.unscaled[stillActive,]=fit$stdev.unscaled
    fit2$sigma[stillActive] = fit$sigma
    fit2$Amean[stillActive] = fit$Amean

    if(dfMethod=="modified"){
      df = rowSums(1-z[stillActive,,drop=FALSE]) - ncol(mmCount)
      fit$df = df
      fit$df.residual = df
    }
    fit2$df[stillActive]    = fit$df
    fit2$df.residual[stillActive]    = fit$df.residual

    countCoef = fit$coefficients
    countMu=tcrossprod(countCoef, mmCount)
    r=sweep((y[stillActive,,drop=FALSE]-countMu),1,fit$sigma,"/")
    residuals[stillActive,]=r

    dat = list(fit = fit2, residuals=residuals)

    return(dat)
  }
}

#' Compute the zero Maximization step.
#'
#' Performs Maximization step calculation for the mixture components. Uses
#' least squares to fit the parameters of the mean of the logistic
#' distribution. $$ pi_j = sum_i^M frac1Mz_ij $$ Maximum-likelihood estimates
#' are approximated using the EM algorithm where we treat mixture membership
#' $delta_ij$ = 1 if $y_ij$ is generated from the zero point mass as latent
#' indicator variables. The density is defined as $f_zig(y_ij = pi_j(S_j) cdot
#' f_0(y_ij) +(1-pi_j (S_j))cdot f_count(y_ij;mu_i,sigma_i^2)$. The
#' log-likelihood in this extended model is $(1-delta_ij) log
#' f_count(y;mu_i,sigma_i^2 )+delta_ij log pi_j(s_j)+(1-delta_ij)log (1-pi_j
#' (sj))$. The responsibilities are defined as $z_ij = pr(delta_ij=1 | data)$.
#'
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a
#' count comes from a spike distribution at 0).
#' @param zeroIndices Index (matrix m x n) of counts that are zero/non-zero.
#' @param mmZero The zero model, the model matrix to account for the change in
#' the number of OTUs observed as a linear effect of the depth of coverage.
#' @return List of the zero fit (zero mean model) coefficients, variance -
#' scale parameter (scalar), and normalized residuals of length
#' sum(zeroIndices).
mvdoZeroMStep <- function(z, zeroIndices, mmZero) {
  pi=sapply(1:ncol(zeroIndices), function(j) {
    if (sum(zeroIndices[,j])==0){
      return(1e-8)
    }
    tmp=mean(z[zeroIndices[,j],j],na.rm=TRUE)
    ifelse(tmp<=1e-8, 1e-8, ifelse(tmp>=1-(1e-8),1-(1e-8),tmp))
  })
  zeroLM=lm.fit(mmZero, qlogis(pi))
  zeroCoef=zeroLM$coef

  r=zeroLM$residuals
  sigma=sd(r)+(1e-3)

  list(zeroLM=zeroLM, zeroCoef=zeroCoef, sigma=sigma, residuals=r/sigma)
}

#' Compute the Expectation step.
#'
#' Estimates the responsibilities $z_ij = fracpi_j cdot I_0(y_ijpi_j cdot
#' I_0(y_ij + (1-pi_j) cdot f_count(y_ij
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The density is defined as
#' $f_zig(y_ij = pi_j(S_j) cdot f_0(y_ij) +(1-pi_j (S_j))cdot
#' f_count(y_ij;mu_i,sigma_i^2)$. The log-likelihood in this extended model is
#' $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (sj))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#'
#' @param countResiduals Residuals from the count model.
#' @param zeroResiduals Residuals from the zero model.
#' @param zeroIndices Index (matrix m x n) of counts that are zero/non-zero.
#' @return Updated matrix (m x n) of estimate responsibilities (probabilities
#' that a count comes from a spike distribution at 0).
mvdoEStep <- function(countResiduals,  zeroResiduals, zeroIndices) {
  pi_prop=mvgetPi(zeroResiduals)
  w1=sweep(zeroIndices, 2, pi_prop, FUN="*")

  countDensity=mvgetCountDensity(countResiduals)
  w2=sweep(countDensity, 2, 1-pi_prop, FUN="*")
  z=w1/(w1+w2)
  z[z>1-1e-6]=1-1e-6
  z[!zeroIndices]=0
  z
}

#' Calculate the current Z estimate responsibilities (posterior probabilities)
#'
#' Calculate the current Z estimate responsibilities (posterior probabilities)
#'
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a
#' count comes from a spike distribution at 0).
#' @param zUsed Matrix (m x n) of estimate responsibilities (probabilities that
#' a count comes from a spike distribution at 0) that are actually used
#' (following convergence).
#' @param stillActive A vector of size M booleans saying if a feature is still
#' active or not.
#' @param nll Vector of size M with the current negative log-likelihoods.
#' @param nllUSED Vector of size M with the converged negative log-likelihoods.
#' @return A list of updated zUsed and nllUSED.
mvgetZ <- function(z,zUsed,stillActive,nll,nllUSED) {
  nllUSED[stillActive] = nll[stillActive]
  k =which(nll< (nllUSED))
  if(length(k)>0){
    zUsed[k,]=z[k,]
    nllUSED[k] = nll[k]
  }
  zUsed[stillActive,] = z[stillActive,]
  dat = list(zUsed = zUsed,nllUSED = nllUSED)
  return(dat);
}

#' Calculate the negative log-likelihoods for the various features given the
#' residuals.
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The log-likelihood in this
#' extended model is $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (sj))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data and current values)$.
#'
#'
#' @param z Matrix (m x n) of estimate responsibilities (probabilities that a
#' count comes from a spike distribution at 0).
#' @param countResiduals Residuals from the count model.
#' @param zeroResiduals Residuals from the zero model.
#' @return Vector of size M of the negative log-likelihoods for the various
#' features.
mvgetNegativeLogLikelihoods <- function(z, countResiduals, zeroResiduals) {
  pi=mvgetPi(zeroResiduals)
  countDensity=mvgetCountDensity(countResiduals, log=TRUE)
  res=(1-z) * countDensity
  res=res+sweep(z, 2, log(pi), FUN="*")
  res=res+sweep(1-z,2,log(1-pi), FUN="*")
  -rowSums(res)
}

#' Calculate the relative difference between iterations of the negative
#' log-likelihoods.
#'
#' Maximum-likelihood estimates are approximated using the EM algorithm where
#' we treat mixture membership $delta_ij$ = 1 if $y_ij$ is generated from the
#' zero point mass as latent indicator variables. The log-likelihood in this
#' extended model is $(1-delta_ij) log f_count(y;mu_i,sigma_i^2 )+delta_ij log
#' pi_j(s_j)+(1-delta_ij)log (1-pi_j (sj))$. The responsibilities are defined
#' as $z_ij = pr(delta_ij=1 | data)$.
#'
#'
#' @param nll Vector of size M with the current negative log-likelihoods.
#' @param nllOld Vector of size M with the previous iterations negative
#' log-likelihoods.
#' @return Vector of size M of the relative differences between the previous
#' and current iteration nll.
mvgetEpsilon <- function(nll, nllOld) {
  eps=(nllOld-nll)/nllOld
  ifelse(!is.finite(nllOld), Inf, eps)
}

#' Function to determine if a feature is still active.
#'
#' In the Expectation Maximization routine features posterior probabilities routinely converge based on a tolerance threshold. This function checks
#' whether or not the feature's negative log-likelihood (measure of the fit) has changed or not.
#'
#' @param eps Vector of size M (features) representing the relative difference between the new nll and old nll.
#' @param tol The threshold tolerance for the difference
#' @param stillActive A vector of size M booleans saying if a feature is still active or not.
#' @param stillActiveNLL A vector of size M recording the negative log-likelihoods of the various features, updated for those still active.
#' @param nll Vector of size M with the current negative log-likelihoods.
#' @return None.
#'
#'
mvisItStillActive <- function(eps, tol,stillActive,stillActiveNLL,nll) {
  stillActive[stillActive]=!is.finite(eps[stillActive]) | eps[stillActive]>tol
  stillActive[which(is.na(eps))]=FALSE

  stillActiveNLL[stillActive]=nll[stillActive]
  dat = list(stillActive=stillActive,stillActiveNLL = stillActiveNLL)
  return(dat)
}

#' Calculate the mixture proportions from the zero model / spike mass model
#' residuals.
#'
#' F(x) = 1 / (1 + exp(-(x-m)/s)) (the CDF of the logistic distribution).
#' Provides the probability that a real-valued random variable X with a given
#' probability distribution will be found at a value less than or equal to x.
#' The output are the mixture proportions for the samples given the residuals
#' from the zero model.
#'
#'
#' @param residuals Residuals from the zero model.
#' @return Mixture proportions for each sample.
mvgetPi <-  function(residuals) {
  plogis(residuals)
}

#' Compute the value of the count density function from the count model
#' residuals.
#'
#' Calculate density values from a normal: $f(x) = 1/(sqrt (2 pi ) sigma )
#' e^-((x - mu )^2/(2 sigma^2))$.  Maximum-likelihood estimates are
#' approximated using the EM algorithm where we treat mixture membership
#' $deta_ij$ = 1 if $y_ij$ is generated from the zero point mass as latent
#' indicator variables. The density is defined as $f_zig(y_ij = pi_j(S_j) cdot
#' f_0(y_ij) +(1-pi_j (S_j))cdot f_count(y_ij;mu_i,sigma_i^2)$. The
#' log-likelihood in this extended model is $(1-delta_ij) log
#' f_count(y;mu_i,sigma_i^2 )+delta_ij log pi_j(s_j)+(1-delta_ij)log (1-pi_j
#' (sj))$. The responsibilities are defined as $z_ij = pr(delta_ij=1 | data)$.
#'
#'
#' @param residuals Residuals from the count model.
#' @param log Whether or not we are calculating from a log-normal distribution.
#' @return Density values from the count model residuals.
mvgetCountDensity <- function(residuals, log=FALSE) {
  dnorm(residuals,log=log)
}

#' Cumulative sum scaling percentile selection
#'
#' Calculates the percentile for which to sum counts up to and scale by.
#' cumNormStat might be deprecated one day. Deviates from methods in Nature Methods paper
#' by making use row means for generating reference.
#'
#' @param obj A matrix or MRexperiment object.
#' @param qFlag Flag to either calculate the proper percentile using
#' R's step-wise quantile function or approximate function.
#' @param pFlag Plot the relative difference of the median deviance from the reference.
#' @param rel Cutoff for the relative difference from one median difference
#' from the reference to the next
#' @param ... Applicable if pFlag == TRUE. Additional plotting parameters.
#' @return Percentile for which to scale data
#' @examples
#'
#' data(mouseData)
#' p = round(cumNormStat(mouseData,pFlag=FALSE),digits=2)
#'
mvcumNormStat <- function(mat,qFlag = TRUE,pFlag = FALSE,rel=.1,...){
  # mat = returnAppropriateObj(obj,FALSE,FALSE)
  if(any(colSums(mat)==0)) stop("Warning empty sample")

  smat = sapply(1:ncol(mat),function(i){sort(mat[,i],decreasing=FALSE)})
  ref  = rowMeans(smat);

  yy = mat;
  yy[yy==0]=NA;

  ncols = ncol(mat);
  refS = sort(ref);

  k = which(refS>0)[1]
  lo = (length(refS)-k+1)

  if(qFlag == TRUE){
    diffr = sapply(1:ncols,function(i){
      refS[k:length(refS)] - quantile(yy[,i],p=seq(0,1,length.out=lo),na.rm=TRUE)
    })
  }
  if(qFlag == FALSE){
    diffr = sapply(1:ncols,function(i){
      refS[k:length(refS)] - approx(sort(yy[,i],decreasing=FALSE),n=lo)$y
    })
  }
  diffr2 = matrixStats::rowMedians(abs(diffr),na.rm=TRUE)
  if(pFlag ==TRUE){
    plot(abs(diff(diffr2[diffr2>0]))/diffr2[diffr2>0][-1],type="h",ylab="Relative difference for reference",xaxt="n",...)
    abline(h=rel)
    axis(1,at=seq(0,length(diffr2),length.out=5),labels = seq(0,1,length.out=5))
  }
  x = which(abs(diff(diffr2))/diffr2[-1]>rel)[1] / length(diffr2)
  if(x<=0.50){
    message("Default value being used.")
    x = 0.50
  }
  if(class(obj)=="MRexperiment"){
    obj@expSummary$cumNormStat = x;
  }
  return(x)
}

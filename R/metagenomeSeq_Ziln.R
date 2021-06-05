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

#' Computes differential abundance analysis using a zero-inflated log-normal model
#'
#' Wrapper to actually run zero-inflated log-normal model given a MRexperiment object
#' and model matrix. User can decide to shrink parameter estimates.
#'
#' @param abun Abundance table with samples as columns and features as rows
#' @param mod The model for the count distribution.
#' @param coef Coefficient of interest to grab log fold-changes.
#' @param B Number of bootstraps to perform if >1. If >1 performs permutation test.
#' @param szero TRUE/FALSE, shrink zero component parameters.
#' @param spos TRUE/FALSE, shrink positive component parameters.
#'
#' @return A list of objects including:
#' \itemize{
#'  \item{call - the call made to fitFeatureModel}
#'  \item{fitZeroLogNormal  - list of parameter estimates for the zero-inflated log normal model}
#'  \item{design - model matrix}
#'  \item{taxa - taxa names}
#'  \item{counts - count matrix}
#'  \item{pvalues - calculated p-values}
#'  \item{permuttedfits - permutted z-score estimates under the null}
#' }
#' @examples
#'
#' data(lungData)
#' lungData = lungData[,-which(is.na(pData(lungData)$SmokingStatus))]
#' lungData=filterData(lungData,present=30,depth=1)
#' lungData <- cumNorm(lungData, p=.5)
#' s <- normFactors(lungData)
#' pd <- pData(lungData)
#' mod <- model.matrix(~1+SmokingStatus, data=pd)
#' lungres1 = fitFeatureModel(lungData,mod)
#'
mvfitFeatureModel <- function(abun,mod,coef=2,B=1,szero=FALSE,spos=TRUE) {
  nf <- unname(unlist(mvcalcNormFactors(abun)))
  if(any(is.na(nf))) stop("At least one NA normalization factors")

  mmCount <- cbind(mod, log(nf/median(nf)))
  colnames(mmCount)[ncol(mmCount)] = "scalingFactor"

  if(ncol(mmCount)>3){ stop("Can't analyze currently.") }
  i = permuttedFits = NULL

  # These pieces get to be a part of the new zero-ln model!
  fitzeroln = mvfitZeroLogNormal(abun,mmCount,coef=coef,szero=szero,spos=spos)

  if(any(is.na(fitzeroln$logFC))) {
    feats = which(is.na(fitzeroln$logFC))
    mat = MRcounts(obj[feats,], norm=TRUE, log=FALSE,sl=median(nf))
    fit = lmFit(log(mat+1),mmCount)
    fit = eBayes(fit)
    fitzeroln$logFC[feats] = coefficients(fit)[,coef]
    fitzeroln$se[feats] = (sqrt(fit$s2.post)*fit$stdev.unscaled)[,coef]
  }
  zscore = fitzeroln$logFC/fitzeroln$se

  if(B>1){
    permutations = replicate(B,sample(mmCount[,coef]))
    mmCountPerm  = mmCount

    permuttedFits = foreach(i = seq(B),.errorhandling="remove",
                            .packages=c("metagenomeSeq","glmnet")) %dopar% {
                              mmCountPerm[,coef] = permutations[,i]
                              permFit = fitZeroLogNormal(obj,mmCountPerm,coef=coef,szero=szero,spos=spos)
                              permFit$logFC/permFit$se
                            }
    zperm = abs(sapply(permuttedFits,function(i)i))
    pvals = rowMeans(zperm>=abs(zscore),na.rm=TRUE)
  } else {
    pvals = 2*(1-pnorm(abs(zscore)))
  }
  # old way of creating results object
  # res = list(call=match.call(),fitZeroLogNormal=fitzeroln,design=mmCount,
  #   taxa=rownames(obj),counts=MRcounts(obj),pvalues=pvals,permuttedFits=permuttedFits)

  # new way with defined results class
  res = new("fitFeatureModelResults", call = match.call(), fitZeroLogNormal=fitzeroln,
            design = mmCount, taxa = rownames(abun), counts = abun,
            pvalues = pvals, permuttedFits = permuttedFits)
  res
}

#' Cumulative sum scaling (css) normalization factors
#'
#' Return a vector of the the sum up to and including a quantile.
#'
#' @param x Abundance table with samples as columns and features as rows
#' @param p The pth quantile.
#'
#' @return Vector of the sum up to and including a sample's pth quantile.
#' @examples
#'
#' data(mouseData)
#' head(calcNormFactors(mouseData))
#'
mvcalcNormFactors <- function(x,p=mvcumNormStatFast(x)){
  # abun must be a matrix with samples as columns and taxa/features as rows
  xx = x
  xx[x == 0] <- NA
  qs = matrixStats::colQuantiles(xx, probs = p, na.rm = TRUE)
  normFactors <- sapply(1:ncol(xx), function(i) {
    xx = (x[, i] - .Machine$double.eps)
    sum(xx[xx <= qs[i]])
  })
  names(normFactors)<-colnames(x)
  as.data.frame(normFactors)
}

#' Cumulative sum scaling percentile selection
#'
#' Calculates the percentile for which to sum counts up to and scale by. Faster
#' version than available in cumNormStat. Deviates from methods described in Nature Methods by
#' making use of ro means for reference.
#'
#' @param mat Abundance table with samples as columns and features as rows
#' @param pFlag Plot the median difference quantiles.
#' @param rel Cutoff for the relative difference from one median difference
#' from the reference to the next.
#' @param ... Applicable if pFlag == TRUE. Additional plotting parameters.
#'
#' @return Percentile for which to scale data
#' @examples
#'
#' data(mouseData)
#' p = round(cumNormStatFast(mouseData,pFlag=FALSE),digits=2)
#'
mvcumNormStatFast <-function(mat,pFlag = FALSE,rel=.1,...){
  smat = lapply(1:ncol(mat), function(i) {
    sort(mat[which(mat[, i]>0),i], decreasing = TRUE)
  })
  leng = max(sapply(smat,length))
  if(any(sapply(smat,length)==1)) stop("Warning sample with one or zero features")

  smat2 = array(NA,dim=c(leng,ncol(mat)))
  for(i in 1:ncol(mat)){
    smat2[leng:(leng-length(smat[[i]])+1),i] = smat[[i]]
  }

  rmat2 = sapply(1:ncol(smat2),function(i){
    quantile(smat2[,i],p=seq(0,1,length.out=nrow(smat2)),na.rm=TRUE)
  })
  smat2[is.na(smat2)] = 0
  ref1 = rowMeans(smat2)

  ncols = ncol(rmat2)
  diffr = sapply(1:ncols, function(i) {
    ref1 - rmat2[,i]
  })
  diffr1=matrixStats::rowMedians(abs(diffr))
  if(pFlag==TRUE){
    plot(abs(diff(diffr1))/diffr1[-1],type="h",...)
    abline(h=rel)
    axis(1,at=seq(0,length(diffr1),length.out=5),labels = seq(0,1,length.out=5))
  }
  x= which(abs(diff(diffr1))/diffr1[-1] > rel)[1]/length(diffr1)
  if(x<=0.50){
    message("Default value being used.")
    x = 0.50
  }

  return(x)
}

#' Compute the log fold-change estimates for the zero-inflated log-normal model
#'
#' Run the zero-inflated log-normal model given a MRexperiment object
#' and model matrix. Not for the average user, assumes structure of the model matrix.
#'
#' @param abun Abundance table with samples as columns and features as rows
#' @param mod The model for the count distribution.
#' @param coef Coefficient of interest to grab log fold-changes.
#' @param szero TRUE/FALSE, shrink zero component parameters.
#' @param spos TRUE/FALSE, shrink positive component parameters.
#'
#' @return A list of objects including:
#' \itemize{
#'  \item{logFC - the log fold-change estimates}
#'  \item{adjFactor  - the adjustment factor based on the zero component}
#'  \item{se - standard error estimates}
#'  \item{fitln - parameters from the log-normal fit}
#'  \item{fitzero - parameters from the logistic fit}
#'  \item{zeroRidge - output from the ridge regression}
#'  \item{posRidge - output from the ridge regression}
#'  \item{tauPos - estimated tau^2 for positive component}
#'  \item{tauZero - estimated tau^2 for zero component}
#'  \item{exclude - features to exclude for various reasons, e.g. all zeros}
#'  \item{zeroExclude - features to exclude for various reasons, e.g. all zeros}
#' }
mvfitZeroLogNormal<-function(abun,mod,coef=2,szero=TRUE,spos=TRUE){
  positiveMod = mod[,-ncol(mod)]
  zeroMod = mod

  nf <- unname(unlist(mvcalcNormFactors(abun)))

  mat <- sweep(abun,2,nf/median(nf),'/')

  posIndices = mat>0

  nr = nrow(mat)
  nc = ncol(mat)
  exclude = zeroExclude = tauZero = tauPos = posRidge = zeroRidge = NULL

  results = array(NA,dim=c(nr,3))
  rownames(results) = rownames(mat)
  colnames(results) = c("logFC","adjFactor","se")

  # calc log-normal component
  fitln = mvcalcPosComponent(mat,positiveMod,posIndices)

  # Don't calculate shrinkage with special cases
  zeros2 = which(fitln[,"s2"]==0)
  rs = rowsum(t(1-(1-posIndices)),positiveMod[,coef])
  exclude = union(which(rs[1,]<=1),which(rs[2,]<=1))
  zeroExclude  = which(colSums(rs)>=(nc-3))
  exclude = union(zeros2,exclude); if(length(exclude)==0) exclude=NULL
  if(length(zeroExclude)==0) zeroExclude=NULL

  sdensity = density(fitln[,"s2"],na.rm=TRUE)
  smode = sdensity$x[which.max(sdensity$y)]
  if(length(zeros2)>0) fitln[zeros2,"s2"] = smode

  # shrink positive
  if(spos==TRUE){
    shrinkPos<-mvcalcShrinkParameters(fitln,coef,smode,exclude)
    tauPos = shrinkPos$tau
    vpost = shrinkPos$v.post
    fitln[,"s2"] = vpost

    posRidge = sapply(seq(nr),function(i){
      k = which(posIndices[i,])
      y = log(mat[i,k])
      x = positiveMod[k,]
      l = vpost[i]/(nrow(x)*tauPos)
      if(i %in% exclude) return(matrix(rep(NA,ncol(positiveMod))))
      ridge = glmnet(y=y,x=x,lambda=l,alpha=0)
      as.matrix(coefficients(ridge)[colnames(positiveMod),])
    })
    posFittedCoefficients = t(posRidge)
    rownames(posFittedCoefficients) = rownames(mat)
    fitln[rownames(posFittedCoefficients),1:ncol(positiveMod)] = posFittedCoefficients
  }
  # calc zero component
  fitzero=mvcalcZeroComponent(mat,zeroMod,posIndices)

  sdensity = density(fitzero[,"s2"],na.rm=TRUE)
  smode = sdensity$x[which.max(sdensity$y)]
  if(length(exclude)>0) fitzero[exclude,"s2"] = smode

  # shrink zero
  if(szero==TRUE){
    shrinkZero<-calcShrinkParameters(fitzero,coef,smode,exclude)
    tauZero = shrinkZero$tau
    vpostZero = shrinkZero$v.post
    fitzero[,"s2"] = vpostZero

    zeroRidge = sapply(1:nr,function(i){
      y = posIndices[i,]
      l = 1/(nc*tauZero)
      if(i %in% c(zeroExclude,exclude)) return(matrix(rep(NA,ncol(zeroMod))))
      ridge = glmnet(y=y,x=zeroMod,lambda=l,family="binomial",alpha=0,
                     penalty.factor = c(rep(1,(ncol(zeroMod)-1)),0))
      as.matrix(coefficients(ridge))[colnames(zeroMod),]
    })
    zeroFittedCoefficients = t(zeroRidge)
    rownames(zeroFittedCoefficients) = rownames(mat)
    fitzero[rownames(zeroFittedCoefficients),1:ncol(zeroMod)] = zeroFittedCoefficients
  }

  # calc se
  se = mvcalcStandardError(zeroMod,fitln,fitzero,coef=coef,exclude=union(exclude,zeroExclude))
  se[zeroExclude] = sqrt(fitln[zeroExclude,"s2"])

  # calc adjFactor
  adjFactor = mvcalcZeroAdjustment(fitln,fitzero,zeroMod,coef,exclude=exclude)
  adjFactor[zeroExclude] = 0

  # calc logFC
  logFC <- fitln[,coef] + adjFactor

  list(logFC=logFC,adjFactor=adjFactor,se=se,
       fitln=fitln,fitzero=fitzero,zeroRidge=zeroRidge,posRidge=posRidge,
       tauPos=tauPos,tauZero=tauZero,exclude=exclude,zeroExclude=zeroExclude)
}
#' Positive component
#'
#' Fit the positive (log-normal) component
#'
#' @param mat A matrix of normalized counts
#' @param mod A model matrix
#' @param weights Weight matrix for samples and counts
mvcalcPosComponent<-function(mat,mod,weights){
  fitln <- lmFit(log(mat),mod,weights=weights)
  b = coefficients(fitln)
  df = fitln$df
  res = residuals(fitln,log(mat))
  s2 = sapply(seq(nrow(res)),function(i){
    sum(res[i,which(weights[i,])]^2,na.rm=TRUE)/df[i]
  })
  fitln<-data.frame(b=b,s2=s2,df=df)
  rownames(fitln) = rownames(mat)
  fitln
}
#' Zero component
#'
#' Fit the zero (logisitic) component
#'
#' @param mat A matrix of normalized counts
#' @param mod A model matrix
#' @param weights Weight matrix for samples and counts
mvcalcZeroComponent<-function(mat,mod,weights){
  fitzero <- sapply(seq(nrow(mat)), function(i) {
    fit <- glm.fit(mod, weights[i,], family=binomial())
    cf = coefficients(fit)
    df = fit$df.residual
    mc = exp(mod %*% cf)
    s2 = sum((weights[i, ] - t(mc/(1 + mc)))^2)/df
    # s2 = sum(residuals(fit)^2)/df
    c(beta= cf, s2 = s2, df = df)
  })
  fitzero <- data.frame(t(fitzero))
  rownames(fitzero) = rownames(mat)
  fitzero
}
#' Calculate shrinkage parameters
#'
#' Calculate the shrunken variances and variance of parameters of interest across features.
#'
#' @param fit A matrix of fits as outputted by calcZeroComponent or calcPosComponent
#' @param coef Coefficient of interest
#' @param mins2 minimum variance estimate
#' @param exclude Vector of features to exclude when shrinking
mvcalcShrinkParameters<-function(fit,coef,mins2,exclude=NULL){

  if(is.null(exclude)){
    shrunkVar <- limma::squeezeVar(fit[,"s2"], fit[,"df"])
    v.post = shrunkVar$var.post
    tau <-var(fit[,coef],na.rm=TRUE)
  } else {
    v.post = rep(mins2,nrow(fit))
    shrunkVar <- limma::squeezeVar(fit[-exclude,"s2"], fit[-exclude,"df"])
    v.post[-exclude] <- shrunkVar$var.post
    tau <- var(fit[-exclude,coef],na.rm=TRUE)
  }
  list(tau=tau,v.post=v.post)
}
#' Calculate the zero-inflated component's adjustment factor
#'
#' Calculate the log ratio of average marginal probabilities for each sample
#' having a positive count. This becomes the adjustment factor for the log
#' fold change.
#'
#' @param fitln A matrix with parameters from the log-normal fit
#' @param fitzero A matrix with parameters from the logistic fit
#' @param mod The zero component model matrix
#' @param coef Coefficient of interest
#' @param exclude List of features to exclude
mvcalcZeroAdjustment<-function(fitln,fitzero,mod,coef,exclude=NULL){
  b = fitln[,1:(ncol(mod)-1)]
  beta = fitzero[,1:ncol(mod)]
  # calculate for zero adjust factor
  mod1 <- mod
  mod1[,coef] <- 1
  theta1 <- mod1 %*% t(beta)
  p1 <- exp(theta1) / (1+exp(theta1))
  p1 <- t(p1)
  if(ncol(b)>2) p1 = p1*exp(t(mod[,3:(ncol(mod)-1)]%*%t(b[,3:ncol(b)])))
  mean_p1 <- rowMeans(p1)

  mod0 <- mod
  mod0[,coef] <- 0
  theta0 <- mod0 %*% t(beta)
  p0 <- exp(theta0) / (1+exp(theta0))
  p0 <- t(p0)
  if(ncol(b)>2) p0 = p0*exp(t(mod[,3:(ncol(mod)-1)]%*%t(b[,3:ncol(b)])))
  mean_p0 <- rowMeans(p0)

  adjFactor <- log(mean_p1/mean_p0)
  if(!is.null(exclude)) adjFactor[exclude] = NA
  adjFactor
}

#' Calculate the zero-inflated log-normal statistic's standard error
#'
#' Calculat the se for the model. Code modified from
#' "Adjusting for covariates in zero-inflated gamma and
#' zero-inflated log-normal models for semicontinuous data", ED Mills
#'
#' @param mod The zero component model matrix
#' @param fitln A matrix with parameters from the log-normal fit
#' @param fitzero A matrix with parameters from the logistic fit
#' @param coef Coefficient of interest
#' @param exclude List of features to exclude
mvcalcStandardError<-function(mod,fitln,fitzero,coef=2,exclude=NULL){
  mod0 = mod1 = mod
  mod1[,coef] <- 1
  mod0[,coef] <- 0
  ve = rep(NA,nrow(fitln))
  features = seq(nrow(fitln))
  if(length(exclude)>0) features = features[-exclude]

  # a) need to speed up
  # b) need to include more covariates

  fullvar = sapply(features,function(i){
    beta = fitzero[i,1:ncol(mod)]
    b = fitln[i,1:(ncol(mod)-1)]
    s = as.numeric(fitln[i,"s2"])

    mu0 = as.vector(exp(mod0[,-ncol(mod)]%*%t(b) + .5*s))
    mu1 = as.vector(exp(mod1[,-ncol(mod)]%*%t(b) + .5*s))

    # calculate for zero adjust factor
    theta <- mod %*% t(beta)
    theta1 <- mod1 %*% t(beta)
    theta0 <- mod0 %*% t(beta)
    p  <- t(exp(theta) / (1+exp(theta)))
    p1 <- t(exp(theta1) / (1+exp(theta1)))
    p0 <- t(exp(theta0) / (1+exp(theta0)))

    checkInverse <- function(m){
      inherits(try(qr.solve(m),silent=T), "matrix")
    }

    Dp2 <- diag(length(p))*as.vector(p*(1-p))
    infz = t(mod)%*%Dp2%*%mod
    Dp <- diag(length(p))*as.vector(p)
    infln = t(mod[,-ncol(mod)])%*%Dp%*%mod[,-ncol(mod)]

    if(checkInverse(infz)) {
      invinf_z <-qr.solve(infz)
    } else {
      return(NA)
    }
    if(checkInverse(infln)) {
      invinf_ln<-as.numeric(s)*qr.solve(infln)
    } else {
      return(NA)
    }
    invInfFull = as.matrix( bdiag(invinf_z,invinf_ln, (2*s^2/sum(p))) )

    logRatioBeta0<- (mean(p1*(1-p1)*mu0)/mean(p1*mu0)) - (mean(p0*(1-p0)*mu0)/mean(p0*mu0))
    logRatioBeta1<-mean(p1*(1-p1)*mu0)/mean(p1*mu0)
    logRatioBeta2<- (mean(mod[,3]*p1*(1-p1)*mu0)/mean(p1*mu0)) - (mean(mod[,3]*p0*(1-p0)*mu0)/mean(p0*mu0))
    # logRatioB2<- (mean(mod[,3]*t(p1)*exp(mod0%*%t(b)))/mean(t(p1)*exp(mod0%*%t(b))))-
    #    (mean(mod[,3]*t(p0)*exp(mod0%*%t(b)))/mean(t(p0)*exp(mod0%*%t(b))))
    # logRatioFull = t(c(logRatioBeta0,logRatioBeta1,logRatioBeta2,0,1,logRatioB2,0))
    logRatioFull = t(c(logRatioBeta0,logRatioBeta1,logRatioBeta2,0,1,0))
    logRatioVar = logRatioFull%*%invInfFull%*%t(logRatioFull)
    logRatioVar
  })
  if(!is.null(exclude)){
    if(length(features)>0){
      ve[features] = fullvar
    }
  } else {
    ve = fullvar
  }
  sqrt(ve)
}

#' Table of top-ranked features from fitZig or fitFeatureModel
#'
#' Extract a table of the top-ranked features from a linear model fit. This
#' function will be updated soon to provide better flexibility similar to
#' limma's topTable.
#'
#'
#' @param obj Output of fitFeatureModel or fitZig.
#' @param by Column number or column name specifying which coefficient or
#' contrast of the linear model is of interest.
#' @param coef Column number(s) or column name(s) specifying which coefficient
#' or contrast of the linear model to display.
#' @param number The number of bacterial features to pick out.
#' @param taxa Taxa list.
#' @param uniqueNames Number the various taxa.
#' @param adjustMethod Method to adjust p-values by. Default is "FDR". Options
#' include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' "none". See p.adjust for more details. Additionally, options using
#' independent hypothesis weighting (IHW) are available. See MRihw for more
#' details.
#' @param alpha Value for p-value significance threshold when running IHW.
#' The default is set to 0.1
#' @param group One of five choices, 0,1,2,3,4. 0: the sort is ordered by a
#' decreasing absolute value coefficient fit. 1: the sort is ordered by the raw
#' coefficient fit in decreasing order. 2: the sort is ordered by the raw
#' coefficient fit in increasing order. 3: the sort is ordered by the p-value
#' of the coefficient fit in increasing order. 4: no sorting.
#' @param eff Filter features to have at least a "eff" quantile or number of effective samples.
#' @param numberEff Boolean, whether eff should represent quantile (default/FALSE) or number.
#' @param counts Filter features to have at least 'counts' counts.
#' @param file Name of output file, including location, to save the table.
#' @return Table of the top-ranked features determined by the linear fit's
#' coefficient.
#' @examples
#'
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim=filterData(lungTrim,present=30)
#' lungTrim=cumNorm(lungTrim,p=0.5)
#' smokingStatus = pData(lungTrim)$SmokingStatus
#' mod = model.matrix(~smokingStatus)
#' fit = fitZig(obj = lungTrim,mod=mod)
#' head(MRcoefs(fit))
#' ####
#' fit = fitFeatureModel(obj = lungTrim,mod=mod)
#' head(MRcoefs(fit))
#'
mvMRcoefs<-function(obj,by=2,coef=NULL,number=10,taxa=obj@taxa,
                  uniqueNames=FALSE,adjustMethod="fdr",alpha=0.1,
                  group=0,eff=0,numberEff=FALSE,counts=0,file=NULL){

  if(length(grep("fitFeatureModel",obj@call))){
    groups = factor(obj@design[,by])
    by = "logFC"; coef = 1:2;
    tb = data.frame(logFC=obj@fitZeroLogNormal$logFC,se=obj@fitZeroLogNormal$se)
    p  = obj@pvalues
  } else {
    tb = obj@fit$coefficients
    if(is.null(coef)){
      coef = 1:ncol(tb)
    }
    p=obj@eb$p.value[,by]
    groups = factor(obj@fit$design[,by])
    if(eff>0){
      effectiveSamples = calculateEffectiveSamples(obj)
      if(numberEff == FALSE){
        valid = which(effectiveSamples>=quantile(effectiveSamples,p=eff,na.rm=TRUE))
      } else {
        valid = which(effectiveSamples>=eff)
      }
    }
  }

  tx = as.character(taxa)
  if(uniqueNames==TRUE){
    for (nm in unique(tx)) {
      ii=which(tx==nm)
      tx[ii]=paste(tx[ii],seq_along(ii),sep=":")
    }
  }

  # adding 'ihw' as pvalue adjustment method
  if (adjustMethod == "ihw-ubiquity" | adjustMethod == "ihw-abundance") {
    # use IHW to adjust pvalues
    padj = MRihw(obj, p, adjustMethod, alpha)
  } else {
    # use classic pvalue adjusment method
    padj = p.adjust(p, method = adjustMethod)
  }

  if(group==0){
    srt = order(abs(tb[,by]),decreasing=TRUE)
  } else if(group==1){
    srt = order((tb[,by]),decreasing=TRUE)
  } else if(group==2){
    srt = order((tb[,by]),decreasing=FALSE)
  } else if(group==3){
    srt = order(p,decreasing=FALSE)
  } else {
    srt = 1:length(padj);
  }

  valid = 1:length(padj);
  if(counts>0){
    np=rowSums(obj@counts);
    valid = intersect(valid,which(np>=counts));
  }
  srt = srt[which(srt%in%valid)][1:min(number,nrow(tb))];

  mat = cbind(tb[,coef],p)
  mat = cbind(mat,padj)
  rownames(mat) = tx;
  mat = mat[srt,]

  nm = c(colnames(tb)[coef],"pvalues","adjPvalues")
  colnames(mat) = nm

  if(!is.null(file)){
    nm = c("Taxa",nm)
    mat2 = cbind(rownames(mat),mat)
    mat2 = rbind(nm,mat2)
    write(t(mat2),ncolumns=ncol(mat2),file=file,sep="\t")
  }
  return(as.data.frame(mat))
}

#' Table of top microbial marker gene from linear model fit including sequence
#' information
#'
#' Extract a table of the top-ranked features from a linear model fit. This
#' function will be updated soon to provide better flexibility similar to
#' limma's topTable. This function differs from \code{link{MRcoefs}} in that it
#' provides other information about the presence or absence of features to help
#' ensure significant features called are moderately present.
#'
#'
#' @param obj Output of fitFeatureModel or fitZig.
#' @param by Column number or column name specifying which coefficient or
#' contrast of the linear model is of interest.
#' @param coef Column number(s) or column name(s) specifying which coefficient
#' or contrast of the linear model to display.
#' @param number The number of bacterial features to pick out.
#' @param taxa Taxa list.
#' @param uniqueNames Number the various taxa.
#' @param adjustMethod Method to adjust p-values by. Default is "FDR". Options
#' include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' "none". See p.adjust for more details.
#' @param group One of five choices, 0,1,2,3,4. 0: the sort is ordered by a
#' decreasing absolute value coefficient fit. 1: the sort is ordered by the raw
#' coefficient fit in decreasing order. 2: the sort is ordered by the raw
#' coefficient fit in increasing order. 3: the sort is ordered by the p-value
#' of the coefficient fit in increasing order. 4: no sorting.
#' @param eff Filter features to have at least a "eff" quantile or number of effective samples.
#' @param numberEff Boolean, whether eff should represent quantile (default/FALSE) or number.
#' @param ncounts Filter features to have at least 'counts' of counts.
#' @param file Name of file, including location, to save the table.
#' @return Table of the top-ranked features determined by the linear fit's
#' coefficient.
#' @examples
#'
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim=filterData(lungTrim,present=30)
#' lungTrim=cumNorm(lungTrim,p=0.5)
#' smokingStatus = pData(lungTrim)$SmokingStatus
#' mod = model.matrix(~smokingStatus)
#' fit = fitZig(obj = lungTrim,mod=mod)
#' head(MRtable(fit))
#' ####
#' fit = fitFeatureModel(obj = lungTrim,mod=mod)
#' head(MRtable(fit))
#'
mvMRtable<-function(obj,by=2,coef=NULL,number=10,taxa=obj@taxa,
                  uniqueNames=FALSE,adjustMethod="fdr",group=0,eff=0,numberEff=FALSE,ncounts=0,file=NULL){

  if(length(grep("fitFeatureModel",obj@call))){
    groups = factor(obj@design[,by])
    by = "logFC"; coef = 1:2;
    tb = data.frame(logFC=obj@fitZeroLogNormal$logFC,se=obj@fitZeroLogNormal$se)
    p  = obj@pvalues
  } else {
    tb = obj@fit$coefficients
    if(is.null(coef)){
      coef = 1:ncol(tb)
    }
    p=obj@eb$p.value[,by]
    groups = factor(obj@fit$design[,by])
    if(eff>0){
      effectiveSamples = calculateEffectiveSamples(obj)
      if(numberEff == FALSE){
        valid = which(effectiveSamples>=quantile(effectiveSamples,p=eff,na.rm=TRUE))
      } else {
        valid = which(effectiveSamples>=eff)
      }
    }
  }

  tx = as.character(taxa)
  if(uniqueNames==TRUE){
    for (nm in unique(tx)) {
      ii=which(tx==nm)
      tx[ii]=paste(tx[ii],seq_along(ii),sep=":")
    }
  }
  padj = p.adjust(p,method=adjustMethod)
  cnts = obj@counts
  posIndices = cnts>0

  np0 = rowSums(posIndices[,groups==0])
  np1 = rowSums(posIndices[,groups==1])

  nc0 = rowSums(cnts[,groups==0])
  nc1 = rowSums(cnts[,groups==1])

  if(group==0){
    srt = order(abs(tb[,by]),decreasing=TRUE)
  } else if(group==1){
    srt = order((tb[,by]),decreasing=TRUE)
  } else if(group==2){
    srt = order((tb[,by]),decreasing=FALSE)
  } else if(group==3){
    srt = order(p,decreasing=FALSE)
  } else {
    srt = 1:length(padj)
  }

  valid = 1:length(padj)
  if(ncounts>0){
    np=rowSums(cbind(np0,np1))
    valid = intersect(valid,which(np>=ncounts))
  }
  srt = srt[which(srt%in%valid)][1:min(number,nrow(tb))]

  mat = cbind(np0,np1)
  mat = cbind(mat,nc0)
  mat = cbind(mat,nc1)
  mat = cbind(mat,tb[,coef])
  mat = cbind(mat,p)
  mat = cbind(mat,padj)
  rownames(mat) = tx
  mat = mat[srt,]

  nm = c("+samples in group 0","+samples in group 1","counts in group 0",
         "counts in group 1",colnames(tb)[coef],"pvalues","adjPvalues")
  colnames(mat) = nm

  if(!is.null(file)){
    nm = c("Taxa",nm)
    mat2 = cbind(rownames(mat),mat)
    mat2 = rbind(nm,mat2)
    write(t(mat2),ncolumns=ncol(mat2),file=file,sep="\t")
  }
  return(as.data.frame(mat))
}

#' Table of top microbial marker gene from linear model fit including sequence
#' information
#'
#' Extract a table of the top-ranked features from a linear model fit. This
#' function will be updated soon to provide better flexibility similar to
#' limma's topTable. This function differs from \code{link{MRcoefs}} in that it
#' provides other information about the presence or absence of features to help
#' ensure significant features called are moderately present.
#'
#'
#' @param obj  Output of fitFeatureModel or fitZig.
#' @param by Column number or column name specifying which coefficient or
#' contrast of the linear model is of interest.
#' @param coef Column number(s) or column name(s) specifying which coefficient
#' or contrast of the linear model to display.
#' @param number The number of bacterial features to pick out.
#' @param taxa Taxa list.
#' @param uniqueNames Number the various taxa.
#' @param adjustMethod Method to adjust p-values by. Default is "FDR". Options
#' include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr",
#' "none". See p.adjust for more details.
#' @param group One of five choices: 0,1,2,3,4. 0: the sort is ordered by a
#' decreasing absolute value coefficient fit. 1: the sort is ordered by the raw
#' coefficient fit in decreasing order. 2: the sort is ordered by the raw
#' coefficient fit in increasing order. 3: the sort is ordered by the p-value
#' of the coefficient fit in increasing order. 4: no sorting.
#' @param eff Filter features to have at least a "eff" quantile or number of effective samples.
#' @param numberEff Boolean, whether eff should represent quantile (default/FALSE) or number.
#' @param ncounts Filter features to those with at least 'counts' counts.
#' @param file Name of output file, including location, to save the table.
#' @return Table of the top-ranked features determined by the linear fit's
#' coefficient.
#' @examples
#'
#' data(lungData)
#' k = grep("Extraction.Control",pData(lungData)$SampleType)
#' lungTrim = lungData[,-k]
#' lungTrim=filterData(lungTrim,present=30)
#' lungTrim=cumNorm(lungTrim,p=0.5)
#' smokingStatus = pData(lungTrim)$SmokingStatus
#' mod = model.matrix(~smokingStatus)
#' fit = fitZig(obj = lungTrim,mod=mod)
#' head(MRfulltable(fit))
#' ####
#' fit = fitFeatureModel(obj = lungTrim,mod=mod)
#' head(MRfulltable(fit))
#'
mvMRfulltable<-function(obj,by=2,coef=NULL,number=10,taxa=obj@taxa,
                      uniqueNames=FALSE,adjustMethod="fdr",group=0,eff=0,numberEff=FALSE,ncounts=0,file=NULL){

  if(length(grep("fitFeatureModel",obj@call))){
    groups = factor(obj@design[,by])
    by = "logFC"; coef = 1:2;
    tb = data.frame(logFC=obj@fitZeroLogNormal$logFC,se=obj@fitZeroLogNormal$se)
    p  = obj@pvalues
  } else {
    tb = obj@fit$coefficients
    if(is.null(coef)){
      coef = 1:ncol(tb)
    }
    p=obj@eb$p.value[,by]
    groups = factor(obj@fit$design[,by])
    if(eff>0){
      effectiveSamples = calculateEffectiveSamples(obj)
      if(numberEff == FALSE){
        valid = which(effectiveSamples>=quantile(effectiveSamples,p=eff,na.rm=TRUE))
      } else {
        valid = which(effectiveSamples>=eff)
      }
    }
  }

  tx = as.character(taxa)
  if(uniqueNames==TRUE){
    for (nm in unique(tx)) {
      ii=which(tx==nm)
      tx[ii]=paste(tx[ii],seq_along(ii),sep=":")
    }
  }
  padj = p.adjust(p,method=adjustMethod)
  cnts = obj@counts
  yy = cnts>0

  pa = matrix(unlist(fitPA(obj@counts,groups)),ncol=5)

  np0 = rowSums(yy[,groups==0])
  np1 = rowSums(yy[,groups==1])

  nc0 = rowSums(cnts[,groups==0])
  nc1 = rowSums(cnts[,groups==1])

  if(group==0){
    srt = order(abs(tb[,by]),decreasing=TRUE)
  } else if(group==1){
    srt = order((tb[,by]),decreasing=TRUE)
  } else if(group==2){
    srt = order((tb[,by]),decreasing=FALSE)
  } else if(group==3){
    srt = order(p,decreasing=FALSE)
  } else {
    srt = 1:length(padj)
  }

  valid = 1:length(padj)
  if(ncounts>0){
    np=rowSums(cbind(np0,np1))
    valid = intersect(valid,which(np>=ncounts))
  }
  srt = srt[which(srt%in%valid)][1:min(number,nrow(tb))]

  mat = cbind(np0,np1)
  mat = cbind(mat,nc0)
  mat = cbind(mat,nc1)
  mat = cbind(mat,pa)
  mat = cbind(mat,tb[,coef])
  mat = cbind(mat,p)
  mat = cbind(mat,padj)
  rownames(mat) = tx
  mat = mat[srt,]

  nm = c("+samples in group 0","+samples in group 1","counts in group 0",
         "counts in group 1",c("oddsRatio","lower","upper","fisherP","fisherAdjP"),
         colnames(tb)[coef],"pvalues","adjPvalues")
  colnames(mat) = nm

  if(!is.null(file)){
    nm = c("Taxa",nm)
    mat2 = cbind(rownames(mat),mat)
    mat2 = rbind(nm,mat2)
    write(t(mat2),ncolumns=ncol(mat2),file=file,sep="\t")
  }
  return(as.data.frame(mat))
}

#' CClasso Network Analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param r_cutoff Coefficient cutoff. Defaults to 0 (none)
#' @param features1 (Optional) Features to consider for set 1
#' @param features2 (Optional) Features to consider for set 2
#'
#' @return Correlation matrices
#' @export
#'
mvCClasso <- function(dataset=NULL,r_cutoff=0,features1=NULL,features2=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- setFVar(dataset)

  mdcolnum <- ncol(dataset$metadata)

  melted <- mvmelt(dataset)
  melted$Other <- NULL

  cor_list <- list()
  for(grp in factor$subset) {
    abd <- melted[melted[[factor$name]]==grp,(mdcolnum+1):ncol(melted)]
    cor.tmp <- cclasso(abd)$cor.w
    cor.tmp[abs(cor.tmp)<r_cutoff] <- 0

    features1 <- features1[features1 %in% colnames(cor.tmp)]
    features2 <- features2[features2 %in% colnames(cor.tmp)]
    if(!is.null(features1) & !is.null(features2)) cor.tmp <- cor.tmp[features1,features2]

    cor_list[[grp]] <- cor.tmp
  }

  na_fts1 <- c()
  na_fts2 <- c()
  for(grp in names(cor_list)) {
    na_fts1 <- union(na_fts1,rownames(cor_list[[grp]])[is.na(rowSums(cor_list[[grp]]))])
    na_fts2 <- union(na_fts2,colnames(cor_list[[grp]])[is.na(colSums(cor_list[[grp]]))])
  }

  hide_fts1 <- rownames(cor_list[[1]])
  hide_fts2 <- colnames(cor_list[[1]])
  for(grp in names(cor_list)) {
    zero_fts1 <- rownames(cor_list[[grp]])[rowSums(cor_list[[grp]],na.rm = T)==0]
    hide_fts1 <- intersect(hide_fts1,c(zero_fts1,na_fts1))

    zero_fts2 <- colnames(cor_list[[grp]])[colSums(cor_list[[grp]],na.rm = T)==0]
    hide_fts2 <- intersect(hide_fts2,c(zero_fts2,na_fts2))
  }

  hide_fts1 <- union(hide_fts1,na_fts1)
  hide_fts2 <- union(hide_fts2,na_fts2)

  keep_fts1 <- features1[!(features1 %in% hide_fts1)]
  keep_fts2 <- features2[!(features2 %in% hide_fts2)]

  for(grp in names(cor_list)) cor_list[[grp]] <- cor_list[[grp]][keep_fts1,keep_fts2]

  return(cor_list)
  # p <- list()
  # for(grp in names(cor_list)) {
  #   p[[grp]] <- ggcorrplot(cor_list[[grp]])+
  #     labs(title = grp)+
  #     theme(title = element_text(hjust = 0.5))
  #   plot(p[[grp]])
  # }
  #
  # return(p)
}

cclasso <- function(x, counts = F, pseudo = 0.5, sig = NULL,
                    lams = 10^(seq(0, -8, by = -0.01)),
                    K = 3, kmax = 5000) {
  # data dimension
  p <- ncol(x);
  n <- nrow(x);

  # # Counts or Fractions?
  # if(counts) {
  #   x <- x + pseudo;
  #   x <- x / rowSums(x);
  # }
  # # log transformation
  # xlog <- log(x);
  xlog <- log((x+sqrt(x^2+4))/2)

  # use all data
  vx <- stats::var(xlog);
  # initial value
  res <- list();
  if(is.null(sig)) {
    res$sig <- diag(rep(1, p));
  }
  else {
    res$sig <- sig;
  }
  #-------------------------------------------------------------------------------
  # weight diagonal for loss
  rmean.vx <- rowMeans(vx);
  wd <- 1 / diag(vx - rmean.vx - rep(rmean.vx, each = p) + mean(rmean.vx));
  wd2 <- sqrt(wd);
  #-------------------------------------------------------------------------------
  # preparation for update sigma in augmented lagrange method
  rho <- 1; # needed
  u.f <- eigen(diag(rep(1, p)) - 1 / p)$vectors; # needed
  wd.u <- (t(u.f) %*% (wd * u.f))[-p, -p];
  diag(wd.u) <- diag(wd.u) + rho;
  wd.u.eig <- eigen(wd.u);
  d0.wd <- 2 / outer(wd.u.eig$values, wd.u.eig$values, "+"); # needed
  u0.wd <- wd.u.eig$vectors; # needed
  #-------------------------------------------------------------------------------
  n_lam <- length(lams);
  tol.zero <- 1e-8;
  #-------------------------------------------------------------------------------
  # cross validation
  if(n_lam == 1) {
    lamA <- lams[1];
  }
  else {
    tol.loss <- 1e-6;
    loss.old <- Inf;
    k.loss <- n_lam;
    n.b <- floor(n / K);
    # loss <- rep(0, n_lam);
    cat('Testing lambda values:')
    for(i in 1:n_lam) {
      cat('',i)
      loss.cur <- 0;
      for(k in 1:K) {
        # testing data and training data
        itest <- (n.b * (k-1) + 1):(n.b * k);
        vxk <- stats::var(xlog[itest, ]);
        #
        vx2k <- stats::var(xlog[-itest, ]);
        # for training
        res <- cclasso.sub(vx = vx2k, wd = wd, lam = lams[i],
                           u.f = u.f, u0.wd = u0.wd, d0.wd = d0.wd,
                           sig = res$sig, rho = rho, kmax = kmax);
        # loss cumulation
        res$sig[abs(res$sig) <= tol.zero] <- 0;
        dsig <- res$sig - vxk;
        rmean.dsig <- rowMeans(dsig);
        half.loss <- dsig * rep(wd2, each = p) -
          outer(rmean.dsig, wd2, "*") +
          rep(  (mean(rmean.dsig) - rmean.dsig) * wd2, each = p);
        # loss[i] <- loss[i] + base::norm(half.loss, "F")^2;
        loss.cur <- loss.cur + base::norm(half.loss, "F")^2;
      }
      if(loss.cur - loss.old >= tol.loss * max(loss.cur, loss.old, 1)) {
        k.loss <- i - 1;
        break;
      }
      else {
        loss.old <- loss.cur;
      }
    }
    cat('\n')
    # select lambda
    # k.loss <- which.min(loss);
    lamA <- lams[k.loss];
    if (k.loss == 1 || k.loss == n_lam) {
      cat("Warning:", "Tuning (", lamA ,") on boundary!\n");
    }
    cat('\nSelected ',lamA,' for lambda')
  }
  #-------------------------------------------------------------------------------
  cat('\nRunning CCLasso\n')
  res <- cclasso.sub(vx = vx, wd = wd, lam = lamA,
                     u.f = u.f, u0.wd = u0.wd, d0.wd = d0.wd,
                     sig = res$sig, rho = rho, kmax = kmax);
  #-------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------
  #
  res$sig[abs(res$sig) <= tol.zero] <- 0;
  if(min(eigen(res$sig)$values) <= tol.zero) {
    sig.sparse <- abs(res$sig) > tol.zero;
    diag(res$sig) <- diag(res$sig) * sign(diag(res$sig));
    res$sig <- as.matrix(nearPD(res$sig)$mat) * sig.sparse;
  }
  # get correlation matrix from covariance matrix
  Is <- sqrt(1 / diag(res$sig));
  cor.w <- Is * res$sig * rep(Is, each = p);
  # remove too small correlation values
  cor.w[abs(cor.w) <= 1e-6] <- 0;
  #
  return(list(cov.w = res$sig, cor.w = cor.w, lam = lamA));
}
# cclasso for only one lambda
cclasso.sub <- function(vx, wd, lam, u.f, u0.wd, d0.wd, sig = NULL,
                        rho = 1, kmax = 5000, x.tol = 1e-6) {
  p <- ncol(vx);
  # initial value
  lam.rho <- lam / rho;
  if(is.null(sig)) {
    sig <- diag(rep(1, p));
  }
  sig2 <- sig;
  LAM <- matrix(0, p, p);
  # loop start
  k <- 0;
  err <- 1;
  while(err > x.tol && k < kmax) {
    # update sigma
    x.sig <- t(u.f) %*% ((sig2  -  vx) - LAM  / rho) %*% u.f;
    x.sig[-p,-p] <- u0.wd %*% ((t(u0.wd) %*% x.sig[-p, -p] %*% u0.wd) *
                                 d0.wd * rho) %*% t(u0.wd);
    sig.new <- vx + u.f %*% x.sig %*% t(u.f);
    # update sigma2
    A <- LAM / rho + sig.new;
    sig2.new <- (A > lam.rho) * (A - lam.rho) +
      (A < -lam.rho) * (A + lam.rho);
    diag(sig2.new) <- diag(A);
    # update Lambda
    LAM <- LAM + rho * (sig.new - sig2.new);
    # calculate error
    err <- max( base::norm(sig.new - sig, "F") / max(1, base::norm(sig)),
                base::norm(sig2.new - sig2, "F") / max(1, base::norm(sig2)));
    # update current value
    sig <- sig.new;
    sig2 <- sig2.new;
    k <- k + 1;
  }
  #
  if(k >= kmax) {
    cat("Warning:", "Maximum ", kmax,
        "iteration while relative error is", err, "\n");
  }
  #
  return(list(sig = sig, k = k));
}

#' @title fastANCOM, a fast method for analysis of composition of microbiomes
#'
#'
#' @description fastANCOM fits the mariginal model and then computing the effect for the joint model.
#' The framework of log-linear regression is used for the log transformed data.
#' The hypothesis is build with the variable of interest for inference at the ecosystem
#' level with the framework of ANCOM.
#'
#'
#' @param Y the abundance table,  count/relative or any smoothed data matrix are supported
#' @param x the variable of interest
#' @param Z the other covariates to be adjusted
#' @param rand the random variable of the mixed effect model
#' @param pseudo the pseudo number for zero smoothing, default is 0.05
#' @param outlier a boolean value, if TRUE the cook's distance is used for outlier detection, the default is FALSE
#' @param sig the level of significance to count the reject number, default is 0.05
#' @param detect.rate a numerical fraction between 0 and 1. Miorobes with proportion of the reject number greater than it will be identified as associated with the variable of interest, i.e., x, default is 0.7
#' @param ref.rate a numerical fraction between 0 and 0.5, the proportion of total microbes with the least reject numbers and these microbes are set as the reference to estimate the effect size in ecosystem, the default is 0.05
#'
#' @return a dataframe with 5 columns,
#' \itemize{
#' \item log2FC, the adjusted effect size of log-linear model
#' \item log2FC.SD, the standard errors (SEs) of log2FC
#' \item log2FC.pval, the p-value of Wald test
#' \item Reject.number, the reject number of each microbe
#' \item REJECT, the final result whether each microbe is associated with the variable of interest or not
#' }
#'
#' @examples
#' data <- matrix(rpois(100*60, 10), 60)
#' b <- rnorm(10, 4)
#' data[31:60,1:10] <- t(t(data[31:60,1:10])*2^b)
#' group <- rep(c(0, 1), each=30)
#' fit <- fastANCOM(Y=data, x=group)
#' head(fit)
#' @export
#'


fastANCOM <- function(Y, x, Z=NULL, rand=NULL, pseudo=0.5, outlier=F,
                        sig=0.05, detect.rate=0.7, ref.rate=0.05){
  n_otu <- ncol(Y)
  n_sample <- nrow(Y)
  logdata <- log2(Y + pseudo)
  logdata[is.infinite(logdata)] <- NA
  if(is.null(rand)) { # repeated the sample ID match OTUdata
    if(is.null(Z)){
      tmp <- lm(logdata ~ x)
      if(outlier) {
        out <- cooks.distance(tmp)
        logdata[out>4*mean(out, na.rm=T)] <- NA
        tmp <- lm(logdata ~ x)
      }
    }else{
      tmp <- lm(logdata ~ x + Z)
      if(outlier) {
        out <- cooks.distance(tmp)
        logdata[out>4*mean(out, na.rm=T)] <- NA
        tmp <- lm(logdata ~ x + Z)
      }
    }
    beta <- tmp$coefficients[2,]
    xtx <- chol2inv(tmp$qr$qr[1:tmp$rank, 1:tmp$rank, drop = FALSE])
    covbeta <- t(tmp$residuals) %*% tmp$residuals * xtx[2,2] /tmp$df.residual
    b12 <- beta - matrix(beta, nrow=n_otu, ncol=n_otu, byrow = T)
    s12 <- matrix(diag(covbeta), nrow=n_otu, ncol=n_otu, byrow = T)
    se12<- sqrt(diag(covbeta) + s12 - 2*covbeta)
    tstat <- -abs(b12/se12)
    ta <- qt(sig/2, df = n_sample-2)
    W1 <- rowSums(tstat < ta, na.rm = T)

  }else{
    logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
    if(is.null(Z)) Z <- 1
    for(ii in 1:(n_otu-1)){
      data.pair <- logdata[,ii:n_otu]
      lrii <- data.pair[,1] - data.pair[,-1]
      logratio.mat[ii,(ii+1):n_otu] <- apply(as.matrix(lrii), 2, function(y) {
        ti <- lme(y ~ x + Z, random=rand); anova(ti)['x','p-value']})
    }
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
    logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
    W1 <- rowSums(logratio.mat<sig, na.rm = T)
  }

  decisionMaking <- function(W, detect=0.7, tau=0.02, theta=0.26){
    W_stat  <- W
    num_OTU <- length(W_stat)
    detected <- detect0.7 <- rep(F, num_OTU)
    if( num_OTU < 10 ){
      detected[which(W_stat > num_OTU )] <- T
    } else{
      ## Detected using a stepwise mode detection
      if( max(W_stat)/num_OTU >= theta ){
        c.start <- max(W_stat)/num_OTU
        cutoff  <- c.start-c(0.05,0.10,0.15,0.20,0.25)
        # cutoff <- cutoff[cutoff>0]

        prop_cut <- rep(0,length(cutoff))
        for(cut in 1:length(cutoff)){
          prop_cut[cut] <- sum(W_stat>=num_OTU*cutoff[cut])/num_OTU
        }

        del <- rep(0,length(cutoff)-1)
        for( ii in 1:(length(cutoff)-1) ){
          del[ii] <- abs(prop_cut[ii]-prop_cut[ii+1])
        }

        if(       del[1]< tau & del[2]<tau & del[3]<tau ){ nu=cutoff[1]
        }else if( del[1]>=tau & del[2]<tau & del[3]<tau ){ nu=cutoff[2]
        }else if( del[2]>=tau & del[3]<tau & del[4]<tau ){ nu=cutoff[3]
        }else{ nu=cutoff[4] }

        # up_point <- min(W_stat[  W_stat >= nu*num_OTU  ])
        # detected[W_stat>up_point] <- T
        detected[W_stat>=nu*num_OTU] <- T

      } else{
        detected   <- F
      }
    }
    detect0.7[W_stat>=0.7*num_OTU] <- T
    results <- data.frame(W=W_stat, detected=detected, detect0.7=detect0.7)
    return(results)
  }
  tmpres <- decisionMaking(W1)
  {
    rfn <- ref.rate*n_otu
    if(rfn<1) rfn <- 1
    if(n_otu<10) {
      refset <- order(W1)[1]
    }else{
      # refset <- which(W1<rfn & tmpres$detected==F)
      # rfn0   <- ifelse(length(refset) > rfn, rfn, length(refset)*0.9)
      refset <- order(W1)[1:rfn]
    }
    rawp <- 2*pt(-abs(beta/diag(covbeta)), df=n_sample-2)
    adj.beta <- beta - mean(beta[refset])
    sers12 <- covbeta[refset, refset]
    se0012 <- covbeta[, refset]
    adj.se12 <- sqrt(diag(covbeta) - 2/n_otu*rowSums(se0012) + sum(sers12)/n_otu^2)
    jpvs <- 2*pt(-abs(adj.beta/adj.se12), df=n_sample-2)
  }
  res <- data.frame(log2FC=adj.beta, log2FC.SD=adj.se12, log2FC.pval=jpvs,
                    Reject.number=W1, REJECT=tmpres$detected)
  res
}

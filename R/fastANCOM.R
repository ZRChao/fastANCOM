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
#' @param zero_cut a numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than zero_cut will be excluded in the analysis. Default is 0.90
#' @param lib_cut a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, i.e. do not filter any sample.
#' @param struc_zero whether to detect structural zeros. Default is FALSE
#' @param pseudo the pseudo number for zero smoothing, default is 0.05
#' @param global whether to perform global test for more than two groups. Default is FALSE
#' @param outlier a boolean value, if TRUE the cook's distance is used for outlier detection, the default is FALSE
#' @param sig the level of significance to count the reject number, default is 0.05
#' @param detect.rate a numerical fraction between 0 and 1. Miorobes with proportion of the reject number greater than it will be identified as associated with the variable of interest, i.e., x, default is 0.7
#' @param ref.rate a numerical fraction between 0 and 0.5, the proportion of total microbes with the least reject numbers and these microbes are set as the reference to estimate the effect size in ecosystem, the default is 0.05
#' @param effect a boolean value, if TRUE then estimate the effect size based on the reference set, the default is T
#'
#' @return a list with two components:
#' \itemize{
#' \item global, if the variable of interest x is continuous or global==FALSE, it returns NULL, otherwise a list with three components
#' \itemize{
#' \item marginal, a list contains the estimated coefficient and the covariance structure of the marginal model
#' \item joint, a list contains 4 elements for the joint model: beta, the coefficients from the joint model; se, the corresponding standard error; statistic, the Wald statistic; and W, the reject number of each component.
#' \item final, a data.frame contains the 5 columns for the final effect size estimation,
#' \itemize{
#' \item log2FC, the adjusted effect size of log-linear model
#' \item log2FC.SD, the standard errors (SEs) of log2FC
#' \item log2FC.pval, the p-value of Wald test
#' \item Reject.number, the reject number of each microbe
#' \item REJECT, the final result whether each microbe is associated with the variable of interest or not
#' }
#' }
#'
#' \item results, a list contains the results for the pairwise group comparison, each element will be same structure as above
#' }
#'
#' @examples
#' data <- matrix(rpois(100*60, 10), 60)
#' b <- rnorm(10, 4)
#' data[31:60,1:10] <- t(t(data[31:60,1:10])*2^b)
#' group <- rep(c(0, 1), each=30)
#' fit <- fastANCOM(Y=data, x=group)
#' summary(fit)
#' head(fit$results$final)
#' @export
#'


fastANCOM <- function(Y, x, Z=NULL, rand=NULL,
                      zero_cut = 0.9, lib_cut = 0, struc_zero = FALSE,
                      pseudo=0.5, global=F, outlier=F,
                      sig=0.05, detect.rate=0.7, ref.rate=0.05, effect=T) {

  feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL,
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb=F){
    feature_table = data.frame(feature_table, check.names = FALSE)
    meta_data = data.frame(meta_data, check.names = FALSE)
    # Drop unused levels
    meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
    # Match sample IDs between metadata and feature table
    sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
    feature_table = feature_table[, sample_ID]
    meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]

    # 1. Identify outliers within each taxon
    if (!is.null(group_var)) {
      group = meta_data[, group_var]
      z = feature_table + 1 # Add pseudo-count (1)
      f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
      f_fit = lm(f ~ group)
      e = rep(0, length(f)); e[!is.na(group)] = residuals(f_fit)
      y = t(t(z) - e)

      outlier_check = function(x){
        # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
        mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
        sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
        pi = 0.75
        n = length(x)
        epsilon = 100
        tol = 1e-5
        score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
        while (epsilon > tol) {
          grp1_ind = (score >= 1)
          mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
          sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
          sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
          pi_new = sum(grp1_ind)/n

          para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
          if(any(is.na(para))) break

          score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
            ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

          epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 +
                           (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
          mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
        }

        if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
          if(pi < out_cut){
            out_ind = grp1_ind
          }else if(pi > 1 - out_cut){
            out_ind = (!grp1_ind)
          }else{
            out_ind = rep(FALSE, n)
          }
        }else{
          out_ind = rep(FALSE, n)
        }
        return(out_ind)
      }
      out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
      out_ind[, !is.na(group)] = t(apply(y, 1, function(i)
        unlist(tapply(i, group, function(j) outlier_check(j)))))

      feature_table[out_ind] = NA
    }

    # 2. Discard taxa with zeros  >=  zero_cut
    zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
    taxa_del = which(zero_prop >= zero_cut)
    if(length(taxa_del) > 0){
      feature_table = feature_table[- taxa_del, ]
    }

    # 3. Discard samples with library size < lib_cut
    lib_size = colSums(feature_table, na.rm = T)
    if(any(lib_size < lib_cut)){
      subj_del = which(lib_size < lib_cut)
      feature_table = feature_table[, - subj_del]
      meta_data = meta_data[- subj_del, ]
    }

    # 4. Identify taxa with structure zeros
    if (!is.null(group_var)) {
      group = factor(meta_data[, group_var])
      present_table = as.matrix(feature_table)
      present_table[is.na(present_table)] = 0
      present_table[present_table != 0] = 1

      p_hat = t(apply(present_table, 1, function(x)
        unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
      samp_size = t(apply(feature_table, 1, function(x)
        unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
      p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

      struc_zero = (p_hat == 0) * 1
      # Whether we need to classify a taxon into structural zero by its negative lower bound?
      if(neg_lb) struc_zero[p_hat_lo <= 0] = 1

      # Entries considered to be structural zeros are set to be 0s
      struc_ind = struc_zero[, group]
      feature_table = feature_table * (1 - struc_ind)

      colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
    }else{
      struc_zero = NULL
    }

    # 5. Return results
    res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
    return(res)
  }

  .fANCOM <- function(Y, x, Z=NULL, rand=NULL, pseudo=0.5, outlier=F,
                      sig=0.05, detect.rate=0.7, ref.rate=0.05, effect=T){
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

      res.m <- list(beta=beta, cov=covbeta)
      res.j <- list(beta=b12, se=se12, statistic=tstat, W=W1)
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
    res <- tmpres <- decisionMaking(W1)
    if(effect){
      rfn <- ref.rate*n_otu
      if(rfn<1) rfn <- 1
      if(n_otu<10) {
        rfn <- 1
        refset <- order(W1)[1]
      }else{
        # refset <- which(W1<rfn & tmpres$detected==F)
        # rfn0   <- ifelse(length(refset) > rfn, rfn, length(refset)*0.9)
        refset <- order(W1)[1:rfn]
      }
      if(rfn<5) warning(paste0('Only take ', rfn, ' components as reference!'))

      rawp <- 2*pt(-abs(beta/diag(covbeta)), df=n_sample-2)
      adj.beta <- beta - mean(beta[refset])
      sers12 <- covbeta[refset, refset]
      se0012 <- as.matrix(covbeta[, refset])
      adj.se12 <- sqrt(diag(covbeta) - 2/rfn*rowSums(se0012) + sum(sers12)/rfn^2)
      jpvs <- 2*pt(-abs(adj.beta/adj.se12), df=n_sample-2)
      jqvs <- p.adjust(jpvs, method = 'fdr')
      finald <- tmpres$detected
      finald[W1<0.5*n_otu] <- F
      finald[W1>=detect.rate*n_otu] <- T
      res <- data.frame(log2FC=adj.beta, log2FC.SD=adj.se12, log2FC.pval=jpvs,
                        log2FC.qval=jqvs, Reject.number=W1, REJECT=finald)
    }
    list(marginal=res.m, joint=res.j, final=res)
  }

  if(struc_zero) {
    if(is.null(rownames(Y))) id <- 1:nrow(Y)
    meta_data <- data.frame(ID=id, x=x)
    if(!is.null(Z)) meta_data <- data.frame(ID=id, x=x, Z=Z)
    tmp <- feature_table_pre_process(feature_table=t(Y), meta_data=meta_data,
                                     sample_var='ID', group_var = 'x',
                                     out_cut = 0.05, zero_cut = zero_cut, lib_cut=lib_cut, neg_lb=F)
    Y <- t(as.matrix(tmp$feature_table))
    x <- unlist(meta_data[, 'x'])
    Z <- NULL
    if(!is.null(Z)) Z <- meta_data[, 3:ncol(meta_data)]
    if(struc_zero) {
      stru.zero <- tmp$structure_zeros
      Y <- Y[rowSums(stru.zero)==0, ]
    }
  }

  lx <- length(unique(x))
  n_sample <- length(x)
  if(lx>n_sample/2){
    message('Please make sure the variable of interest is continuous.')
    x <- as.numeric(as.character(x))
    res1 <- .fANCOM(Y=Y, x=x, Z=Z, rand=rand, pseudo=pseudo, outlier=outlier,
                   sig=sig, detect.rate=detect.rate, ref.rate=ref.rate, effect=effect)
    res <- list(global=NULL,  results=res1)
    return(res)
  }else{
    x.label <- unique(x)
    if(lx>2) {
      if(global) {
        message('Global hypothsis among multiple groups was done!')
        x <- as.numeric(as.character(x))
        res1 <- .fANCOM(Y=Y, x=x, Z=Z, rand=rand, pseudo=pseudo, outlier=outlier,
                        sig=sig, detect.rate=detect.rate, ref.rate=ref.rate, effect=effect)
      }else{ # pariwise comparison
        labx <- c(); res2 <- list(); s <- 0
        for(i in 1:(lx-1)) {
          x1.idx <- which(x.label==x.label[i])
          for(j in (i+1):lx) {
            s <- s+1
            x2.idx <- which(x.label==x.label[j])
            ij <- c(x1.idx, x2.idx)
            zij <- NULL
            if(!is.null(Z)) zij <- Z[ij,]
            pc <- paste0(x.label[i], '-', x.label[j])
            res2[s] <- .fANCOM(Y=Y[ij, ], x=x[ij], Z=zij, rand=rand, pseudo=pseudo, outlier=outlier,
                            sig=sig, detect.rate=detect.rate, ref.rate=ref.rate, effect=effect)
            labx <- c(labx, pc)
            message(pc, ' done!')
          }
        }
        names(res2) <- pc
      }
    }else{
      if(lx<2) stop('Please check the value of variable of interest.')
      res1 <- NULL
      res2 <- .fANCOM(Y=Y, x=x, Z=Z, rand=rand, pseudo=pseudo, outlier=outlier,
                      sig=sig, detect.rate=detect.rate, ref.rate=ref.rate, effect=effect)
    }
    res <- list(global=res1, results=res2)
    return(res)
  }
}




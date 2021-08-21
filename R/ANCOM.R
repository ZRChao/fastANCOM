#' @title ANCOM, the traditional method for analysis of composition of microbiomes
#'
#'
#' @description The traditional one uses non-parameteric methods, Wilcoxon or Kruskal-Wallis rank-sum test for two or more groups, to test each joint model directly.
#' And the emperical distribution of the reject numbers are used to define a cutoff for concluding each microbe to be differential or not. This is same function which implements in QIIME2 with python.
#' If there is other covariates, as implement in https://github.com/FrederickHuangLin/ANCOM who used ANOVA. But it could be used fastANCOM to achieve the same goal with the log-linear model.
#'
#'
#' @param Y the abundance table,  count/relative or any smoothed data matrix are supported
#' @param x the variable of interest
#' @param tfun the test function used, including t.test, wilcox.test, kruskal.test, and a speed version for wilcox.test based on the Mann-Whitney U statistic which denote as 't2'
#' @param pseudo the pseudo number for zero smoothing, default is 0.5
#' @param sig the level of significance to count the reject number, default is 0.05
#' @param zero_cut a numerical fraction between 0 and 1. Taxa with proportion of zeroes greater than zero_cut will be excluded in the analysis. Default is 0.90
#' @param lib_cut a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, i.e. do not filter any sample.
#' @param struc_zero whether to detect structural zeros. Default is FALSE
#'
#' @return a dataframe with 3 columns,
#' \itemize{
#' \item Reject.number, the reject number of each microbe
#' \item REJECT, the final result whether each microbe is associated with the variable of interest or not
#' \item detect0.7 hard threshold 0.7 to identify the differential microbes
#' }
#'
#' @examples
#' data <- matrix(rpois(100*60, 10), 60)
#' b <- rnorm(10, 4)
#' data[31:60,1:10] <- t(t(data[31:60,1:10])*2^b)
#' group <- rep(c(0, 1), each=30)
#' fit <- ANCOM(Y=data, x=group)
#' head(fit)
#' @export
#'


ANCOM <- function(Y, x, tfun='t2', pseudo=0.5, sig=0.05,
                  zero_cut=0.9, lib_cut=0, struc_zero=F){

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

  if(struc_zero) {
    if(is.null(rownames(Y))) id <- 1:nrow(Y)
    meta_data <- data.frame(ID=id, x=x)
    tmp <- feature_table_pre_process(feature_table=t(Y), meta_data=meta_data,
                                     sample_var='ID', group_var = 'x',
                                     out_cut = 0.05, zero_cut = zero_cut, lib_cut=lib_cut, neg_lb=F)
    Y <- t(as.matrix(tmp$feature_table))
    x <- unlist(meta_data[, 'x'])
    if(struc_zero) {
      stru.zero <- tmp$structure_zeros
      Y <- Y[rowSums(stru.zero)==0, ]
    }
  }

  n_otu <- ncol(Y)
  n_sample <- nrow(Y)
  logdata <- log2(Y + pseudo)
  logdata[is.infinite(logdata)] <- NA
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  if(class(tfun) == "function"){
    for(ii in 1:(n_otu-1)){
      if(ii%%500 ==0) cat(ii, 'th for ANCOM done!\n')
      data.pair <- logdata[,ii:n_otu]
      lr <- data.pair[,1] - data.pair[,-1]
      logratio.mat[ii,(ii+1):n_otu] <- apply(as.matrix(lr), 2, function(y) tfun(y~x)$p.value)
    }
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
    logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
    W1 <- rowSums(logratio.mat < sig, na.rm = T)
  }
  if(class(tfun) != "function") { # for two groups
    g1 <- unique(x)[1]
    n1 <- table(x)[1]
    n2 <- n_sample - n1
    sig <- qwilcox(sig/2, m = n1, n = n2, lower.tail = F)
    if(tfun == 't1') {
      for(ii in 1:(n_otu-1)) {
        data.pair <- logdata[,ii:n_otu]
        lr <- as.matrix(data.pair[,1] - data.pair[,-1])
        lr1<- as.matrix(lr[x==g1, ])
        lr2<- as.matrix(lr[x!=g1, ])
        u <- rep(NA, ncol(lr))
        for(j in 1:ncol(lr)) {
          x1 <- lr1[, j]
          t2 <- matrix(lr2[, j], nrow=n1, ncol=n2, byrow = T)
          x2 <- sum(x1>t2) + 0.5*sum(x1==t2)
          u[j] <- ifelse(x2>n1*n2/2, x2, n1*n2 - x2)
        }
        # t1 <- matrix(matrix(rep(t(lr1), n2), ncol=ncol(lr), byrow = T), nrow=n1)
        # t2 <- matrix(lr2, nrow=n1, ncol=n2*ncol(lr), byrow = T)
        # t3 <- colSums(t1 > t2) + colSums(t1==t2)*0.5
        # u <- rowSums(matrix(t3, nrow=ncol(lr), byrow = T))
        # u[u<n1*n2/2] <- n1*n2 - u[u<n1*n2/2]
        logratio.mat[ii, (ii+1):n_otu] <- u
      }
    }
    if(tfun == 't2') {
      tmp1 <- logdata[x==g1, ]
      tmp2 <- logdata[x!=g1, ]
      tmp12 <- matrix(NA, n1*n2, n_otu)
      for(ii in 1:n_otu) {
        tmp12[, ii] <- tmp1[, ii]- rep(tmp2[, ii], each=n1)
      }
      logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
      for(ii in 1:(n_otu-1)) {
        ti2 <- as.matrix(tmp12[, (ii+1):n_otu])
        ti <- colSums(tmp12[, ii] > ti2) + colSums(tmp12[, ii] == ti2)*0.5
        ti[ti<n1*n2/2] <- n1*n2 - ti[ti<n1*n2/2]
        logratio.mat[ii, (ii+1):n_otu] <- ti
      }
    }
    if(tfun == 't3') {
      tmp1 <- logdata[x==g1, ]
      tmp2 <- logdata[x!=g1, ]
      tmp12 <- matrix(NA, n1*n2, n_otu)
      tmp12[,1] <- tmp1[, 1]- rep(tmp2[, 1], each=n1)
      for(ii in 2:n_otu) {
        tmp12[, ii] <- tmp1[, ii]- rep(tmp2[, ii], each=n1)
        tii <- as.matrix(tmp12[, ii] - tmp12[,1:(ii-1)])
        tiu <- colSums(tii>0) + colSums(tii==0)*0.5
        tiu[tiu<n1*n2/2] <- n1*n2 - tiu[tiu<n1*n2/2]
        logratio.mat[1:(ii-1), ii] <- tiu
      }
    }
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
    logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- NA
    W1 <- rowSums(logratio.mat > sig, na.rm = T)
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
  res <- decisionMaking(W1)
  res
}













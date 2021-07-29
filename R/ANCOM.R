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


ANCOM <- function(Y, x, tfun='t2', pseudo=0.5, sig=0.05){
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













# print to command line
printf <- function(...)cat(sprintf(...))

# Downloads files from the web and verifies their checksum. Will use local copy in current directory, if it exists
load.web.file <- function(
  url, md5sum, outfile, zipfile = F) {
  # check if local file exists
  if (file.exists(outfile)) {
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) stop(sprintf("Local file %s has wrong checksum: %s", outfile, realsum))
    # do not delete wrong file, it was already here before
    
  } else {
    if(zipfile){
      # download file
      temp <- tempfile()
      download.file(url,temp)
      unzip(zipfile = temp, files = outfile, exdir = ".")
    } else {
      # download file
      download.file(url, outfile)
    }
    # verify checksum
    realsum <- tools::md5sum(outfile)[[1]]
    if (realsum != md5sum) { 
      # delete wrong file
      unlink(outfile)
      stop(sprintf("Remote file %s has wrong checksum: %s", url, realsum))
    }
  }
}

# load adjacency matrices
get_adja <- function(platform){
  
  if (platform=="UPLC") {
    
    adja <- read.csv("data/CRC_adja_extended.csv", header = T, sep = ";", dec=",")
    adja <- adja[,2:dim(adja)[2]]
    rownames(adja) <- colnames(adja)
    
  } else if (platform == "MALDI") {
    
    adja <- read.table("data/adjacency_matrix_sorted.txt", sep = "\t", header = T)
    
  } else if (platform %in%  c("LC-ESI-MS")) {
    
    adja <- read.csv("data/PathwayNEW.csv", header = TRUE, sep = ";", dec = ",",row.names = 1)
    
  }
  
  adja[is.na(adja)] <- 0
  adja[lower.tri(adja, diag = TRUE)] <- NA
  adja
}

# Total Area Normalization
tanorm <- function(X){
  X <- as.matrix(X)
  ta <- X/rowSums(X)
  return(ta)
}

# Median Scaling
mednorm <- function(X){
  X <- as.matrix(X)
  med <- apply(X, 1, median)
  medn <- X - med
}

# Rank Normalization
ranknorm <- function(X, dir = 2){
  
  if(dir == 2) {
    df_rank <- apply(X, 2, rank, ties.method = "min")
  } else{
    df_rank <- t(apply(X, 1, rank, ties.method = "min"))
  }
  return(df_rank)
}

# Quantile Normalization
quantile_normalisation <- function(df, dir = 2){
  
  if(dir==2) {
    df_rank <- apply(df, 2, rank, ties.method = "min")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)
    
    index_to_mean <- function(my_index, my_mean) {
      return(my_mean[my_index])
    }
    
    df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
    
  } else{
    df_rank <- t(apply(df, 1, rank, ties.method = "min"))
    df_sorted <- t(data.frame(apply(df, 1, sort)))
    df_mean <- apply(df_sorted, 2, mean)
    
    index_to_mean <- function(my_index, my_mean) {
      return(my_mean[my_index])
    }
    
    df_final <- t(apply(df_rank, 1, index_to_mean, my_mean = df_mean))
    
  }
  
  return(df_final)
}

# Probabilistic Quotient Normalization
quotNorm <- function(X, vars=1:dim(X)[2], NAerror=F) {
  
  # check if there are any NAs
  if (sum(is.na(X[,vars]))>0) {
    # throw warning or error?
    if (NAerror) {
      stop('Data matrix contains NAs')
    } else {
      warning('Data matrix contains NAs')
    }
  }
  
  # median reference sample
  ref <- apply(X[,vars],2,function(x)median(x,na.rm=T))
  # get dilution factors
  d <- apply(X[,vars],1,  function(s) median(as.numeric(s/ref),na.rm=T))
  # apply to each sample  (for each row=sample, divide values by median dilution factor)
  Y <- t(sapply(1:dim(X)[1], function(i)X[i,]/d[i]))
  Y <- apply(Y, 2, unlist) %>% as.data.frame()
  
  # return
  list(X=Y,dilution=d)
}

# Adaptation of GeneNet (version  1.2.15) package code to run silently
# source code taken directly from package
# (removed hardcoded, unconditional calls of cat function)
network.test.edges.silent <- function (r.mat, fdr = TRUE, direct = FALSE, plot = TRUE, ...) 
{
  pcor = sm2vec(r.mat)
  indexes = sm.index(r.mat)
  colnames(indexes) = c("node1", "node2")
  w = cbind(pcor, indexes)
  if (fdr == TRUE) {
    # cat("Estimate (local) false discovery rates (partial correlations):\n")
    fdr.out = fdrtool(w[, 1], statistic = "correlation", 
                      plot = plot, ...)
    pval = fdr.out$pval
    qval = fdr.out$qval
    prob = 1 - fdr.out$lfdr
  }
  else {
    pval = rep(NA, length(w[, 1]))
    qval = pval
    prob = pval
  }
  result = cbind(w, pval, qval, prob)
  if (direct == TRUE) {
    spvar = attr(r.mat, "spv")
    if (is.null(spvar)) {
      r.mat.cor = pcor2cor(r.mat)
      spvar = 1/diag(solve(r.mat.cor))
    }
    p = length(spvar)
    r.spvar = (t(spvar %*% t(rep(1, p)))/(spvar %*% t(rep(1, 
                                                          p))))
    log.spvar = log(sm2vec(r.spvar))
    if (fdr == TRUE) {
      if (plot == TRUE) {
        dev.new()
      }
      # cat("Estimate (local) false discovery rates (log ratio of spvars):\n")
      fdr.out = fdrtool(log.spvar, statistic = "normal", 
                        plot = plot, ...)
      pval.dir = fdr.out$pval
      qval.dir = fdr.out$qval
      prob.dir = 1 - fdr.out$lfdr
    }
    else {
      pval.dir = rep(NA, length(w[, 1]))
      qval.dir = pval.dir
      prob.dir = pval.dir
    }
    result = cbind(result, log.spvar, pval.dir, qval.dir, 
                   prob.dir)
  }
  sort.idx = order(-abs(result[, 1]))
  result = as.data.frame(result[sort.idx, ])
  return(result)
}

# Create contingency table
contab <- function(known, corr) {
  ct_a <-
    length(which(known == 1 & corr == 1)) #true positive, upper left
  ct_b <-
    length(which(known == 0 & corr == 1)) #false positive, upper right
  ct_c <-
    length(which(known == 1 & corr == 0)) #false negative, lower left
  ct_d <-
    length(which(known == 0 & corr == 0)) #true negative, lower right
  ct <- cbind(c(ct_a, ct_c), c(ct_b, ct_d))
  return(ct)
}

# Compute 95% confidence interval of the values
confin <- function(fisp) {
  tmean <- mean(fisp,na.rm = T)
  tupper <- quantile(fisp, probs = (0.025), na.rm = T)
  tlower <- quantile(fisp, probs = (0.975), na.rm = T)
  conf <- cbind(tupper, tmean, tlower)
}

# create ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# revert list structure
revert_list_str_4 <- function(ls) { 
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list) 
}
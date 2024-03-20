## generateSummaryStats accepts a vector of numeric values and uses it to
## generate several different summary statistics
generateSummaryStats <- function(var_data, varname = "") {
  
  tor = data.frame(mean   = mean(var_data),
                   sd     = sd(var_data),
                   median = median(var_data),
                   min    = min(var_data),
                   max    = max(var_data),
                   q25    = quantile(var_data, p = 0.25),
                   q75    = quantile(var_data, p = 0.75) )
  
  names(tor) <- paste0(varname, "_", names(tor))
  
  tor
  
}

## localMaximumApproach finds a local maximum mean.pred is vector sequence of
## prediction strength metrics at each k is the (one-sided) width of the search
## space for local maximum function returns the first peak identified if no
## peaks found, 0 will be returned
localMaximumApproach = function (mean.pred, m = 3){
  shape = diff( sign( diff( mean.pred, na.pad = FALSE ) ) )
  if ( sum( shape < 0 ) > 0 )  {
    pks = sapply( which( shape < 0 ), FUN = function(i){
      z = i - m + 1
      z = ifelse(z > 0, z, 1)
      w = i + m + 1
      w = ifelse( w < length(mean.pred), w, length(mean.pred) )
      if(all( mean.pred[c(z : i, (i + 2) : w)] <= mean.pred[i + 1])) return(i + 1) else return(0)
    })
    pks = unlist(pks)
  }
  else
  {
    pks = 0
  }
  return(pks)
}

# maxDiag maximizes the sum of the elements in the diagonal of a square matrix.
# Currently, function does not check whether the matrix is square.
maxDiag <- function(x) {
  n <- ncol(x)
  per <- gtools::permutations(n,n)
  d <- apply(per,1,function(y) sum(diag(x[,y])))
  return(x[,per[which.max(d),]])
}

# av_ts_designed creates an average time series for each designed cluster
avg_ts_designed <- function(wide_ts,sim_clusters){
  wide_clusters <- merge(sim_clusters,wide_ts,by=c("id","day"))
  
  wide_clusters %<>% dplyr::select(-observed_cluster)
  cluster_groupsplit <- dplyr::group_split(wide_clusters %>% 
                                             dplyr::group_by(designed_cluster))
  
  avg_ts <- c()
  rownames_avg_ts <- c()
  for (j in 1:length(cluster_groupsplit)){
    avg_ts <- rbind(avg_ts,colMeans(cluster_groupsplit[[j]][,4:99]))
    rownames_avg_ts <- rbind(rownames_avg_ts,cluster_groupsplit[[j]][1,"designed_cluster"])
  }
  
  rownames(avg_ts) <- rownames_avg_ts$designed_cluster
  return(avg_ts)
}

# av_ts_estimated creates an average time series for each estimated cluster
avg_ts_estimated <- function(wide_ts){
  cluster_groupsplit <- dplyr::group_split(wide_ts %>% 
                                             dplyr::group_by(observed_cluster))
  
  avg_ts <- c()
  rownames_avg_ts <- c()
  for (j in 1:length(cluster_groupsplit)){
    avg_ts <- rbind(avg_ts,colMeans(cluster_groupsplit[[j]][,3:98]))
    rownames_avg_ts <- rbind(rownames_avg_ts,cluster_groupsplit[[j]][1,"observed_cluster"])
  }
  
  rownames(avg_ts) <- rownames_avg_ts$observed_cluster
  return(avg_ts)
}

# dist_function creates a distance matrix between the average designed and
# estimated clusters. The function allows the user to choose amongst the following
# distance metrics: chevyshev, cosine, cross correlation, euclidean, and manhattan
dist_function <- function(ts_observed,ts_estimated,
                          variant=c("chebyshev","cosine","cross_correlation",
                                    "euclidean","manhattan")){
  dist_df <- as.data.frame(matrix(data=NA,
                                  nrow=nrow(ts_observed),
                                  ncol=nrow(ts_observed)))
  rownames(dist_df) <- rownames(ts_observed)
  colnames(dist_df) <- rownames(ts_estimated)
  
  for (j in 1:nrow(dist_df)){
    for (k in 1:ncol(dist_df)){
      dist_df[j,k] <- switch(variant,
                             chebyshev = abdiv::chebyshev(ts_observed[j,],ts_estimated[k,]),
                             cosine = 1-lsa::cosine(ts_observed[j,],ts_estimated[k,]),
                             cross_correlation = TSdist::CCorDistance(ts_observed[j,],ts_estimated[k,]),
                             euclidean = TSdist::EuclideanDistance(ts_observed[j,],ts_estimated[k,]),
                             manhattan = abdiv::manhattan(ts_observed[j,],ts_estimated[k,]))
    }
  }
  return(dist_df)
}

# resample implements a multi-stage bootstrap recursively. That is, levels of
# hierarchy are traversed by nested calls to resample. The dat argument is a
# dataframe with factor fields for each level of hierarchy (e.g., id, day),
# and a numeric field of measured values. The cluster argument is a character
# vector that identifies the hierarchy in order from top to bottom
# (e.g., c('id','day')). The replace argument is a logical vector that indicates
# whether sampling should be with replacement at the corresponding level of
# hierarchy (e.g., c(TRUE,FALSE)).
resample <- function(dat, cluster, replace) {
  
  # exit early for trivial data
  if(nrow(dat) == 1 || all(replace==FALSE))
    return(dat)
  
  # sample the clustering factor
  cls <- sample(unique(dat[[cluster[1]]]), replace=replace[1])
  
  # subset on the sampled clustering factors
  sub <- lapply(cls, function(b) subset(dat, dat[[cluster[1]]]==b))
  
  # sample lower levels of hierarchy (if any)
  if(length(cluster) > 1)
    sub <- lapply(sub, resample, cluster=cluster[-1], replace=replace[-1])
  
  # join and return samples
  do.call(rbind, sub)
  
}

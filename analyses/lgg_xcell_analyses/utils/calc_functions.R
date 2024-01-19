suppressPackageStartupMessages({
  library(diptest)
})

get_gene_list_by_diptest <- function(count_matrix, min_n=1000){
  # define a dataframe to store the pval
  count_matrix_pval <- data.frame()
  for(j in 1:nrow(count_matrix)){
    # calculate dip.test for each gene
    diptest_each <- count_matrix[i,] %>% 
      as.numeric() %>%
      dip.test(simulate.p.value = FALSE, B = 2000)
    # gather the pvalue for each gene
    pval_each <- diptest_each$p.value
    # add another column to store pval
    count_matrix_pval_each <- count_matrix[i,] %>%
      mutate(pval = pval_each)
    # cocgine each line back to the dataframe 
    count_matrix_pval <- bind_rows(count_matrix_pval,  count_matrix_pval_each)
  }
  
  # define rownames
  rownames(count_matrix_pval) <- rownames(count_matrix)
  
  # Filter to expression that only has pval<0.05 if over min_n genes satisfy this or take top min_n
  # see how many genes have pval<0.5
  n <- count_matrix_pval %>% 
    dplyr::filter(pval <0.05) %>%
    nrow()
  
  if(n>min_n){
    count_matrix_pval <-count_matrix_pval %>% 
      dplyr::filter(pval <0.05) %>% 
      dplyr::select(-pval)
  } else {
    count_matrix_pval <-count_matrix_pval %>% 
      dplyr::arrange(pval, descending = FALSE) %>% 
      head(min_n) %>%
      dplyr::select(-pval)
  }
  return(rownames(count_matrix_pval))
}


# function to compute v.test
compute.v.test <- function(x, clustering_col){
  
  # get unique clusters i.e. 1-n clusters
  clusters <- unique(x[[clustering_col]])
  
  out <- data.frame()
  # iterate over each cluster
  for(i in 1:length(clusters)){
    
    # subset x to cluster i
    y = x %>% 
      dplyr::filter(get(clustering_col) == clusters[i])
    
    # mean score per cluster
    mean.celltype.cluster = unique(y$cluster_celltype_mean_fraction)
    
    # mean score overall (i.e. global mean)
    mean.fraction = unique(y$celltype_mean_fraction)
    
    # calculate numerator
    num = mean.celltype.cluster - mean.fraction
    
    # total sample size
    n = nrow(x)
    
    # cluster sample size
    ng = nrow(y)
    
    # variance of cell type score (i.e. global variance)
    var.celltype = unique(x$celltype_variance_fraction) 
    
    # calculate denominator
    denom = (n-ng/n-1)*(var.celltype/ng)
    denom = sqrt(denom)
    
    # calculate vscore
    v = num/denom
    out[i,'cluster'] <- clusters[i]
    out[i,'v_score'] <- v
  }
  return(out)
}

###### function that does zscore transformation 
zscore_transform <- function(df){
  df_log <- log2(df+1)
  # calculate mean
  matrix_means <-rowMeans(df_log, na.rm = TRUE)
  # calculate sd
  matrix_sd <- apply(df_log, 1, sd, na.rm = TRUE)
  # subtract mean
  df_log_minus_mean <- sweep(df_log, 1, matrix_means, FUN = "-")
  # divide by SD remove NAs and Inf values from zscore for genes with 0 in normData
  df_zscored <- sweep(df_log_minus_mean, 1,matrix_sd, FUN = "/") %>% 
    na_if(Inf) %>% na.omit()
  
  return(df_zscored)
}

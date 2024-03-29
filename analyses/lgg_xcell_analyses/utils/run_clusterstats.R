#' Function: Calculate cluster stats for each k
#'
#' @author Komal S. Rathi
#' 
#' import dplyr
#' import IntNMF
#'
#' @param dat a list of matrices with features as rows and samples as columns
#' @param wt weight assigned to each matrix or modality
#' @param output_dir output directory to write results
#' @param k_value max value of k to be evaluated
#'
#' @return a list containing information on the best fit "k" and corresponding cluster assigments
#'
#' @export

run_clusterstats <- function(dat, wt = NULL, output_dir, k_value){
  
  # if weight is not assigned, use default
  if(is.null(wt)){
    wt = if(is.list(dat)) rep(1,length(dat)) else 1
  }
  
  # do this for each cluster
  # run Nonnegative Matrix Factorization of Multiple data using Nonnegative Alternating Least Square
  nmf_fname <- file.path(output_dir, "intnmf_fit_all.rds")
  if(!file.exists(nmf_fname)){
    nmf_output <- list() # nmf mnnals output
    for(i in 2:k_value){
      print(i)
      var <- paste0("k", i)
      nmf_output[[var]] <- nmf.mnnals(dat = dat, 
                                      maxiter = 200, # default
                                      st.count = 20, # default
                                      wt = wt,
                                      k = i, 
                                      seed = TRUE)
    }
    # save final nmf output for all k-values
    saveRDS(nmf_output, file = nmf_fname)
  } else {
    nmf_output <- readRDS(nmf_fname)
  }
  
  # do for each cluster
  output_fname <- file.path(output_dir, "intnmf_clusterstats.tsv")
  output_df <- data.frame() # output df for clusterstats
  for(j in 1:length(dat)){
    print(j)
    input_dat <- dat[[j]]
    for(i in 2:k_value){
      print(i)
      var <- paste0("k", i)
      df <- data.frame(sample_id = names(nmf_output[[var]]$clusters), cluster = nmf_output[[var]]$clusters)
      count_matrix_dist <- factoextra::get_dist(input_dat, method = 'pearson')
      count_matrix_dist[is.na(count_matrix_dist)] <- 0
      clus_stats <- fpc::cluster.stats(d = count_matrix_dist, clustering = df$cluster)
      clus_stats$avg_sil <- mean(clus_stats$clus.avg.silwidths) # calculate mean 
      clus_stats <- lapply(clus_stats, FUN = function(x){
        if(length(x) > 1 & !is.null(x)) {
          x <- toString(round(x, 2))
        } else {
          x <- round(as.numeric(x), digits = 2)
        }
      })
      clus_stats <- t(as.data.frame(unlist(clus_stats)))
      clus_stats <- clus_stats %>%
        as.data.frame() %>%
        dplyr::mutate(k = var,
                      dataset = j) %>%
        dplyr::select(k, dataset, cluster.size, sindex, noisen, diameter, average.distance, median.distance, average.between, average.within, within.cluster.ss, clus.avg.silwidths, avg_sil, dunn, dunn2, entropy, wb.ratio)
      # combine with other k values
      output_df <- rbind(output_df, clus_stats)
    }
  }
  
  # save dataframe with cluster stats for all k-values
  write_tsv(output_df, file = output_fname)
  
  # # read cpi output for all k
  # cpi <- readRDS(file = file.path(output_dir, "cpi_output.rds"))
  
  # cluster stats output, compute max k by taking the difference between average between and average within
  # max_k = output_df[which.max(as.numeric(output_df$average.between)-as.numeric(output_df$average.within)),] %>% pull(k)
  # sum the average.between and average.within across the modes of data at each k, 
  # then take the difference of those two summed values at each k to select the optimal cluster number
  max_k <- output_df %>% 
    group_by(k) %>% 
    dplyr::summarise(n = sum(as.numeric(average.between)) - sum(as.numeric(average.within))) %>%
    filter(n == max(n)) %>%
    pull(k)
  
  # # if there are ties, choose the one with max CPI
  # if(length(max_k) > 1){
  #   cpi_subset <- cpi[rownames(cpi) %in% max_k,]
  #   max_k <- names(which.max(apply(cpi_subset, MARGIN = 1, mean)))
  # }
  
  # if there are ties, choose the one with max intNMF's silhouette index
  sil_score_df <- data.frame()
  if(length(max_k) > 1){
    for(i in 1:length(max_k)){
      fit <- nmf_output[[max_k[i]]]
      distance <- 1 - fit$consensus
      summary <- cluster::silhouette(fit$clusters, dmatrix = distance) %>% summary
      average_width <- summary$avg.width
      tmp <- data.frame(k = max_k[i], sil_score = average_width)
      sil_score_df <- rbind(sil_score_df, tmp)
    }
    # get the k with max average silhouette width
    max_k <- sil_score_df[which.max(sil_score_df$sil_score),"k"]
  }
  
  # add modality names instead of H1, H2... for better clarity in downstream plotting
  names(nmf_output[[max_k]]$H) <- names(dat)
  
  # save best fit to file
  saveRDS(object = nmf_output[[max_k]], file = file.path(output_dir, "intnmf_best_fit.rds"))
  
  # return the nmf output corresponding to the most optimal k for downstream analyses
  return(nmf_output[[max_k]])
}

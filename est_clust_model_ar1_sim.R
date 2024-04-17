est_clust_model_ar1_sim <- function(feat_dur,feat_label,feat_stat,num_clust,
                                    dist_variant,beta_vector,
                                    ar1_phi,ar1_sig_z){
  ## Grabbing file names
  sim_files <- list.files(path="/Users/ajensen/Library/CloudStorage/Box-Box/Missing Data in Digital Health Studies/Clustering Manuscript/simulated_data/c7521802f265555c61be8d28687f4b56_raw_cohort/",
                          all.files=FALSE,full.names=TRUE)
  
  ## AR(1) error structure
  sig_eps <- ar1_sig_z/(1-(ar1_phi^2))
  Sigma <- sig_eps * matrix(c(1,ar1_phi,ar1_phi^2,ar1_phi^3,ar1_phi^4,
                              ar1_phi,1,ar1_phi,ar1_phi^2,ar1_phi^3,
                              ar1_phi^2,ar1_phi,1,ar1_phi,ar1_phi^2,
                              ar1_phi^3,ar1_phi^2,ar1_phi,1,ar1_phi,
                              ar1_phi^4,ar1_phi^3,ar1_phi^2,ar1_phi,1),5,5)
  mu <- c(0,0,0,0,0)
  
  ## Creating the data frame to store the fitted parameters
  estlabel_betas <- data.frame(matrix(data=NA,nrow=length(sim_files),ncol=22))
  colnames(estlabel_betas) <- c("Cohort","prop_correct",
                                "clust1_correct","clust2_correct","clust3_correct",
                                "clust4_correct","clust5_correct","clust6_correct",
                                "phi","MSE",
                                "Beta_0","Beta_Clust2","Beta_Clust3",
                                "Beta_Clust4","Beta_Clust5","Beta_Clust6",
                                "Beta_0_SE","Beta_Clust2_SE","Beta_Clust3_SE",
                                "Beta_Clust4_SE","Beta_Clust5_SE","Beta_Clust6_SE")
  estlabel_betas$Cohort <- seq(from=1,to=length(sim_files),by=1)
  
  ## Simulation loop
  for (i in 1:length(sim_files)){
    # reading the file
    mims_data <- readRDS(sim_files[i])
    
    # slightly reformat the data and create a vector magnitude column
    mims_data_aug <- mims_data %>%
      mutate(timestamp = HEADER_TIME_STAMP,
             date = day,
             bucket_index = feature_bucket,
             vecmag = sqrt(MIMS_UNIT_X^2 + MIMS_UNIT_Y^2 + MIMS_UNIT_Z^2))
    
    # collapse second-epoch data into feature windows
    bucketed_data <- mims_data_aug %>%
      split( f = list(.$id, .$day, .$bucket_index) ) %>%
      map_dfr(function(bucket){
        
        tor <- data.frame(id = bucket$id[1],
                          day = bucket$day[1],
                          designed_cluster = bucket$cluster[1],
                          # pool_id   = bucket$pool_id[1],
                          # pool_date = bucket$pool_date[1],
                          # pool_bucket   = bucket$bucket_index[1],
                          bucket_width_seconds = feat_dur,
                          bucket_number = bucket$bucket_index[1],
                          n_seconds = nrow(bucket))
        
        x_stats      <- generateSummaryStats(bucket$MIMS_UNIT_X, paste0(feat_label, "_x") )
        y_stats      <- generateSummaryStats(bucket$MIMS_UNIT_Y, paste0(feat_label, "_y"))
        z_stats      <- generateSummaryStats(bucket$MIMS_UNIT_Z, paste0(feat_label, "_z"))
        vecmag_stats <- generateSummaryStats(bucket$vecmag,      paste0(feat_label, "_vecmag"))
        
        tor <- bind_cols(tor,
                         x_stats,
                         y_stats,
                         z_stats,
                         vecmag_stats) %>%
          mutate(source_file = sim_files[i])
        
        row.names(tor) <- NULL
        
        tor
        
      })

    # note simulated cluster membership
    simulated_clusters <- bucketed_data %>%
      dplyr::select(id, day, designed_cluster) %>%
      distinct
    
    # transform buckets to wide
    wide_buckets <- bucketed_data %>%
      mutate(bucket_day = bucket_number) %>%
      dplyr::select(id,
                    day,
                    bucket_day,
                    matches(feat_stat)) %>%
      dcast(id + day ~ bucket_day)
    
    # transform into a matrix for prediction.strength()
    wide_buckets_matrix <- wide_buckets %>%
      dplyr::select( -id, -day) %>%
      as.matrix
    
    # generate pam model using a priori known number of clusters
    optimal_clusters <- pam(x = wide_buckets_matrix,
                            k = num_clust,
                            diss = FALSE,
                            metric = "euclidean")
    
    # merge pam cluster assignments back in with wide data
    wide_buckets$observed_cluster <- optimal_clusters$clustering
    
    cluster_match <- wide_buckets %>%
      merge( y = simulated_clusters ,
             by = c("id", "day")) %>%
      dplyr::select(id, day, designed_cluster, observed_cluster)
    
    # making the average time series for designed and estimated clusters
    avg_ts_design <- avg_ts_designed(wide_buckets,simulated_clusters)
    avg_ts_estimate <- avg_ts_estimated(wide_buckets)
    
    # mapping the designed cluster names to the estimated clusters based on a
    # distance function pre-specified by the user
    # des_est_dist <- dist_function(avg_ts_design,avg_ts_estimate,
    #                               variant=dist_variant)
    # des_est_mapping <- as.data.frame(cbind(rownames(des_est_dist),colnames(des_est_dist)[apply(des_est_dist, 1, which.min)])) %>% rename(Designed_Cluster=V1,Estimated_Cluster=V2)
    
    # mapping the designed cluster names to the estimated clusters based on a
    # Hungarian matching algorithm using a distance funciton pre-specified by
    # the user 
    des_est_dist <- dist_function(avg_ts_design,avg_ts_estimate,
                                  variant=dist_variant)
    
    hung_match <- clue::solve_LSAP(as.matrix(des_est_dist), maximum=FALSE)
    
    des_est_mapping <- as.data.frame(cbind(rownames(des_est_dist),colnames(des_est_dist)[hung_match[1:num_clust]])) %>% rename(Designed_Cluster=V1,Estimated_Cluster=V2)
    
    # adding the mapped cluster names to the cluster_match dataset
    cluster_match$obs_cluster_label <- dplyr::case_when(cluster_match$observed_cluster==1 ~ des_est_mapping[des_est_mapping$Estimated_Cluster==1,]$Designed_Cluster,
                                                        cluster_match$observed_cluster==2 ~ des_est_mapping[des_est_mapping$Estimated_Cluster==2,]$Designed_Cluster,
                                                        cluster_match$observed_cluster==3 ~ des_est_mapping[des_est_mapping$Estimated_Cluster==3,]$Designed_Cluster,
                                                        cluster_match$observed_cluster==4 ~ des_est_mapping[des_est_mapping$Estimated_Cluster==4,]$Designed_Cluster,
                                                        cluster_match$observed_cluster==5 ~ des_est_mapping[des_est_mapping$Estimated_Cluster==5,]$Designed_Cluster,
                                                        cluster_match$observed_cluster==6 ~ des_est_mapping[des_est_mapping$Estimated_Cluster==6,]$Designed_Cluster)
    
    # getting proportion of correct predicted clusters
    cont_table <- as.matrix(table(cluster_match$designed_cluster,
                                  cluster_match$obs_cluster_label))
    
    cont_table_sorted <- cont_table[rownames(cont_table),]
    
    estlabel_betas[i,"prop_correct"] <- sum(diag(cont_table_sorted))/sum(cont_table_sorted)
    estlabel_betas[i,"clust1_correct"] <- diag(cont_table_sorted)["0-24 sedentary"]/rowSums(cont_table_sorted)["0-24 sedentary"]
    estlabel_betas[i,"clust2_correct"] <- diag(cont_table_sorted)["0-6 active"]/rowSums(cont_table_sorted)["0-6 active"]
    estlabel_betas[i,"clust3_correct"] <- diag(cont_table_sorted)["6-12 active"]/rowSums(cont_table_sorted)["6-12 active"]
    estlabel_betas[i,"clust4_correct"] <- diag(cont_table_sorted)["12-18 active"]/rowSums(cont_table_sorted)["12-18 active"]
    estlabel_betas[i,"clust5_correct"] <- diag(cont_table_sorted)["18-24 active"]/rowSums(cont_table_sorted)["18-24 active"]
    estlabel_betas[i,"clust6_correct"] <- diag(cont_table_sorted)["0-24 active"]/rowSums(cont_table_sorted)["0-24 active"]
    
    # grabbing estimated cluster membership
    estimated_clusters <- cluster_match %>%
      dplyr::select(id, day, obs_cluster_label) %>%
      distinct
    
    # data management of the estimated cluster memberships
    dat <- estimated_clusters %>%
      arrange(id,day) 
    
    dat <- merge(dat,simulated_clusters,by=c("id","day"))
    
    dat %<>% mutate(cluster1 = ifelse(obs_cluster_label=="0-24 sedentary",1,0),
                    cluster2 = ifelse(obs_cluster_label=="0-6 active",1,0),
                    cluster3 = ifelse(obs_cluster_label=="6-12 active",1,0),
                    cluster4 = ifelse(obs_cluster_label=="12-18 active",1,0),
                    cluster5 = ifelse(obs_cluster_label=="18-24 active",1,0),
                    cluster6 = ifelse(obs_cluster_label=="0-24 active",1,0)) %>%
      mutate(des_cluster1 = ifelse(designed_cluster=="0-24 sedentary",1,0),
             des_cluster2 = ifelse(designed_cluster=="0-6 active",1,0),
             des_cluster3 = ifelse(designed_cluster=="6-12 active",1,0),
             des_cluster4 = ifelse(designed_cluster=="12-18 active",1,0),
             des_cluster5 = ifelse(designed_cluster=="18-24 active",1,0),
             des_cluster6 = ifelse(designed_cluster=="0-24 active",1,0))
    
    n_subj = length(unique(dat$id)) #number of subjects
    n_days = length(unique(dat$day)) #number of days each subject is observed

    # simulated AR(1) error structure
    ar1_error <- as.data.frame(MASS::mvrnorm(n_subj,mu,Sigma))
    ar1_error %<>% mutate(id=unique(dat$id)) %>%
      relocate(id,.before="V1")
    
    ar1_error_long <- ar1_error %>%
      pivot_longer(cols=`V1`:`V5`,
                   names_to="V_Num",
                   values_to="AR1_Value") %>%
      mutate(day=dat$day) %>%
      dplyr::select(-c("V_Num")) %>%
      relocate(day,.after="id")
    
    # simulating the triglyceride outcome
    dat$TG_nodecomp <- beta_vector[1] + 
      beta_vector[2]*dat$des_cluster2 +
      beta_vector[3]*dat$des_cluster3 +
      beta_vector[4]*dat$des_cluster4 +
      beta_vector[5]*dat$des_cluster5 +
      beta_vector[6]*dat$des_cluster6 +
      ar1_error_long$AR1_Value
    
    # mixed model
    mixed_nodecomp <- nlme::gls(TG_nodecomp ~ cluster2 + cluster3 + cluster4 + 
                                  cluster5 + cluster6,
                                data=dat,corr=corAR1(, form= ~ 1 | id))
    
    fixed_effects <- as.numeric(broom.mixed::tidy(mixed_nodecomp)$estimate)
    
    # adding the phi estimate to the dataframe
    estlabel_betas[i,"phi"] <- coef(mixed_nodecomp$modelStruct$corStruct,unconstrained=FALSE)
    estlabel_betas[i,"MSE"] <- mean(mixed_nodecomp$residuals^2)
    
    # adding the parameter estimate to the dataframe
    for (j in 1:length(fixed_effects)){
      estlabel_betas[i,j+10] <- broom.mixed::tidy(mixed_nodecomp)[j,"estimate"]
      estlabel_betas[i,j+(10+length(fixed_effects))] <- broom.mixed::tidy(mixed_nodecomp)[j,"std.error"]
    }
  }
  return(estlabel_betas)
}
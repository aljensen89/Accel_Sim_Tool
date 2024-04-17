noisy_true_model_iid_sim <- function(beta_vector,iid_sigma){
  ## Grabbing file names
  sim_files <- list.files(path="/Users/ajensen/Library/CloudStorage/Box-Box/Missing Data in Digital Health Studies/Clustering Manuscript/simulated_data/1ee63015eba36f0e4aa00549233089f3_raw_cohort/",
                          all.files=FALSE,full.names=TRUE)
  
  ## Creating the data frame to store the fitted parameters
  truelabel_betas <- data.frame(matrix(data=NA,nrow=length(sim_files),ncol=15))
  colnames(truelabel_betas) <- c("Cohort","iid_sigma","MSE",
                                 "Beta_0","Beta_Clust2","Beta_Clust3",
                                 "Beta_Clust4","Beta_Clust5","Beta_Clust6",
                                 "Beta_0_SE","Beta_Clust2_SE","Beta_Clust3_SE",
                                 "Beta_Clust4_SE","Beta_Clust5_SE","Beta_Clust6_SE")
  truelabel_betas$Cohort <- seq(from=1,to=length(sim_files),by=1)
  
  ## Simulation loop
  for (i in 1:length(sim_files)){
    # reading the file
    mims_data <- readRDS(sim_files[i])
    
    # grabbing simulated cluster membership
    simulated_clusters <- mims_data %>%
      dplyr::select(id, day, cluster) %>%
      distinct
    
    # data management of the simulated cluster memberships
    dat <- simulated_clusters %>%
      arrange(id,day) 
    
    dat %<>% mutate(cluster1 = ifelse(cluster=="0-24 sedentary",1,0),
                    cluster2 = ifelse(cluster=="0-6 active",1,0),
                    cluster3 = ifelse(cluster=="6-12 active",1,0),
                    cluster4 = ifelse(cluster=="12-18 active",1,0),
                    cluster5 = ifelse(cluster=="18-24 active",1,0),
                    cluster6 = ifelse(cluster=="0-24 active",1,0))
    
    n_subj = length(unique(dat$id)) #number of subjects
    n_days = length(unique(dat$day)) #number of days each subject is observed
    
    # simulated iid error structure
    eps = rnorm(n = n_subj*n_days, mean = 0, sd = iid_sigma)
    
    # simulating the triglyceride outcome
    dat$TG_nodecomp <- beta_vector[1] + 
      beta_vector[2]*dat$cluster2 +
      beta_vector[3]*dat$cluster3 +
      beta_vector[4]*dat$cluster4 +
      beta_vector[5]*dat$cluster5 +
      beta_vector[6]*dat$cluster6 +
      eps
    
    # linear model
    linear_nodecomp<- lm(TG_nodecomp ~ cluster2 + cluster3 + cluster4 + 
                           cluster5 + cluster6, data=dat)
    
    lm_results <- nrow(broom::tidy(linear_nodecomp))
    
    # adding the random intercept and phi estimates to the dataframe
    truelabel_betas[i,"iid_sigma"] <- summary(linear_nodecomp)$sigma
    truelabel_betas[i,"MSE"] <- mean(linear_nodecomp$residuals^2)
    
    # adding the parameter estimate to the dataframe
    for (j in 1:lm_results){
      truelabel_betas[i,j+3] <- broom::tidy(linear_nodecomp)[j,"estimate"]
      truelabel_betas[i,j+(3+lm_results)] <- broom::tidy(linear_nodecomp)[j,"std.error"]
    }
  }
  return(truelabel_betas)
}
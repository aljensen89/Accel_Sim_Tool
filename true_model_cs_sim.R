true_model_cs_sim <- function(beta_vector,cs_rho){
  ## Grabbing file names
  sim_files <- list.files(path="/Users/ajensen/Library/CloudStorage/Box-Box/Missing Data in Digital Health Studies/Clustering Manuscript/simulated_data/c7521802f265555c61be8d28687f4b56_raw_cohort/",
                          all.files=FALSE,full.names=TRUE)
  
  ## CS error structure
  Sigma <- matrix(c(1,cs_rho,cs_rho,cs_rho,cs_rho,
                    cs_rho,1,cs_rho,cs_rho,cs_rho,
                    cs_rho,cs_rho,1,cs_rho,cs_rho,
                    cs_rho,cs_rho,cs_rho,1,cs_rho,
                    cs_rho,cs_rho,cs_rho,cs_rho,1),5,5)
                    
  mu <- c(0,0,0,0,0)
  
  ## Creating the data frame to store the fitted parameters
  truelabel_betas <- data.frame(matrix(data=NA,nrow=length(sim_files),ncol=15))
  colnames(truelabel_betas) <- c("Cohort","cs_rho","MSE",
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
    
    # simulated CS error structure
    cs_error <- as.data.frame(MASS::mvrnorm(n_subj,mu,Sigma))
    cs_error %<>% mutate(id=unique(dat$id)) %>%
      relocate(id,.before="V1")
    
    cs_error_long <- cs_error %>%
      pivot_longer(cols=`V1`:`V5`,
                   names_to="V_Num",
                   values_to="CS_Value") %>%
      mutate(day=dat$day) %>%
      dplyr::select(-c("V_Num")) %>%
      relocate(day,.after="id")
    
    # simulating the triglyceride outcome
    dat$TG_nodecomp <- beta_vector[1] + 
      beta_vector[2]*dat$cluster2 +
      beta_vector[3]*dat$cluster3 +
      beta_vector[4]*dat$cluster4 +
      beta_vector[5]*dat$cluster5 +
      beta_vector[6]*dat$cluster6 +
      cs_error_long$CS_Value
    
    # mixed model
    mixed_nodecomp<- nlme::gls(TG_nodecomp ~ cluster2 + cluster3 + cluster4 + 
                                 cluster5 + cluster6,
                               method="REML",data=dat,
                               correlation=corCompSymm(form=~1|id))
    
    fixed_effects <- nrow(broom.mixed::tidy(mixed_nodecomp))
    
    # adding the random intercept and phi estimates to the dataframe
    truelabel_betas[i,"cs_rho"] <- coef(mixed_nodecomp$modelStruct$corStruct,unconstrained=FALSE)
    truelabel_betas[i,"MSE"] <- mean(mixed_nodecomp$residuals^2)
    
    # adding the parameter estimate to the dataframe
    for (j in 1:fixed_effects){
      truelabel_betas[i,j+3] <- broom.mixed::tidy(mixed_nodecomp)[j,"estimate"]
      truelabel_betas[i,j+(3+fixed_effects)] <- broom.mixed::tidy(mixed_nodecomp)[j,"std.error"]
    }
  }
  return(truelabel_betas)
}
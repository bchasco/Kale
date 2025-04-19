#---------------------------------------------------------------------------------------------------------- -
#                                                   ----
#---------------------------------------------------------------------------------------------------------- -
# install_or_load_pack <- function(packages) {
#   new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
#   if (length(new_packages)) {
#     install.packages(new_packages, dependencies = TRUE)
#   }
  # sapply(packages, require, character.only = TRUE)
# }

# Load functions
# sapply(FUN = source, paste(getwd(), "functions", list.files("functions"), sep="/"))
logit<-function(x){log(x/(1-x))}

# Install/Load packages
package_list<-c("here", "gplots", "tidyverse", "dataRetrieval", "RColorBrewer", "openxlsx"
                , "RMark" , "R2jags", "tidybayes", "MCMCvis", "patchwork", "ggthemes", "glue", "coda")

# install_or_load_pack(package_list)
lapply(package_list, try(library), character.only = TRUE)

# Specify model files of interest
models<-c(
   "sst_model.r"
)

# Load M-R dataset via saved .rds file
location<-c("NFLewis")
species<-c("Chinook")
JS_stats_ALL<-readRDS(paste0("Data/Originals/NF_Lewis/Lewis_",year,"_JS_data.rds"))

# Function to collapse third dimension

collapse_third_dimension <- function(data_list) {
  # Initialize an empty list to store long format data frames
  long_dfs <- list()
  # Loop over each element in the list
  for (level_name in names(data_list)) {
    long_df <- data_list[[level_name]] %>%
      mutate(type = level_name) %>%
      pivot_longer(-c(dates, type), names_to = 'variable', values_to = 'value') %>%
      mutate(value = as.numeric(as.character(value)))
    long_dfs[[level_name]] <- long_df
  }
  # Combine all long data frames
  combined_long <- bind_rows(long_dfs)
  # Summarize to collapse the levels
  collapsed_df <- combined_long %>%
    group_by(dates, variable) %>%
    summarize(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = variable, values_from = value) |>
    mutate(dates = as.Date(dates))
  return(collapsed_df)
}
# Collapse 3rd dimension of JS_stats into a single grouping using "collapse_third_dimension" function
(collapsed_df <- collapse_third_dimension(JS_stats_ALL))
# collapsed_df <- (JS_stats_ALL$Male)
# for(i in 2:7)
#   collapsed_df[,i] <- as.numeric(as.character(collapsed_df[,i]))

(year<-year(min(collapsed_df$dates)))

# Define "good" JS vectors of data based on rule-sets
good_s_T_mat_1 = apply(as.matrix(collapsed_df$zi + collapsed_df$mi),2,function(x) c(which(x>0),rep(NA,length(which(x==0)))))
good_s_T_mat_2 = apply(apply(good_s_T_mat_1,2,function(x) ifelse(x%in%c(2:(nrow(collapsed_df)-1)),x,NA)),2,function(x) sort(x,na.last=TRUE))
good_s_T_final = apply(good_s_T_mat_2,2,function(x) length(x[!is.na(x)]))

good_s_R_mat_1<-apply(as.matrix(collapsed_df$Ri),2,function(x) c(which(x>0),rep(NA,length(which(x==0)))))
good_s_R_mat_2=apply(apply(good_s_R_mat_1,2,function(x) ifelse(x%in%c(1:(nrow(collapsed_df)-1)),x,NA)),2,function(x) sort(x,na.last=TRUE))
good_s_R_final<-apply(good_s_R_mat_2,2,function(x) length(x[!is.na(x)]))

# Create data object for JAGS model
MR_data<-
  list(
    a=0.5
    , b=0.5
    , alpha = 0.1
    , LL = sum(collapsed_df$ui)
    , UL = sum(collapsed_df$ui) * 20
    , s = nrow(collapsed_df)
    #, num_strata = dim(n)[2]
    , uTot = sum(collapsed_df$ui)
    , n = collapsed_df$ni
    , R = collapsed_df$Ri
    , u = collapsed_df$ui
    , r = collapsed_df$ri
    ,m = collapsed_df$mi
    , T = collapsed_df$zi + collapsed_df$mi
    , days=c(NA,as.numeric(diff(lubridate::ymd(collapsed_df$dates)))) #number of days between last period and current period...for first period, value is NA
    , good_s_T_mat = good_s_T_mat_2
    , good_s_R_mat = good_s_R_mat_2
    , good_s_T = good_s_T_final
    , good_s_R = good_s_R_final
    , logit_phi_1_mu_prior = logit(exp(log(0.5)/mean(c(NA,as.numeric(diff(lubridate::ymd(collapsed_df$dates)))),na.rm=TRUE)))
  )

if(is.null(dim(MR_data$R)[2])==TRUE){
  MR_data[lapply(MR_data,class)=="matrix"]<-lapply(MR_data[lapply(MR_data,class)=="matrix"],as.vector)
  MR_data$alpha<-rep(1,MR_data$s)
}

# Create inits for uTot
uTot<-MR_data$uTot
inits_func<-function(){
  inits<-list(
    "Ntot" = uTot * runif(1,1,5)
  )
}

# Run JS models of interest
for (i in 1:length(models)){
  start.time<-Sys.time()
  print(start.time)
  jag_out<-
    jags.parallel(
      data = MR_data
      , model.file = glue("models/{models[i]}")
      , n.chains= 4
      , n.thin= 40
      , n.burnin= 5000
      , n.iter = 20000
      , parameters.to.save=c("p",
                             "phi",
                             "temp",
                             "logit_phi",
                             "lambda",
                             "pent_tmp",
                             "psi_tmp",
                             # "pent",
                             # "log_delta",
                             # "delta",
                             "psiPtot",
                             "multP",
                             "v",
                             "tau",
                             "est_uTot",
                             "Nsuper",
                             "psi",
                             "delta"
                             # "Bstar",
                             # "rho_phi",
                             # "rho_p",
                             # "rho_pent",
                             # "sigma_phi_process",
                             # "sigma_p_process",
                             # "sigma_pent_process"
      )
      #, inits = inits_func
      # , n.chains = n_chains
      # , n.iter = n_iter
      # , n.burnin = n_burn
      # , n.thin = n_thin
      # , parameters.to.save=params_to_save
    )
}


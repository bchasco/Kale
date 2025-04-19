tmb_list <- list()
jags_list <- list()

yrs <- 2013:2024
uTot_yr <- matrix(0,2,length(yrs))

icnt <- 1
st <- Sys.time()

for(yr in yrs){
  year <- yr
  source("r/sst_JAGS.r")
  source("r/sst_RTMB_ind.r")
  jags_list[[as.character(yr)]] <- jag_out
  tmb_list[[as.character(yr)]] <- list(tmb_rep = tmb_rep, tmb_sdr = tmb_sdr)
  uTot_yr[1,icnt] <- MR_data$uTot
  uTot_yr[2,icnt] <- sum(d$n)
  icnt <- icnt  + 1
}
st <- Sys.time() - st
print(st)

library(dplyr)
library(tibble)
library(ggplot2)

summary_df <- tibble()

for (yr in yrs) {
  yr <- as.character(yr)

  # JAGS posterior samples
  jag_post <- jags_list[[yr]]$BUGSoutput$sims.list$Nsuper

  print(mean(jag_post))
  # RTMB estimate and SE
  tmb_mean <- as.list(tmb_list[[yr]]$tmb_sdr, "Estimate", report = TRUE)[["Nsuper"]]
  tmb_sd   <- as.list(tmb_list[[yr]]$tmb_sdr, "Std. Error", report = TRUE)[["Nsuper"]]

  dens <- density(jag_post)
  # Combine summaries
  summary_df <- bind_rows(
    summary_df,
    tibble(
      year = yr,
      method = "jags",
      mean = dens$x[which.max(dens$y)],
      # mean = mean(jag_post),
      lower = quantile(jag_post, 0.025),
      upper = quantile(jag_post, 0.975)
    ),
    tibble(
      year = yr,
      method = "tmb",
      mean = tmb_mean,
      lower = tmb_mean - 1.96 * tmb_sd,
      upper = tmb_mean + 1.96 * tmb_sd
    )
  )
}

# Plot
ggplot(summary_df, aes(x = year, y = mean, color = method)) +
  geom_point(position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
                position = position_dodge(width = 0.4)) +
  labs(
    y = expression(N[super]),
    x = "Year",
    title = "Estimated Nsuper by Year and Method"
  ) +
  theme_bw()


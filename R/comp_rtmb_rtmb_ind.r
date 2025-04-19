
year <- 2024

source("r/sst_JAGS.r")
source("r/sst_RTMB_ind.r")

tmb <- tmb_rep$out$v_w
jag <- colMeans(jag_out$BUGSoutput$sims.list$v)[1:length(tmb)]
obs <- (MR_data$R/MR_data$n)[1:length(tmb)]
obs_tmb <- d %>%
  group_by(t_k) %>%
  summarize(r = sum(tag*n), R = sum(n)) %>% #sum(tag * (!is.na(r_wk)))) %>%
  mutate(obs_v = r/R) %>%
  select(obs_v)

df_v <- data.frame(jag = jag,tmb = tmb, obs = obs, obs_tmb = obs_tmb$obs_v[1:length(tmb)]) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb, obs, obs_tmb), names_to = "method", values_to = "value")
v <- ggplot(df_v, aes(x = id, y = value, color = method)) +
  geom_line(data = filter(df_v, !(method %in% c("obs","obs_tmb"))), size = 1.2, alpha = 0.7) +
  # geom_point(data = filter(df_v, method != "obs"), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_v, method %in% c("obs","obs_tmb")), size = 3) +
  labs(x = "Strata", y = "Marking rate") +
  theme_bw()
print(v)

tmb <- tmb_rep$out$tau
jag <- colMeans(jag_out$BUGSoutput$sims.list$tau)[1:length(tmb)]
obs <- (MR_data$m/MR_data$T)[1:length(tmb)]
# Filter for tagged fish only
tagged <- d %>%
  filter(tag == TRUE)

m <- d %>%
  filter(!is.na(r_k)) %>%
  group_by(t_k) %>%
  summarize(m = sum(n), .groups = "drop") %>%
  complete(t_k = 1:s, fill = list(m = 0)) %>%
  arrange(t_k)

# Expand each group to all available periods between t_k + 1 and r_k (or s if not recaptured)
expanded_T <- tagged %>%
  mutate(
    t_start = t_k,
    t_end = ifelse(is.na(r_k), s, r_k)
  ) %>%
  rowwise() %>%
  mutate(
    available_times = list(seq(t_start, t_end))
  ) %>%
  unnest(available_times) %>%
  filter(available_times <= s) %>%
  group_by(available_times) %>%
  summarize(T = sum(n), .groups = "drop") %>%
  complete(available_times = 1:s, fill = list(T = 0)) %>%
  arrange(available_times)

# Rename for clarity
m$T_vec <- c(expanded_T$T)
obs_tmb <- m$m/m$T_vec;

df_tau <- data.frame(jag = jag, tmb = tmb, obs = c(obs)) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb, obs), names_to = "method", values_to = "value")
tau <- ggplot(df_tau, aes(x = id, y = value, color = method)) +
  geom_line(data = filter(df_tau, !(method %in% c("obs","obs_tmb"))), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_tau, (method %in% c("obs","obs_tmb"))), size = 3) +
  labs(x = "Strata", y = "Recapture probability (tau)") +
  theme_bw()
print(tau)


tmb <- tmb_rep$out$lambda
jag <- colMeans(jag_out$BUGSoutput$sims.list$lambda)[1:length(tmb)]
obs <- (MR_data$r/MR_data$R)[1:length(tmb)]
df_lambda <- data.frame(jag = jag, tmb = tmb, obs = obs) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb, obs), names_to = "method", values_to = "value")
lambda <- ggplot(df_lambda, aes(x = id, y = value, color = method)) +
  geom_line(data = filter(df_lambda, method != "obs"), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_lambda, method == "obs"), size = 3) +
  labs(x = "Strata", y = "Survival + recruitment (Lambda)") +
  theme_bw()
print(lambda)


tmb <- tmb_rep$out$multP
jag <- colMeans(jag_out$BUGSoutput$sims.list$multP)[1:length(tmb)]
obs <- (MR_data$u/MR_data$uTot)[1:length(tmb)]
df_multP <- data.frame(jag = jag, tmb = tmb, obs = obs) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb, obs), names_to = "method", values_to = "value")
multP <- ggplot(df_multP, (aes(x = id, y = value, color = method))) +
  geom_line(data = filter(df_multP, method != "obs"), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_multP, method == "obs"), size = 3) +
  ylab("multP") +
  xlab("Strata") +
  theme_bw()
print(multP)

tmb <- tmb_rep$out$delta_w
tmb <- tmb/sum(tmb)
jag <- colMeans(jag_out$BUGSoutput$sims.list$delta)[1:length(tmb)]
jag <- jag/sum(jag)
df_delta <- data.frame(jag = jag, tmb = tmb) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb), names_to = "method", values_to = "value")
delta <- ggplot(df_delta, aes(x = id, y = value, color = method)) +
  geom_line(data = filter(df_delta, method != "obs"), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_delta, method == "obs"), size = 3) +
  labs(x = "Strata", y = "New entrants (pent = f(delta)") +
  theme_bw()
print(delta)

tmb <- tmb_rep$out$temp
# tmb <- tmb/sum(tmb)
jag <- colMeans(jag_out$BUGSoutput$sims.list$temp)[1:length(tmb)]
# jag <- jag/sum(jag)
df_temp <- data.frame(jag = jag, tmb = tmb) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb), names_to = "method", values_to = "value")
temp <- ggplot(df_temp, aes(x = id, y = value, color = method)) +
  geom_line(data = filter(df_temp, method != "obs"), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_temp, method == "obs"), size = 3) +
  labs(x = "Strata", y = "New entrants (pent = f(temp)") +
  theme_bw()
print(temp)

tmb <- tmb_rep$out$psi
jag <- colMeans(jag_out$BUGSoutput$sims.list$psi)[1:length(tmb)]
obs <- (MR_data$u/MR_data$uTot)[1:length(tmb)]
df_multP <- data.frame(jag = jag, tmb = tmb) %>%
  rowid_to_column("id") %>%
  pivot_longer(cols = c(jag, tmb), names_to = "method", values_to = "value")
psi <- ggplot(df_multP, (aes(x = id, y = value, color = method))) +
  geom_line(data = filter(df_multP, method != "obs"), size = 1.2, alpha = 0.7) +
  geom_point(data = filter(df_multP, method == "obs"), size = 3) +
  ylab("Probability an individual is present (psi = f(phi,p,pent))") +
  xlab("Strata") +
  theme_bw()
print(psi)

# tmb <- tmb_rep$out$psi_tmp
# jag <- colMeans(jag_out$BUGSoutput$sims.list$psi_tmp)[1:length(tmb)]
# obs <- (MR_data$u/MR_data$uTot)[1:length(tmb)]
# df_pent_temp <- data.frame(jag = jag, tmb = tmb) %>%
#   rowid_to_column("id") %>%
#   pivot_longer(cols = c(jag, tmb), names_to = "method", values_to = "value")
# pent_tmp <- ggplot(df_pent_temp, (aes(x = id, y = value, color = method))) +
#   geom_line(data = filter(df_pent_temp, method != "obs"), size = 1.2, alpha = 0.7) +
#   geom_point(data = filter(df_pent_temp, method == "obs"), size = 3) +
#   ylab("Pent_temp") +
#   xlab("Strata") +
#   theme_bw()
# print(pent_tmp)


ndraw <- length(jag_out$BUGSoutput$sims.list$phi)
jag <- list(phi = c(jag_out$BUGSoutput$sims.list$phi),
                    p = c(jag_out$BUGSoutput$sims.list$p),
                    Nsuper = c(jag_out$BUGSoutput$sims.list$Nsuper),
                    psiPtot = c(jag_out$BUGSoutput$sims.list$psiPtot)) %>%
  reshape2::melt() %>%
  mutate(method = "jag")


g_par <- list(phi = exp(rnorm(ndraw, sdr_est$log_phi - 0.5*sdr_se$log_phi^2, sdr_se$log_phi)),
                     p = exp(rnorm(ndraw, sdr_est$log_p - 0.5*sdr_se$log_p^2, sdr_se$log_p)),
                     Nsuper = rnorm(ndraw, sdr_est$Nsuper, sdr_se$Nsuper),
                     psiPtot = exp(rnorm(ndraw, sdr_est$log_psiPtot - 0.5*sdr_se$log_phi^2, sdr_se$log_psiPtot))) %>%
  reshape2::melt() %>%
  mutate(method = "tmb") %>%
  rbind(jag) %>%
  ggplot(aes(x = value, color = method, fill = method)) +
  facet_wrap( ~ L1, scales = "free") +
  geom_histogram(position = "identity", alpha = 0.3) +
  theme_bw()

print(g_par)


g <- ggpubr::ggarrange(multP,tau,lambda,v, g_par) +
  theme(text = element_text(size = 14))
print(g)


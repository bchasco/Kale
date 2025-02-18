pvd <- read.csv("data/simpleData2.csv") %>%
  mutate(t_wk = lubridate::week(lubridate::mdy(TagDate)),
         r_wk = lubridate::week(lubridate::mdy(RecapDate))) %>%
  filter(t_wk > 10) %>%
  mutate(t_k = t_wk - min(t_wk) + 1,
         r_k = r_wk - min(t_wk) + 1) %>%
  mutate(t_l = TagState,
         r_l = RecapState) %>%
  filter(!(Tag1=="")) %>%
  count(t_l, r_l) %>%
  pivot_wider(names_from = r_l, values_from = n, values_fill = 0)


pobs <- round(pvd[,2:7]/rowSums(pvd),2)

est <- as.list(sd,"Estimate",report=TRUE) %>%
  reshape2::melt()
se <- as.list(sd,"Std. Error",report=TRUE) %>%
  reshape2::melt()

g <- data.frame(est,se = se$value, obs = c(unlist(pobs))) %>%
  filter(plogis(value)>0) %>%
  ggplot(aes(x=Var2, y = plogis(value))) +
  geom_col(aes(x=Var2,y = obs), color = "lightgrey", alpha = 0.2) +
  geom_point(color = "red", alpha = 0.5, size = 2) +
  geom_errorbar(aes(ymin = plogis(value - 1.96 * se), ymax = plogis(value + 1.96 * se)), width = 0.2, color = "red") +
  facet_wrap(~Var1, ncol = 2) +
  theme_classic() +
  theme(text = element_text(size = 16))
print(g)
ggsave(filename = "output/plot_the_data.png", g, height = 8, width = 8)

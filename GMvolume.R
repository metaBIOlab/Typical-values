## AvgThickness ---
m.1 <-
  lme4::lmer(GMvolume ~ Age * ROI + (1|Sample:ROI), data = dados_datasetscomp_rate)

re <- as_tibble(ranef(m.1)) %>%
  filter(grpvar == "Sample:ROI") %>%
  mutate(
    T_shift = condval,
    sd_T_shift = condsd,
    Sample = str_split(grp, pattern = ":", simplify = TRUE)[, 1],
    ROI = str_split(grp, pattern = ":", simplify = TRUE)[, 2]
  ) %>%
  dplyr::select(-c(condval, grpvar, term, condsd, grp))

Age.trend <- as_tibble(lstrends(m.1, ~ ROI, var = "Age")) %>%
  mutate(Age.trend_T = Age.trend) %>%
  dplyr::select(c(ROI, Age.trend_T))

dados_datasetscomp <- full_join(dados_datasetscomp, Age.trend) %>%
  full_join(re) %>%
  mutate(
    AvgThickness_shiftc = AvgThickness - T_shift,
    AvgThickness_age_decay = AvgThickness - Age.trend_T * (Age - Age.cor),
    AvgThickness_age_decay_shiftc = AvgThickness - T_shift - Age.trend_T *
      (Age - Age.cor),
    logAvgThickness_shiftc = log10(AvgThickness_shiftc),
    logAvgThickness_age_decay = log10(AvgThickness_age_decay),
    logAvgThickness_age_decay_shiftc = log10(AvgThickness_age_decay_shiftc)
  )


ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"), aes(x = Age, y = GMvolume, color = Sample, alpha = 0.4))+
  geom_point() +
  theme_pubr() + 
  guides(alpha = "none") +
  scale_x_continuous(limits = c(0,100))

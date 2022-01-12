# DEAGING + HARMONIZATION FULL SAMPLE ---- , Age > 18, Sample == "NKI"

dados_datasetscomp_rate <-
  filter(dados_datasetscomp, Diagnostic == "CTL")

#dados_datasetscomp_rate$Diagnostic <- factor(dados_datasetscomp_rate$Diagnostic, levels = c("CTL", "MCI","AD"))
dados_datasetscomp_rate$Sample <-
  as.factor(dados_datasetscomp_rate$Sample)
dados_datasetscomp_rate$ROI <-
  factor(dados_datasetscomp_rate$ROI,
         levels = c("hemisphere", "F", "O", "P", "T"))

## AvgThickness ---
m.1 <-
  lme4::lmer(AvgThickness ~ Age * ROI + (1|Sample:ROI), data = dados_datasetscomp_rate)

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

## TotalArea ---
m.1 <-
  lme4::lmer(TotalArea_corrected ~ Age * ROI + (1 | Sample:ROI), data = dados_datasetscomp_rate)

re <- as_tibble(ranef(m.1)) %>%
  filter(grpvar == "Sample:ROI") %>%
  mutate(
    AT_shift = condval,
    sd_AT_shift = condsd,
    Sample = str_split(grp, pattern = ":", simplify = TRUE)[, 1],
    ROI = str_split(grp, pattern = ":", simplify = TRUE)[, 2]
  ) %>%
  dplyr::select(-c(condval, grpvar, term, condsd, grp))

Age.trend <-
  as_tibble(lstrends(m.1, ~ ROI, var = "Age")) %>%
  mutate(Age.trend_AT = Age.trend) %>%
  dplyr::select(c(ROI, Age.trend_AT))

dados_datasetscomp <- full_join(dados_datasetscomp, re) %>%
  full_join(Age.trend) %>%
  mutate(
    TotalArea_shiftc = TotalArea_corrected - AT_shift,
    TotalArea_age_decay = TotalArea_corrected - Age.trend_AT * (Age - Age.cor),
    TotalArea_age_decay_shiftc = TotalArea_corrected - AT_shift - Age.trend_AT *
      (Age - Age.cor),
    logTotalArea_shiftc = log10(TotalArea_shiftc),
    logTotalArea_age_decay = log10(TotalArea_age_decay),
    logTotalArea_age_decay_shiftc = log10(TotalArea_age_decay_shiftc)
  )

## ExposedArea ---
m.1 <-
  lme4::lmer(
    ExposedArea_corrected  ~ Age * ROI + (1 | Sample:ROI), data = dados_datasetscomp_rate)

re <- as_tibble(ranef(m.1)) %>%
  filter(grpvar == "Sample:ROI") %>%
  mutate(
    AE_shift = condval,
    sd_AE_shift = condsd,
    Sample = str_split(grp, pattern = ":", simplify = TRUE)[, 1],
    ROI = str_split(grp, pattern = ":", simplify = TRUE)[, 2]
  ) %>%
  dplyr::select(-c(condval, grpvar, term, condsd, grp))

Age.trend <-
  as_tibble(lstrends(m.1, ~ ROI, var = "Age")) %>%
  mutate(Age.trend_AE = Age.trend) %>%
  dplyr::select(c(ROI, Age.trend_AE))

dados_datasetscomp <- full_join(dados_datasetscomp, re) %>%
  full_join(Age.trend) %>%
  mutate(
    ExposedArea_shiftc = ExposedArea_corrected - AE_shift,
    ExposedArea_age_decay = ExposedArea_corrected - Age.trend_AE * (Age - Age.cor),
    ExposedArea_age_decay_shiftc = ExposedArea_corrected - AE_shift - Age.trend_AE * (Age - Age.cor),
    logExposedArea_shiftc = log10(ExposedArea_shiftc),
    logExposedArea_age_decay = log10(ExposedArea_age_decay),
    logExposedArea_age_decay_shiftc = log10(ExposedArea_age_decay_shiftc)
  )

## AvgThickness ---
m.1 <-
  lme4::lmer(c ~ Age * ROI + (1|Sample:ROI), data = dados_datasetscomp_rate)

re <- as_tibble(ranef(m.1)) %>%
  filter(grpvar == "Sample:ROI") %>%
  mutate(
    c_shift = condval,
    sd_c_shift = condsd,
    Sample = str_split(grp, pattern = ":", simplify = TRUE)[, 1],
    ROI = str_split(grp, pattern = ":", simplify = TRUE)[, 2]
  ) %>%
  dplyr::select(-c(condval, grpvar, term, condsd, grp))

Age.trend <- as_tibble(lstrends(m.1, ~ ROI, var = "Age")) %>%
  mutate(Age.trend_c = Age.trend) %>%
  dplyr::select(c(ROI, Age.trend_c))

dados_datasetscomp <- full_join(dados_datasetscomp, Age.trend) %>%
  full_join(re) %>%
  mutate(
    c_shiftc = c - c_shift,
    c_age_decay = c - Age.trend_c * (Age - Age.cor),
    c_age_decay_shiftc = c - c_shift - Age.trend_c *
      (Age - Age.cor)
  )

# ----
dados_datasetscomp <- dados_datasetscomp %>%
  mutate(
    localGI_age_decay = TotalArea_age_decay/ExposedArea_age_decay,
    localGI_shiftc = TotalArea_shiftc/ExposedArea_shiftc,
    localGI_age_decay_shiftc = TotalArea_age_decay_shiftc/ExposedArea_age_decay_shiftc,
    K_age_decay = log10(TotalArea_age_decay) + 1/4*log10(AvgThickness_age_decay^2) - 5/4*log10(ExposedArea_age_decay),
    K_shiftc = log10(TotalArea_shiftc) + 1/4*log10(AvgThickness_shiftc^2) - 5/4*log10(ExposedArea_shiftc),
    K_age_decay_shiftc = log10(TotalArea_age_decay_shiftc) + 1/4*log10(AvgThickness_age_decay_shiftc^2) - 5/4*log10(ExposedArea_age_decay_shiftc),
    I_age_decay = log10(TotalArea_age_decay) + log10(ExposedArea_age_decay) + log10(AvgThickness_age_decay^2),
    I_shiftc = log10(TotalArea_shiftc) + log10(ExposedArea_shiftc) + log10(AvgThickness_shiftc^2),
    I_age_decay_shiftc = log10(TotalArea_age_decay_shiftc) + log10(ExposedArea_age_decay_shiftc) + log10(AvgThickness_age_decay_shiftc^2),
    S_age_decay = 3/2*log10(TotalArea_age_decay) + 3/4*log10(ExposedArea_age_decay) - 9/4*log10(AvgThickness_age_decay^2),
    S_shiftc = 3/2*log10(TotalArea_shiftc) + 3/4*log10(ExposedArea_shiftc) - 9/4*log10(AvgThickness_shiftc^2),
    S_age_decay_shiftc = 3/2*log10(TotalArea_age_decay_shiftc) + 3/4*log10(ExposedArea_age_decay_shiftc) - 9/4*log10(AvgThickness_age_decay_shiftc^2),
    Knorm_shiftc = K_shiftc / sqrt(1 + (1 / 4) ^ 2 + (5 / 2) ^ 2),
    Snorm_shiftc = S_shiftc / sqrt((3 / 2) ^ 2 + (3 / 4) ^ 2 + (9 / 4) ^ 2),
    Inorm_shiftc = I_shiftc / sqrt(1 ^ 2 + 1 ^ 2 + 1 ^ 2),
    Knorm_age_decay = K_age_decay / sqrt(1 + (1 / 4) ^ 2 + (5 / 2) ^ 2),
    Snorm_age_decay = S_age_decay / sqrt((3 / 2) ^ 2 + (3 / 4) ^ 2 + (9 / 4) ^ 2),
    Inorm_age_decay = I_age_decay / sqrt(1 ^ 2 + 1 ^ 2 + 1 ^ 2),
    Knorm_age_decay_shiftc = K_age_decay_shiftc / sqrt(1 + (1 / 4) ^ 2 + (5 / 2) ^ 2),
    Snorm_age_decay_shiftc = S_age_decay_shiftc / sqrt((3 / 2) ^ 2 +  (3 / 4) ^ 2 + (9 / 4) ^ 2),
    Inorm_age_decay_shiftc = I_age_decay_shiftc / sqrt(1 ^ 2 + 1 ^ 2 + 1 ^ 2)
  )


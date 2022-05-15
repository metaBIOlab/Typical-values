---
title: "Manuscript 2 from Fernanda Hansen P. de Moraes Thesis - Cortical folding alterations in humans due to aging and diseases"
author: "Fernanda Hansen Pacheco de Moraes"
date: "18 jan 2022"
output:
  html_document: 
    fig_caption: yes
    fig_width: 8
    number_sections: yes
    theme: paper
    toc: yes
editor_options: 
  chunk_output_type: inline
---

Description of the procedures and analysis present in Manuscript 2,
**Establishing a baseline for human cortical folding morphological variables: a multicenter study**, at the Doctorate Thesis presented
to the Programa de Pós-Graduação em Ciências Médicas at the Instituto
D'Or de Pesquisa e Ensino as a partial requirement to obtain the
Doctorate Degree.

Part of the data used here cannot be shared due to restrictions of the
Ethic Committee. Data can be shared upon reasonable request to the
corresponding author. To fulfill these limitation, we will generate
random data to simulate the results.

Get in touch with us
([fernandahmoraes\@gmail.com](mailto:fernandahmoraes@gmail.com){.email})
in case any help is needed, our aim is to improve the code as needed!

# RMARKDOWN AND R SETUP



```r
setwd("D:/GitHub/Typical-values")
```


```r
## define functions

# test angular coeficinet versus theoretical value
test_coef <- function(reg, coefnum, val){
  co <- coef(summary(reg))
  tstat <- (co[coefnum,1] - val)/co[coefnum,2]
  2 * pt(abs(tstat), reg$df.residual, lower.tail = FALSE)
}

# wrap text
wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")
```


```r
library(readr)
library(tidyverse)
library(lubridate)
library(ggpubr)
library(kableExtra)
library(broom)
library(lme4)
library(lsmeans)
library(MuMIn)
library(arm)
library(effects)
library(compute.es)
```


```r
# COLOR BLIND PALETTE WITH BLACK
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette2 <- c("#56B4E9", "#E69F00", "#CC79A7", "#56B4E9", "#D55E00", "#0072B2", "#009E73", "#F0E442")
```

# set seed for random process


```r
set.seed(1)
```


```r
# dados_datasetscomp <- read_csv("dados_datasetscomp2.csv")
dados_datasetscomp <- read_csv("data_typical_values.csv")

dados_datasetscomp <- dados_datasetscomp %>%
  filter(Sample != "HCP500r", Diagnostic != "MCI")
```

# DATA PREPARATION

```r
# estimate cortical folding variables
dados_datasetscomp <- dados_datasetscomp %>%
  mutate(
    # create new variables
    GMvolume = ifelse(!is.na(GMvolume),GMvolume,AvgThickness*TotalArea),
    logAvgThickness = log10(AvgThickness),
    logTotalArea = log10(TotalArea),
    logExposedArea = log10(ExposedArea),
    localGI = TotalArea / ExposedArea,
    k = sqrt(AvgThickness) * TotalArea / (ExposedArea ^ 1.25),
    K = 1 / 4 * log10(AvgThickness ^ 2)  + log10(TotalArea) - 5 / 4 * log10(ExposedArea),
    S = 3 / 2 * log10(TotalArea) + 3 / 4 * log10(ExposedArea) - 9 /  4 * log10(AvgThickness ^
                                                                                 2) ,
    I = log10(TotalArea) + log10(ExposedArea) + log10(AvgThickness ^ 2),
    # c = as.double(ifelse(
    #   ROI == "hemisphere", NA, 4 * pi / GaussianCurvature
    # )),
    TotalArea_corrected = ifelse(ROI == "hemisphere", TotalArea, TotalArea * c),
    ExposedArea_corrected = ifelse(ROI == "hemisphere", ExposedArea, ExposedArea * c),
    logTotalArea_corrected = log10(TotalArea_corrected),
    logExposedArea_corrected = log10(ExposedArea_corrected),
    localGI_corrected = ifelse(
      ROI == "hemisphere",
      TotalArea / ExposedArea,
      TotalArea_corrected / ExposedArea_corrected
    ),
    k_corrected = ifelse(
      ROI == "hemisphere",
      sqrt(AvgThickness) * log10(TotalArea) / (log10(ExposedArea) ^ 1.25),
      sqrt(AvgThickness) * log10(TotalArea_corrected) / (log10(ExposedArea_corrected ^
                                                                 1.25))
    ),
    K_corrected =  ifelse(
      ROI == "hemisphere",
      1 / 4 * log10(AvgThickness ^ 2) + log10(TotalArea) - 5 / 4 * log10(ExposedArea),
      1 / 4 * log10(AvgThickness ^ 2) + log10(TotalArea_corrected) - 5 / 4 * log10(ExposedArea_corrected)
    ),
    I_corrected = ifelse(
      ROI == "hemisphere",
      log10(TotalArea) + log10(ExposedArea) + log10(AvgThickness ^ 2) ,
      log10(TotalArea_corrected) + log10(ExposedArea_corrected) + log10(AvgThickness ^ 2)
    ),
    S_corrected = ifelse(
      ROI == "hemisphere",
      3 / 2 * log10(TotalArea) + 3 / 4 * log10(ExposedArea) - 9 /  4 * log10(AvgThickness ^ 2) ,
      3 / 2 * log10(TotalArea_corrected) + 3 / 4 * log10(ExposedArea_corrected) - 9 /  4 * log10(AvgThickness ^ 2)
    ),
    Knorm = K_corrected / sqrt(1 + (1 / 4) ^ 2 + (5 / 4) ^ 2),
    Snorm = S_corrected / sqrt((3 / 2) ^ 2 + (3 / 4) ^ 2 + (9 / 4) ^ 2),
    Inorm = I_corrected / sqrt(1 ^ 2 + 1 ^ 2 + 1 ^ 1)
  )

# create age intervals
dados_datasetscomp$Age_interval <- cut(dados_datasetscomp$Age,
                                       breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100),
                                       right = FALSE,
                                       include.lowest = TRUE)

dados_datasetscomp$Age_interval10 <- cut(dados_datasetscomp$Age,
                                         breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                                         right = FALSE,
                                         include.lowest = TRUE)
```


```r
dados_all <- dados_datasetscomp %>% filter(!is.na(logAvgThickness), ExposedArea != 0 | !is.na(localGI), !is.infinite(logExposedArea)) %>% 
  droplevels()

dados_datasetscomp <- dados_all

dados_datasetscomp$Diagnostic <- as.factor(dados_datasetscomp$Diagnostic)
dados_datasetscomp$Diagnostic <- relevel(dados_datasetscomp$Diagnostic, ref = "CTL")
```

# Deaging


```r
# define age for deaging
Age.cor = 25

# DEAGING + HARMONIZATION FULL DATA ----

dados_datasetscomp_rate <-
  filter(dados_datasetscomp, Diagnostic == "CTL")

dados_datasetscomp_rate$Sample <-
  as.factor(dados_datasetscomp_rate$Sample)

dados_datasetscomp_rate$ROI <-
  factor(dados_datasetscomp_rate$ROI,
         levels = c("hemisphere", "F", "O", "P", "T"))
```

## GM volume


```r
m.1 <-
  lme4::lmer(GMvolume ~ Age * ROI + (1|Sample:ROI), data = dados_datasetscomp_rate)

re <- as_tibble(ranef(m.1)) %>%
  filter(grpvar == "Sample:ROI") %>%
  mutate(
    GM_shift = condval,
    sd_GM_shift = condsd,
    Sample = str_split(grp, pattern = ":", simplify = TRUE)[, 1],
    ROI = str_split(grp, pattern = ":", simplify = TRUE)[, 2]
  ) %>%
  dplyr::select(-c(condval, grpvar, term, condsd, grp))

Age.trend <- as_tibble(lstrends(m.1, ~ ROI, var = "Age")) %>%
  mutate(Age.trend_GM = Age.trend) %>%
  dplyr::select(c(ROI, Age.trend_GM))

dados_datasetscomp <- full_join(dados_datasetscomp, Age.trend) %>%
  full_join(re) %>%
  mutate(
    GMvolume_shiftc = GMvolume - GM_shift,
    GMvolume_age_decay = GMvolume - Age.trend_GM * (Age - Age.cor),
    GMvolume_age_decay_shiftc = GMvolume - GM_shift - Age.trend_GM *
      (Age - Age.cor)
  )
```

## AvgThickness


```r
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
```

## TotalArea


```r
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
```

## ExposedArea


```r
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
```

## Lobes correction factor


```r
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
```


```r
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
```

# RESULTS

## Data description
**Table 1**
<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:left;"> Diagnostic </th>
   <th style="text-align:right;"> N </th>
   <th style="text-align:left;"> age </th>
   <th style="text-align:left;"> age_range </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ADNI </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 868 </td>
   <td style="text-align:left;"> 75 ± 6.5 </td>
   <td style="text-align:left;"> 56 ;  96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ADNI </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:right;"> 542 </td>
   <td style="text-align:left;"> 75 ± 8 </td>
   <td style="text-align:left;"> 56 ;  92 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHEAD </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:left;"> 42 ± 19 </td>
   <td style="text-align:left;"> 24 ;  76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AOMICPIOP1 </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 208 </td>
   <td style="text-align:left;"> 22 ± 1.8 </td>
   <td style="text-align:left;"> 18 ;  26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AOMICPIOP2 </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 224 </td>
   <td style="text-align:left;"> 22 ± 1.8 </td>
   <td style="text-align:left;"> 18 ;  26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCPr900 </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 881 </td>
   <td style="text-align:left;"> 29 ± 3.6 </td>
   <td style="text-align:left;"> 24 ;  37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IDOR </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 77 </td>
   <td style="text-align:left;"> 66 ± 8.4 </td>
   <td style="text-align:left;"> 43 ;  80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IDOR </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:left;"> 77 ± 6.1 </td>
   <td style="text-align:left;"> 63 ;  86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IXI-Guys </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 314 </td>
   <td style="text-align:left;"> 51 ± 16 </td>
   <td style="text-align:left;"> 20 ;  86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IXI-HH </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 181 </td>
   <td style="text-align:left;"> 47 ± 17 </td>
   <td style="text-align:left;"> 20 ;  82 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IXI-IOP </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:left;"> 42 ± 17 </td>
   <td style="text-align:left;"> 20 ;  86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NKI </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:left;"> 34 ± 19 </td>
   <td style="text-align:left;"> 4 ;  85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OASIS </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 312 </td>
   <td style="text-align:left;"> 45 ± 24 </td>
   <td style="text-align:left;"> 18 ;  94 </td>
  </tr>
</tbody>
</table>

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Diagnostic </th>
   <th style="text-align:right;"> N </th>
   <th style="text-align:left;"> age </th>
   <th style="text-align:left;"> age_range </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 3095 </td>
   <td style="text-align:left;"> 46 ± 23 </td>
   <td style="text-align:left;"> 4 ;  96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:right;"> 555 </td>
   <td style="text-align:left;"> 75 ± 8 </td>
   <td style="text-align:left;"> 56 ;  92 </td>
  </tr>
</tbody>
</table>

## Correlation between K, S, and I with Age for the Healthy Control group

### K 


```r
figS6a <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"),
       aes(x = Age , y = K)) +
  geom_point(aes(color = Sample, alpha = 0.4)) +
  geom_smooth(color = "black", method = "lm") +
  # stat_cor() +
  theme_pubr() +
  guides(alpha = "none")+ 
  labs(x = "Age [years]") +
  scale_x_continuous(limits = c(0,100))

cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$K, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$K and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -100.83, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.7834543 -0.7643992
## sample estimates:
##        cor 
## -0.7741021
```

### S


```r
figS6b <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"),
       aes(x = Age, y = S)) +
  geom_point(aes(color = Sample, alpha = 0.4)) +
  geom_smooth(color = "black", method = "lm") +
  # stat_cor() +
  theme_pubr() +
  guides(alpha = "none", color = FALSE,fill = FALSE)+ 
  labs(x = "Age [years]") +
  theme(legend.position= "none") +
  scale_x_continuous(limits = c(0,100))

cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$S, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$S and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = 37.406, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.3931982 0.4326213
## sample estimates:
##       cor 
## 0.4131033
```

### I


```r
figS6c <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"),
       aes(x = Age, y = I)) +
  geom_point(aes(color = Sample, alpha = 0.4)) + 
  geom_smooth(color = "black", method = "lm") +
  # stat_cor() +
  theme_pubr() +
  guides(alpha = "none", color = FALSE,fill = FALSE)+ 
  labs(x = "Age [years]") +
  theme(legend.position= "none") +
  scale_x_continuous(limits = c(0,100))

cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$I, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$I and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -70.698, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.6643648 -0.6369635
## sample estimates:
##        cor 
## -0.6508761
```


```r
figS6 <- ggarrange(figS6a, ggarrange(figS6b, figS6c,labels = c("B","C"), nrow = 1, ncol = 2, font.label = list(size = 11)),labels = c("A"), nrow = 2, ncol = 1, font.label = list(size = 11), common.legend = TRUE, legend = "top")

ggsave("figS6.pdf", plot = figS6, width = 18, height = 22, units = "cm", device = "pdf")
```

**FIGURE 1**

```r
figS6
```

<div class="figure">
<img src="mainresults_review_files/figure-html/figure1-1.png" alt="\label{fig:figure1}Through age, every Healthy Control subject for the independent morphological component, K, S, and I. For HCPr900 and AHEAD, ages are determined by an interval of years, thereby considering the mean age. The solid line represents a linear regression with a 95% confidence interval." width="681.6" />
<p class="caption">\label{fig:figure1}Through age, every Healthy Control subject for the independent morphological component, K, S, and I. For HCPr900 and AHEAD, ages are determined by an interval of years, thereby considering the mean age. The solid line represents a linear regression with a 95% confidence interval.</p>
</div>

## Slope


```r
summary(lm(
  1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
  data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere",
       Diagnostic == "CTL"),
  na.action = na.omit
))
```

```
## 
## Call:
## lm(formula = 1/2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea), 
##     data = filter(dados_datasetscomp, Sample != "Mota&Houzel2015", 
##         ROI == "hemisphere", Diagnostic == "CTL"), na.action = na.omit)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.154606 -0.025469  0.000308  0.026792  0.098626 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        -1.11930    0.05546  -20.18   <2e-16 ***
## log10(ExposedArea)  1.37800    0.01208  114.07   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.03218 on 6800 degrees of freedom
## Multiple R-squared:  0.6568,	Adjusted R-squared:  0.6567 
## F-statistic: 1.301e+04 on 1 and 6800 DF,  p-value: < 2.2e-16
```

```r
summary(lm(
  1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
  data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere", Age > 18 & Age < 40,
       Diagnostic == "CTL"),
  na.action = na.omit
))
```

```
## 
## Call:
## lm(formula = 1/2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea), 
##     data = filter(dados_datasetscomp, Sample != "Mota&Houzel2015", 
##         ROI == "hemisphere", Age > 18 & Age < 40, Diagnostic == 
##             "CTL"), na.action = na.omit)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.075674 -0.012016  0.002633  0.013807  0.065150 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        -0.55119    0.05033  -10.95   <2e-16 ***
## log10(ExposedArea)  1.25912    0.01095  115.02   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.02063 on 3534 degrees of freedom
## Multiple R-squared:  0.7892,	Adjusted R-squared:  0.7891 
## F-statistic: 1.323e+04 on 1 and 3534 DF,  p-value: < 2.2e-16
```

```r
## Displays confidence interval
# tidy(lm(
#   1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
#   data = dados_hemi_v1,
#   na.action = na.omit
# ), conf.int = TRUE)

paste(
  "Student's t test comapring slope with theoretical value 1.25. t = ",
  signif(abs(coef(summary(
    lm(
      1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
      data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere",
       Diagnostic == "CTL"),
      na.action = na.omit
    )
  ))[2, 1] - 1.25) / coef(summary(
    lm(
      1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
      data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere",
       Diagnostic == "CTL"),
      na.action = na.omit
    )
  ))[2, 2], 3) 
)
```

```
## [1] "Student's t test comapring slope with theoretical value 1.25. t =  10.6"
```

```r
paste(
  "Student's t test comapring slope with theoretical value 1.25. p value = ",
  signif(test_coef(
    lm(
      1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
      data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere",
       Diagnostic == "CTL"),
      na.action = na.omit
    ),
    2,
    1.25
  ),3)
)
```

```
## [1] "Student's t test comapring slope with theoretical value 1.25. p value =  4.98e-26"
```

### Correlation 


```r
lm_Age <-
  filter(
    dados_datasetscomp,
    Sample != "Mota&Houzel2015",
    ROI == "hemisphere",
    Diagnostic == "CTL",
    Age_interval != "[0,5)",
    Age_interval != "[5,10)",
    Age_interval != "[10,15)"
  ) %>%
  group_by(Age_interval10) %>%
  do(fit_Age = tidy(
    lm(
      1 / 2 * log10(AvgThickness) +  log10(TotalArea) ~  log10(ExposedArea),
      data = .,
      na.action = na.omit
    ),
    conf.int = TRUE
  )) %>%
  unnest(cols = c(fit_Age))

lm_Age <- lm_Age %>% mutate(Age_interval10 = as.double((str_sub(Age_interval10,2,3))))

cor.test(filter(lm_Age, term == "log10(ExposedArea)")$estimate, filter(lm_Age, term == "log10(ExposedArea)")$Age_interval10)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(lm_Age, term == "log10(ExposedArea)")$estimate and filter(lm_Age, term == "log10(ExposedArea)")$Age_interval10
## t = -3.2156, df = 7, p-value = 0.01474
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.9494277 -0.2218865
## sample estimates:
##        cor 
## -0.7722149
```

## Uncertainties estimation for baseline

**Part of Table 2 - Natural fluctuation and Type B errors**

### GM volume \~ Age


```r
m.1 <- lme4::lmer(GMvolume ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: GMvolume ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: 183968.1
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.8395 -0.6727 -0.0320  0.6298  4.5351 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev.
##  Sample:Diagnostic (Intercept)  27722513  5265   
##  Sample            (Intercept) 366049150 19132   
##  Residual                      733577651 27085   
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                   Estimate Std. Error t value
## (Intercept)      303365.74    6119.08  49.577
## Age               -1061.86      27.79 -38.216
## DiagnosticAD     -67730.51    9855.71  -6.872
## Age:DiagnosticAD    725.19     105.75   6.858
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.196              
## DiagnostcAD -0.083  0.195       
## Ag:DgnstcAD  0.052 -0.263 -0.809
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.3957362 0.6067994
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##    lsmean
##     <dbl>
## 1 250218.
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 27000
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##       19000
```

### AvgThickness \~ Age

Correlation:


```r
cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$AvgThickness, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$AvgThickness and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -79.602, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.7066199 -0.6820092
## sample estimates:
##        cor 
## -0.6945177
```

Effect size of correlation: 


```r
res(
  cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$AvgThickness,filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)$estimate,
  var.r = NULL,(cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$AvgThickness,filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)$parameter + 2)/2,
  level = 95,
  dig = 2,
  verbose = TRUE
)
```

```
## Mean Differences ES: 
##  
##  d [ 95 %CI] = -1.93 [ -2.02 , -1.84 ] 
##   var(d) = 0 
##   p-value(d) = 0 
##   U3(d) = 2.68 % 
##   CLES(d) = 8.61 % 
##   Cliff's Delta = -0.83 
##  
##  Correlation ES: 
##  
##  r [ 95 %CI] = -0.69 [ -0.71 , -0.68 ] 
##   var(r) = 0 
##   p-value(r) = 0 
##  
##  z [ 95 %CI] = -0.86 [ -0.89 , -0.82 ] 
##   var(z) = 0 
##   p-value(z) = 0 
##  
##  Odds Ratio ES: 
##  
##  OR [ 95 %CI] = 0.03 [ 0.03 , 0.04 ] 
##   p-value(OR) = 0 
##  
##  Log OR [ 95 %CI] = -3.5 [ -3.67 , -3.33 ] 
##   var(lOR) = 0.01 
##   p-value(Log OR) = 0 
##  
##  Other: 
##  
##  NNT = -5.07 
##  Total N = 3401
```

R squared (correlation):


```r
glance(lm(AvgThickness ~ Age, data = filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")))$r.squared
```

```
## [1] 0.4823548
```

Linear Model - statement the harmonization is possible trhough a linear model:


```r
summary(lm(AvgThickness ~ Age, data = filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")))
```

```
## 
## Call:
## lm(formula = AvgThickness ~ Age, data = filter(dados_datasetscomp, 
##     ROI == "hemisphere", Diagnostic == "CTL"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.73402 -0.07348  0.01216  0.08971  0.36903 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  2.761e+00  3.584e-03   770.6   <2e-16 ***
## Age         -5.563e-03  6.989e-05   -79.6   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.1311 on 6800 degrees of freedom
## Multiple R-squared:  0.4824,	Adjusted R-squared:  0.4823 
## F-statistic:  6336 on 1 and 6800 DF,  p-value: < 2.2e-16
```

Model:


```r
m.1 <- lme4::lmer(AvgThickness ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## AvgThickness ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -13111.9
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -7.8038 -0.6167  0.0104  0.6495  3.4140 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev.
##  Sample:Diagnostic (Intercept) 0.0006026 0.02455 
##  Sample            (Intercept) 0.0071880 0.08478 
##  Residual                      0.0110146 0.10495 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)       2.7038037  0.0270734  99.869
## Age              -0.0044384  0.0001077 -41.194
## DiagnosticAD     -0.3421007  0.0405587  -8.435
## Age:DiagnosticAD  0.0029937  0.0004098   7.305
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.172              
## DiagnostcAD -0.084  0.183       
## Ag:DgnstcAD  0.045 -0.263 -0.762
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.4615061 0.6845916
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   2.48
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.1
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##       0.085
```

### TotalArea \~ Age

Correlation:


```r
cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$TotalArea, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$TotalArea and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -43.615, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.4859053 -0.4487599
## sample estimates:
##        cor 
## -0.4675389
```

Effect size of correlation: 


```r
res(
  cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$TotalArea,filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)$estimate,
  var.r = NULL,(cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$TotalArea,filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)$parameter + 2)/2,
  level = 95,
  dig = 2,
  verbose = TRUE
)
```

```
## Mean Differences ES: 
##  
##  d [ 95 %CI] = -1.06 [ -1.13 , -0.98 ] 
##   var(d) = 0 
##   p-value(d) = 0 
##   U3(d) = 14.51 % 
##   CLES(d) = 22.72 % 
##   Cliff's Delta = -0.55 
##  
##  Correlation ES: 
##  
##  r [ 95 %CI] = -0.47 [ -0.49 , -0.44 ] 
##   var(r) = 0 
##   p-value(r) = 0 
##  
##  z [ 95 %CI] = -0.51 [ -0.54 , -0.47 ] 
##   var(z) = 0 
##   p-value(z) = 0 
##  
##  Odds Ratio ES: 
##  
##  OR [ 95 %CI] = 0.15 [ 0.13 , 0.17 ] 
##   p-value(OR) = 0 
##  
##  Log OR [ 95 %CI] = -1.92 [ -2.06 , -1.78 ] 
##   var(lOR) = 0 
##   p-value(Log OR) = 0 
##  
##  Other: 
##  
##  NNT = -5.84 
##  Total N = 3401
```

R squared (correlation):


```r
glance(lm(TotalArea ~ Age, data = filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")))$r.squared
```

```
## [1] 0.2185926
```

Linear Model - statement the harmonization is possible trhough a linear model:


```r
summary(lm(TotalArea ~ Age, data = filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")))
```

```
## 
## Call:
## lm(formula = TotalArea ~ Age, data = filter(dados_datasetscomp, 
##     ROI == "hemisphere", Diagnostic == "CTL"))
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -37034  -7034   -374   6733  42332 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 113305.640    273.130  414.84   <2e-16 ***
## Age           -232.300      5.326  -43.62   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 9991 on 6800 degrees of freedom
## Multiple R-squared:  0.2186,	Adjusted R-squared:  0.2185 
## F-statistic:  1902 on 1 and 6800 DF,  p-value: < 2.2e-16
```

Model:


```r
m.1 <- lme4::lmer(TotalArea ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: TotalArea ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: 167642
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.9066 -0.7162 -0.0476  0.6633  4.2777 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept)        0    0    
##  Sample            (Intercept) 23918246 4891    
##  Residual                      93180143 9653    
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)      113025.061   1543.436  73.229
## Age                -244.189      9.866 -24.750
## DiagnosticAD     -13184.422   2852.379  -4.622
## Age:DiagnosticAD    152.485     37.675   4.047
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.277              
## DiagnostcAD -0.073  0.259       
## Ag:DgnstcAD  0.072 -0.263 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.2353489 0.3915348
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##    lsmean
##     <dbl>
## 1 100803.
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 9700
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##        4900
```

### ExposedArea \~ Age

Correlation:


```r
cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$ExposedArea, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$ExposedArea and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -17.297, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.2279449 -0.1824164
## sample estimates:
##        cor 
## -0.2052917
```

Effect size of correlation: 


```r
res(
  cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$ExposedArea,filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)$estimate,
  var.r = NULL,(cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$ExposedArea,filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)$parameter + 2)/2,
  level = 95,
  dig = 2,
  verbose = TRUE
)
```

```
## Mean Differences ES: 
##  
##  d [ 95 %CI] = -0.42 [ -0.49 , -0.35 ] 
##   var(d) = 0 
##   p-value(d) = 0 
##   U3(d) = 33.74 % 
##   CLES(d) = 38.34 % 
##   Cliff's Delta = -0.23 
##  
##  Correlation ES: 
##  
##  r [ 95 %CI] = -0.21 [ -0.24 , -0.17 ] 
##   var(r) = 0 
##   p-value(r) = 0 
##  
##  z [ 95 %CI] = -0.21 [ -0.24 , -0.17 ] 
##   var(z) = 0 
##   p-value(z) = 0 
##  
##  Odds Ratio ES: 
##  
##  OR [ 95 %CI] = 0.47 [ 0.41 , 0.53 ] 
##   p-value(OR) = 0 
##  
##  Log OR [ 95 %CI] = -0.76 [ -0.89 , -0.64 ] 
##   var(lOR) = 0 
##   p-value(Log OR) = 0 
##  
##  Other: 
##  
##  NNT = -10.38 
##  Total N = 3401
```

R squared (correlation):


```r
glance(lm(ExposedArea ~ Age, data = filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")))$r.squared
```

```
## [1] 0.04214469
```

Linear Model - statement the harmonization is possible trhough a linear model:


```r
summary(lm(ExposedArea ~ Age, data = filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")))
```

```
## 
## Call:
## lm(formula = ExposedArea ~ Age, data = filter(dados_datasetscomp, 
##     ROI == "hemisphere", Diagnostic == "CTL"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -8592.5 -2019.7  -161.5  1927.6 10994.8 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) 40342.093     78.049   516.9   <2e-16 ***
## Age           -26.326      1.522   -17.3   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 2855 on 6800 degrees of freedom
## Multiple R-squared:  0.04214,	Adjusted R-squared:  0.042 
## F-statistic: 299.2 on 1 and 6800 DF,  p-value: < 2.2e-16
```

Model:


```r
m.1 <- lme4::lmer(ExposedArea ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## ExposedArea ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: 148576.9
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.8868 -0.7093 -0.0618  0.6706  3.8432 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept)       0     0.0  
##  Sample            (Intercept)  279141   528.3  
##  Residual                      8382585  2895.3  
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                   Estimate Std. Error t value
## (Intercept)      40714.770    207.022 196.669
## Age                -40.602      2.859 -14.204
## DiagnosticAD      -262.931    853.269  -0.308
## Age:DiagnosticAD     1.535     11.275   0.136
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.597              
## DiagnostcAD -0.151  0.249       
## Ag:DgnstcAD  0.152 -0.255 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
r.squaredGLMM(m.1)
```

```
##             R2m       R2c
## [1,] 0.09982676 0.1288366
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1 38683.
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 2900
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##         530
```

### GI \~ Age


```r
m.1 <- lme4::lmer(localGI ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: localGI ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -14689.7
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.7049 -0.6845 -0.0074  0.6784  6.4692 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept) 0.000000 0.00000 
##  Sample            (Intercept) 0.012446 0.11156 
##  Residual                      0.009017 0.09496 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)       2.7763677  0.0339367   81.81
## Age              -0.0034689  0.0000975  -35.58
## DiagnosticAD     -0.3178661  0.0280689  -11.32
## Age:DiagnosticAD  0.0037476  0.0003707   10.11
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.124              
## DiagnostcAD -0.033  0.260       
## Ag:DgnstcAD  0.033 -0.264 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.2623413 0.6900971
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   2.60
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.095
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##        0.11
```

### K \~ Age


```r
m.1 <- lme4::lmer(K ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: K ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -42555
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.9733 -0.6065  0.0381  0.6446  5.2725 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev.
##  Sample:Diagnostic (Intercept) 4.607e-06 0.002146
##  Sample            (Intercept) 4.410e-04 0.021001
##  Residual                      2.658e-04 0.016304
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)      -4.920e-01  6.412e-03  -76.73
## Age              -8.612e-04  1.675e-05  -51.40
## DiagnosticAD     -8.439e-02  5.444e-03  -15.50
## Age:DiagnosticAD  8.855e-04  6.366e-05   13.91
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.113              
## DiagnostcAD -0.039  0.221       
## Ag:DgnstcAD  0.030 -0.264 -0.882
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.4372324 0.7897408
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1 -0.535
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.016
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##       0.021
```

### S \~ Age


```r
m.1 <- lme4::lmer(S ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: S ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -10789.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.5837 -0.6903 -0.0315  0.6406  5.5519 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev.
##  Sample:Diagnostic (Intercept) 0.0003168 0.01780 
##  Sample            (Intercept) 0.0057953 0.07613 
##  Residual                      0.0147855 0.12160 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)       9.0840260  0.0242645 374.375
## Age               0.0017030  0.0001246  13.664
## DiagnosticAD      0.1969053  0.0413152   4.766
## Age:DiagnosticAD -0.0013879  0.0004747  -2.924
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.222              
## DiagnostcAD -0.079  0.212       
## Ag:DgnstcAD  0.058 -0.263 -0.866
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m      R2c
## [1,] 0.1515722 0.399721
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   9.17
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.12
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##       0.076
```

### I \~ Age


```r
m.1 <- lme4::lmer(I ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))
```

Model summary and R squared:

```r
summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: I ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -17050
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.9306 -0.6593 -0.0005  0.6762  3.2113 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev.
##  Sample:Diagnostic (Intercept) 0.0001859 0.01364 
##  Sample            (Intercept) 0.0013945 0.03734 
##  Residual                      0.0067038 0.08188 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)       1.053e+01  1.259e-02 836.069
## Age              -3.118e-03  8.371e-05 -37.254
## DiagnosticAD     -1.829e-01  2.848e-02  -6.421
## Age:DiagnosticAD  1.747e-03  3.196e-04   5.467
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.287              
## DiagnostcAD -0.109  0.200       
## Ag:DgnstcAD  0.075 -0.262 -0.844
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.4505911 0.5554086
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   10.4
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.082
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##       0.037
```

## Multisite harmonization

<img src="mainresults_review_files/figure-html/unnamed-chunk-68-1.png" width="768" />

<img src="mainresults_review_files/figure-html/unnamed-chunk-69-1.png" width="768" />


```r
S_lifespan_b <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"), aes(x = Age, y = S_shiftc, alpha = 0.4))+
    geom_point(aes(color = Sample, fill = Sample, alpha = 0.1)) +
    geom_smooth(method = "lm", color = "black") +
  geom_smooth(aes(color = Sample, alpha = 0.1), method = "lm") +
    guides(alpha = "none") + 
  theme_pubr() +
  labs(y = "S (harmonized)") +
  scale_x_continuous(limits = c(0,100))

I_lifespan_b <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"), aes(x = Age, y = I_shiftc, alpha = 0.4))+
    geom_point(aes(color = Sample, fill = Sample, alpha = 0.1)) +
    geom_smooth(method = "lm", color = "black") +
    geom_smooth(aes(color = Sample, alpha = 0.1), method = "lm") +
    guides(alpha = "none") + 
  theme_pubr() +
  labs(y = "I (harmonized)") +
  scale_x_continuous(limits = c(0,100))

lifespan_shift <- ggarrange(K_lifespan_b, S_lifespan_b, I_lifespan_b,labels = c("A","B", "C"), nrow = 3, ncol = 1, font.label = list(size = 11), common.legend = TRUE, legend = "right", vjust = 1.2, align = "hv")

ggsave("lifespan_shift.pdf", plot = lifespan_shift, width = 18, height = 16, units = "cm", device = "pdf")
ggsave("lifespan_shift.png", plot = lifespan_shift, width = 18, height = 16, units = "cm", device = "png")

lifespan_shift
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-70-1.png" width="768" />


**FIGURE 2**

```r
K_lifespan
```

<div class="figure">
<img src="mainresults_review_files/figure-html/figure2-1.png" alt="\label{fig:figure2}K through age for Healthy Controls (A) raw data, (B) after harmonizing (removing the estimated residual from the Linear Mixed Model) from T, AT, and AE resulting in the harmonized values of K." width="681.6" />
<p class="caption">\label{fig:figure2}K through age for Healthy Controls (A) raw data, (B) after harmonizing (removing the estimated residual from the Linear Mixed Model) from T, AT, and AE resulting in the harmonized values of K.</p>
</div>

## Baseline and rate estimates and Diagnostic discrimination based on changing rates


```r
# dados_datasetscomp$Diagnostic <- as.factor(dados_datasetscomp$Diagnostic)
# dados_datasetscomp$Diagnostic <- relevel(dados_datasetscomp$Diagnostic, ref = "CTL")
```


### Slope


```r
summary(lm(
  1 / 2 * log10(AvgThickness_shiftc) + log10(TotalArea_shiftc) ~ log10(ExposedArea_shiftc),
  data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere",
       Diagnostic == "CTL"),
  na.action = na.omit
))
```

```
## 
## Call:
## lm(formula = 1/2 * log10(AvgThickness_shiftc) + log10(TotalArea_shiftc) ~ 
##     log10(ExposedArea_shiftc), data = filter(dados_datasetscomp, 
##     Sample != "Mota&Houzel2015", ROI == "hemisphere", Diagnostic == 
##         "CTL"), na.action = na.omit)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.126471 -0.016037  0.002491  0.017389  0.093314 
## 
## Coefficients:
##                            Estimate Std. Error t value Pr(>|t|)    
## (Intercept)               -1.195278   0.040456  -29.55   <2e-16 ***
## log10(ExposedArea_shiftc)  1.394751   0.008818  158.18   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.02417 on 6800 degrees of freedom
## Multiple R-squared:  0.7863,	Adjusted R-squared:  0.7863 
## F-statistic: 2.502e+04 on 1 and 6800 DF,  p-value: < 2.2e-16
```

```r
## Displays confidence interval
# tidy(lm(
#   1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
#   data = dados_hemi_v1,
#   na.action = na.omit
# ), conf.int = TRUE)
```


```r
lm_Age <-
  filter(
    dados_datasetscomp,
    Sample != "Mota&Houzel2015",
    ROI == "hemisphere",
    Diagnostic == "CTL",
    Age_interval != "[0,5)",
    Age_interval != "[5,10)",
    Age_interval != "[10,15)"
  ) %>%
  group_by(Age_interval10) %>%
  do(fit_Age = tidy(
    lm(
      1 / 2 * log10(AvgThickness_shiftc) +  log10(TotalArea_shiftc) ~  log10(ExposedArea_shiftc),
      data = .,
      na.action = na.omit
    ),
    conf.int = TRUE
  )) %>%
  unnest(cols = c(fit_Age))

lm_Age %>%
  kable(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Age_interval10 </th>
   <th style="text-align:left;"> term </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> std.error </th>
   <th style="text-align:right;"> statistic </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:right;"> conf.low </th>
   <th style="text-align:right;"> conf.high </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> [10,20) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.89 </td>
   <td style="text-align:right;"> 0.19 </td>
   <td style="text-align:right;"> -4.76 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> -1.26 </td>
   <td style="text-align:right;"> -0.52 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [10,20) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.33 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 32.78 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.25 </td>
   <td style="text-align:right;"> 1.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [20,30) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.56 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> -13.85 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> -0.64 </td>
   <td style="text-align:right;"> -0.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [20,30) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.26 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 144.01 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.24 </td>
   <td style="text-align:right;"> 1.28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [30,40) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.53 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:right;"> -8.81 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> -0.65 </td>
   <td style="text-align:right;"> -0.41 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [30,40) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.25 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 95.40 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.23 </td>
   <td style="text-align:right;"> 1.28 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [40,50) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.18 </td>
   <td style="text-align:right;"> 0.11 </td>
   <td style="text-align:right;"> -1.74 </td>
   <td style="text-align:right;"> 0.08 </td>
   <td style="text-align:right;"> -0.39 </td>
   <td style="text-align:right;"> 0.02 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [40,50) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.17 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 51.29 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.13 </td>
   <td style="text-align:right;"> 1.22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [50,60) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.37 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:right;"> -3.72 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> -0.56 </td>
   <td style="text-align:right;"> -0.17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [50,60) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.21 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 55.97 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.17 </td>
   <td style="text-align:right;"> 1.26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [60,70) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> 0.21 </td>
   <td style="text-align:right;"> 0.09 </td>
   <td style="text-align:right;"> 2.46 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:right;"> 0.38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [60,70) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.08 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 57.57 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.05 </td>
   <td style="text-align:right;"> 1.12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [70,80) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.05 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:right;"> -0.69 </td>
   <td style="text-align:right;"> 0.49 </td>
   <td style="text-align:right;"> -0.19 </td>
   <td style="text-align:right;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [70,80) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.14 </td>
   <td style="text-align:right;"> 0.02 </td>
   <td style="text-align:right;"> 72.06 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.11 </td>
   <td style="text-align:right;"> 1.17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [80,90) </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> -0.16 </td>
   <td style="text-align:right;"> 0.12 </td>
   <td style="text-align:right;"> -1.34 </td>
   <td style="text-align:right;"> 0.18 </td>
   <td style="text-align:right;"> -0.40 </td>
   <td style="text-align:right;"> 0.08 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [80,90) </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.16 </td>
   <td style="text-align:right;"> 0.03 </td>
   <td style="text-align:right;"> 43.60 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1.11 </td>
   <td style="text-align:right;"> 1.22 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [90,100] </td>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:right;"> 0.16 </td>
   <td style="text-align:right;"> 0.39 </td>
   <td style="text-align:right;"> 0.40 </td>
   <td style="text-align:right;"> 0.69 </td>
   <td style="text-align:right;"> -0.64 </td>
   <td style="text-align:right;"> 0.95 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> [90,100] </td>
   <td style="text-align:left;"> log10(ExposedArea_shiftc) </td>
   <td style="text-align:right;"> 1.09 </td>
   <td style="text-align:right;"> 0.08 </td>
   <td style="text-align:right;"> 12.84 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.92 </td>
   <td style="text-align:right;"> 1.26 </td>
  </tr>
</tbody>
</table>

```r
lm_Age <- lm_Age %>% mutate(Age_interval10 = as.double((str_sub(Age_interval10,2,3))))

cor.test(filter(lm_Age, term == "log10(ExposedArea_shiftc)")$estimate, filter(lm_Age, term == "log10(ExposedArea_shiftc)")$Age_interval10)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(lm_Age, term == "log10(ExposedArea_shiftc)")$estimate and filter(lm_Age, term == "log10(ExposedArea_shiftc)")$Age_interval10
## t = -4.7078, df = 7, p-value = 0.002189
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.9727213 -0.4931597
## sample estimates:
##        cor 
## -0.8717632
```


```r
amostra_Coef <- filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL", Age_interval != "[0,5)", Age_interval != "[95,100]") %>%
    group_by(Age_interval10) %>%
    do(fit_amostra = tidy(lm(1/2 * log10(AvgThickness_shiftc) + log10(TotalArea_shiftc) ~ log10(ExposedArea_shiftc), data = ., na.action = na.omit), conf.int = TRUE, conf.level = 0.95)) %>% 
    unnest(fit_amostra)

amostras_Coef_age <- filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL", Age_interval != "[0,5)", Age_interval != "[95,100]") %>%
    group_by(Age_interval10) %>%
    summarise(
        N = n_distinct(SUBJ),
        min_age = min(Age),
        max_age = max(Age))

amostra_Coef <- full_join(amostra_Coef, amostras_Coef_age) %>% filter(term == "log10(ExposedArea_shiftc)")

fig_slope_ageinterval <- ggplot(amostra_Coef, aes(y = estimate, x = Age_interval10)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    geom_point() +
  geom_line(group = 1) + 
  geom_hline(yintercept = 1.25, linetype = "dashed") +
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9)
    ) + labs(y ="Slope (after harmonization)", x = "Age interval") +  geom_text(aes(label = N), nudge_y = 0.2)

ggsave("fig_slope_ageinterval.pdf", plot = fig_slope_ageinterval, width = 18, height = 14, units = "cm", device = "pdf")
```

**FIGURE 3**

```r
fig_slope_ageinterval
```

<div class="figure">
<img src="mainresults_review_files/figure-html/figure3-1.png" alt="\label{fig:figure3}Slope alpha derived from the Cortical Folding model from Mota &amp; Houzel through age for Healthy Controls. Data points from (0, 5] and [95, 100] years old were omitted due to the reduced data. The numbers on top of each point indicates the number of subjects at each interval. Mean value for all ages were 1.38+-0.012 and correlation of alpha and Age is Pearson's r~=~-0.69, p~=~0.0031 from [15,20) years old." width="681.6" />
<p class="caption">\label{fig:figure3}Slope alpha derived from the Cortical Folding model from Mota & Houzel through age for Healthy Controls. Data points from (0, 5] and [95, 100] years old were omitted due to the reduced data. The numbers on top of each point indicates the number of subjects at each interval. Mean value for all ages were 1.38+-0.012 and correlation of alpha and Age is Pearson's r~=~-0.69, p~=~0.0031 from [15,20) years old.</p>
</div>


```r
summary(lm(
  1 / 2 * log10(AvgThickness_shiftc) + log10(TotalArea_shiftc) ~ log10(ExposedArea_shiftc),
  data = filter(dados_datasetscomp,
       Sample != "Mota&Houzel2015",
       ROI == "hemisphere",
       Diagnostic == "AD"),
  na.action = na.omit
))
```

```
## 
## Call:
## lm(formula = 1/2 * log10(AvgThickness_shiftc) + log10(TotalArea_shiftc) ~ 
##     log10(ExposedArea_shiftc), data = filter(dados_datasetscomp, 
##     Sample != "Mota&Houzel2015", ROI == "hemisphere", Diagnostic == 
##         "AD"), na.action = na.omit)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.105523 -0.015849  0.001627  0.017815  0.052604 
## 
## Coefficients:
##                           Estimate Std. Error t value Pr(>|t|)    
## (Intercept)               -0.49777    0.08615  -5.778 9.81e-09 ***
## log10(ExposedArea_shiftc)  1.23286    0.01884  65.434  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.02369 on 1108 degrees of freedom
## Multiple R-squared:  0.7944,	Adjusted R-squared:  0.7942 
## F-statistic:  4282 on 1 and 1108 DF,  p-value: < 2.2e-16
```

```r
## Displays confidence interval
# tidy(lm(
#   1 / 2 * log10(AvgThickness) + log10(TotalArea) ~ log10(ExposedArea),
#   data = dados_hemi_v1,
#   na.action = na.omit
# ), conf.int = TRUE)
```


```r
amostra_Coef <- filter(dados_datasetscomp, ROI == "hemisphere", Age_interval != "[0,5)", Age_interval != "[95,100]") %>%
    group_by(Diagnostic, Age_interval10) %>%
    do(fit_amostra = tidy(lm(1/2 * log10(AvgThickness_shiftc) + log10(TotalArea_shiftc) ~ log10(ExposedArea_shiftc), data = ., na.action = na.omit), conf.int = TRUE, conf.level = 0.95)) %>% 
    unnest(fit_amostra)

amostras_Coef_age <- filter(dados_datasetscomp, ROI == "hemisphere", Age_interval != "[0,5)", Age_interval != "[95,100]") %>%
    group_by(Diagnostic, Age_interval10) %>%
    summarise(
        N = n_distinct(SUBJ),
        min_age = min(Age),
        max_age = max(Age))

amostra_Coef <- full_join(amostra_Coef, amostras_Coef_age) %>% filter(term == "log10(ExposedArea_shiftc)")

fig_slope_ageinterval_AD <- ggplot(amostra_Coef, aes(y = estimate, x = Age_interval10, color = Diagnostic)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    geom_point() +
  geom_line(aes(group=Diagnostic)) + 
  geom_hline(yintercept = 1.25, linetype = "dashed") +
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9)
    ) + labs(y ="Slope (after harmonization)", x = "Age interval") +  geom_text(aes(label = N), nudge_y = 0.2) +
  scale_fill_manual(values=cbbPalette2) +
  scale_colour_manual(values=cbbPalette2)

fig_slope_ageinterval_AD
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-76-1.png" width="768" />

```r
ggsave("fig_slope_ageinterval_AD.pdf", plot = fig_slope_ageinterval_AD, width = 18, height = 14, units = "cm", device = "pdf")
```

### GM volume \~ Age


```r
m.1 <- lme4::lmer(GMvolume_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## GMvolume_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: 183907.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.8431 -0.6738 -0.0306  0.6297  4.5365 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept) 0.00e+00     0   
##  Sample            (Intercept) 0.00e+00     0   
##  Residual                      7.33e+08 27074   
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                   Estimate Std. Error t value
## (Intercept)      303346.92     740.16 409.837
## Age               -1062.07      14.43 -73.583
## DiagnosticAD     -72189.07    7734.48  -9.333
## Age:DiagnosticAD    732.24     102.94   7.113
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")
pairs <- pairs(m.lstk) 

mean <- as_tibble(mean)
m.lstk <- as_tibble(m.lstk)


lst <- m.lstk %>%
  mutate(TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableGM <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "GM_shiftc")
 
coefs <- data.frame(coef(summary(m.1)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
```

```
##                     Estimate Std..Error    t.value          p.z
## (Intercept)      303346.9208  740.16472 409.837041 0.000000e+00
## Age               -1062.0694   14.43361 -73.583075 0.000000e+00
## DiagnosticAD     -72189.0707 7734.47570  -9.333415 0.000000e+00
## Age:DiagnosticAD    732.2418  102.94214   7.113140 1.134426e-12
```

```r
pairs %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> -732.2418 </td>
   <td style="text-align:right;"> 102.9421 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -7.11314 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
ggplot(data = as_tibble(pairs), aes(
    x = contrast,
    y = estimate,
    ymin = estimate - SE,
    ymax = estimate + SE)) +
    geom_hline(yintercept = 0,
               linetype = "11",
               colour = "grey60") +
    geom_pointrange( position = position_dodge(width = 0.2)) +
   geom_text(aes(label = str_c("p adj ", signif(p.value, digits = 2))), nudge_x = 0.3, nudge_y = 0.0003) +
    coord_flip() + 
    labs(y =  expression('Differences in slopes ('*Delta*'GM/'*Delta*'Age)'), x = "Contrast") +
  theme_pubr() + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10), text = element_text(size = 10))
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-77-1.png" width="768" />

```r
e <- effect("Age:Diagnostic", m.1)
e <- as.data.frame(e)

ggplot(e, aes(x = Age, y = fit,ymin=lower, ymax=upper, shape = Diagnostic)) + geom_pointrange() + geom_line() + theme_pubr() + labs(x = "Age [years]", y = "K", shape = "Diagnostic") + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10), text = element_text(size = 10))
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-77-2.png" width="768" />

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##    lsmean
##     <dbl>
## 1 250189.
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 27000
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### AvgThickness \~ Age


```r
m.1 <- lme4::lmer(AvgThickness_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: AvgThickness_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 |  
##     Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -13173.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -7.8059 -0.6167  0.0103  0.6500  3.4058 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept) 0.00000  0.0000  
##  Sample            (Intercept) 0.00000  0.0000  
##  Residual                      0.01101  0.1049  
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)       2.704e+00  2.869e-03 942.562
## Age              -4.440e-03  5.594e-05 -79.376
## DiagnosticAD     -3.695e-01  2.998e-02 -12.325
## Age:DiagnosticAD  3.029e-03  3.990e-04   7.592
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")

mean <- as_tibble(mean)
m.lstk <- as_tibble(m.lstk)

lst <- m.lstk %>%
  mutate(TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableT <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "T_shiftc")
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   2.48
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.1
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### TotalArea \~ Age


```r
m.1 <- lme4::lmer(TotalArea_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## TotalArea_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: 167586.5
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.9088 -0.7169 -0.0470  0.6641  4.2805 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept)        0    0    
##  Sample            (Intercept)        0    0    
##  Residual                      93066204 9647    
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)      113025.783    263.740 428.550
## Age                -244.375      5.143 -47.515
## DiagnosticAD     -13212.654   2755.994  -4.794
## Age:DiagnosticAD    152.779     36.681   4.165
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
sigma.hat(m.1)
```

```
## $sigma
## $sigma$data
## [1] 9647.083
## 
## $sigma$`Sample:Diagnostic`
## (Intercept) 
##           0 
## 
## $sigma$Sample
## (Intercept) 
##           0 
## 
## 
## $cors
## $cors$data
## [1] NA
## 
## $cors$`Sample:Diagnostic`
## [1] NA
## 
## $cors$Sample
## [1] NA
```

```r
r.squaredGLMM(m.1)
```

```
##            R2m       R2c
## [1,] 0.2795186 0.2795186
```

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")

mean <- as_tibble(mean)
m.lstk <- as_tibble(m.lstk)

 lst <- m.lstk %>%
   mutate(
     TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableAT <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "AT_shiftc")
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##    lsmean
##     <dbl>
## 1 100795.
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 9600
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### ExposedArea \~ Age


```r
m.1 <- lme4::lmer(ExposedArea_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## ExposedArea_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: 148540.8
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -2.8881 -0.7093 -0.0605  0.6726  3.8448 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev. 
##  Sample:Diagnostic (Intercept) 1.753e-13 4.187e-07
##  Sample            (Intercept) 3.654e-13 6.045e-07
##  Residual                      8.372e+06 2.893e+03
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                   Estimate Std. Error t value
## (Intercept)      40729.774     79.104 514.890
## Age                -41.214      1.543 -26.718
## DiagnosticAD      -321.787    826.609  -0.389
## Age:DiagnosticAD     2.208     11.002   0.201
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")

mean <- as_tibble(mean)
m.lstk <- as_tibble(m.lstk)

 lst <- m.lstk %>% mutate(
     TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableAE <- full_join(mean, lst, by = c("Diagnostic")) %>% mutate(percent = Age.trend*100/abs(lsmean), Variable = "AE_shiftc")
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1 38667.
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 2900
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##       6e-07
```

### GI \~ Age


```r
m.1 <- lme4::lmer(localGI_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## localGI_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -14560.8
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.7703 -0.6798 -0.0048  0.6759  6.4327 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept) 0.000000 0.00000 
##  Sample            (Intercept) 0.000000 0.00000 
##  Residual                      0.009238 0.09611 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error t value
## (Intercept)       2.778e+00  2.628e-03 1057.14
## Age              -3.509e-03  5.124e-05  -68.49
## DiagnosticAD     -3.243e-01  2.746e-02  -11.81
## Age:DiagnosticAD  3.811e-03  3.655e-04   10.43
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")

mean <- as_tibble(mean)
m.lstk <- as_tibble(m.lstk)

lst <- m.lstk %>%
  mutate(
    TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableGI <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "localGI_shiftc")
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   2.60
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.096
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### K \~ Age


```r
m.1 <- lme4::lmer(K_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: K_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -42446
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.5322 -0.6041  0.0391  0.6431  5.2736 
## 
## Random effects:
##  Groups            Name        Variance  Std.Dev.
##  Sample:Diagnostic (Intercept) 0.0000000 0.00000 
##  Sample            (Intercept) 0.0000000 0.00000 
##  Residual                      0.0002717 0.01648 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error  t value
## (Intercept)      -4.916e-01  4.507e-04 -1090.92
## Age              -8.600e-04  8.788e-06   -97.86
## DiagnosticAD     -8.706e-02  4.709e-03   -18.49
## Age:DiagnosticAD  8.940e-04  6.268e-05    14.26
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
coefs <- data.frame(coef(summary(m.1)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
```

```
##                       Estimate   Std..Error     t.value p.z
## (Intercept)      -0.4916413892 4.506654e-04 -1090.92322   0
## Age              -0.0008599766 8.788218e-06   -97.85563   0
## DiagnosticAD     -0.0870584281 4.709304e-03   -18.48647   0
## Age:DiagnosticAD  0.0008939612 6.267857e-05    14.26263   0
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1 -0.535
```

Mean values:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
pairs_means <- pairs(mean) 
mean <- as_tibble(mean)

pairs_means %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> 0.0423146 </td>
   <td style="text-align:right;"> 0.0016447 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 25.72712 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
paste("Difference in pct ", as_tibble(pairs_means)[2]*100/mean$lsmean[1])
```

```
## [1] "Difference in pct  -7.91394490809058"
```

Trends:

```r
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")
pairs_trends <- pairs(m.lstk) 
m.lstk <- as_tibble(m.lstk)

pairs_trends %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> -0.000894 </td>
   <td style="text-align:right;"> 6.27e-05 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -14.26263 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
lst <- m.lstk %>%
  mutate(
    TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableK <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "K_shiftc")
```

Figures:

```r
fig_TX_K_age_c <- ggplot(data = as_tibble(pairs), aes(
    x = contrast,
    y = estimate,
    ymin = estimate - SE,
    ymax = estimate + SE)) +
    geom_hline(yintercept = 0,
               linetype = "11",
               colour = "grey60") +
    geom_pointrange( position = position_dodge(width = 0.2)) +
   geom_text(aes(label = str_c("p adj ", signif(p.value, digits = 2))), nudge_x = 0.3, nudge_y = 0.0003) +
    coord_flip() + 
    labs(y =  expression('Differences in slopes ('*Delta*'K/'*Delta*'Age)'), x = "Contrast") +
  theme_pubr() + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10), text = element_text(size = 10))

e <- effect("Age:Diagnostic", m.1)
e <- as.data.frame(e)

fig_K_age_c <- ggplot(e, aes(x = Age, y = fit,ymin=lower, ymax=upper, shape = Diagnostic)) + geom_pointrange() + geom_line() + theme_pubr() + labs(x = "Age [years]", y = "K", shape = "Diagnostic") + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10), text = element_text(size = 10))
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.016
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### S \~ Age


```r
m.1 <- lme4::lmer(S_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: S_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -10771.6
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -3.4726 -0.6906 -0.0334  0.6362  5.8254 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept) 0.00000  0.0000  
##  Sample            (Intercept) 0.00000  0.0000  
##  Residual                      0.01492  0.1221  
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error  t value
## (Intercept)       9.088e+00  3.339e-03 2721.724
## Age               1.606e-03  6.511e-05   24.671
## DiagnosticAD      2.073e-01  3.489e-02    5.943
## Age:DiagnosticAD -1.339e-03  4.644e-04   -2.884
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
coefs <- data.frame(coef(summary(m.1)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
```

```
##                      Estimate   Std..Error     t.value          p.z
## (Intercept)       9.087638690 3.338928e-03 2721.723676 0.000000e+00
## Age               0.001606360 6.511088e-05   24.671142 0.000000e+00
## DiagnosticAD      0.207341356 3.489069e-02    5.942599 2.805383e-09
## Age:DiagnosticAD -0.001339213 4.643782e-04   -2.883885 3.928021e-03
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   9.17
```

Mean values:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
pairs_means <- pairs(mean) 
mean <- as_tibble(mean)

pairs_means %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> -0.1403122 </td>
   <td style="text-align:right;"> 0.0121857 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -11.51445 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
paste("Difference in pct ", as_tibble(pairs_means)[2]*100/mean$lsmean[1])
```

```
## [1] "Difference in pct  -1.53044891669301"
```

Trends:

```r
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")
pairs_trends <- pairs(m.lstk) 
m.lstk <- as_tibble(m.lstk)

pairs_trends %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> 0.0013392 </td>
   <td style="text-align:right;"> 0.0004644 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 2.883885 </td>
   <td style="text-align:right;"> 0.003928 </td>
  </tr>
</tbody>
</table>

```r
lst <- m.lstk %>% mutate(
    TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableS <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "S_shiftc")
```

Figures:


```r
fig_TX_S_age_c <- ggplot(data = as_tibble(pairs), aes(
    x = contrast,
    y = estimate,
    ymin = estimate - SE,
    ymax = estimate + SE)) +
    geom_hline(yintercept = 0,
               linetype = "11",
               colour = "grey60") +
    geom_pointrange( position = position_dodge(width = 0.2)) +
   geom_text(aes(label = str_c("p adj ", signif(p.value, digits = 2))), nudge_x = 0.3, nudge_y = 0.0003) +
    coord_flip() + 
    labs(y =  expression('Differences in slopes ('*Delta*'S/'*Delta*'Age)'), x = "Contrast") +
  theme_pubr() + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10), text = element_text(size = 10))

e <- effect("Age:Diagnostic", m.1)
e <- as.data.frame(e)

fig_S_age_c <- ggplot(e, aes(x = Age, y = fit,ymin=lower, ymax=upper, shape = Diagnostic)) + geom_pointrange() + geom_line() + theme_pubr() + labs(x = "Age [years]", y = "S", shape = "Diagnostic") + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    text = element_text(size = 10), legend.position = "none")
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.12
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### I \~ Age


```r
m.1 <- lme4::lmer(I_shiftc ~ Age * Diagnostic +(1|Sample) + (1|Sample:Diagnostic) , data = filter(dados_datasetscomp, ROI == "hemisphere"))

summary(m.1)
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: I_shiftc ~ Age * Diagnostic + (1 | Sample) + (1 | Sample:Diagnostic)
##    Data: filter(dados_datasetscomp, ROI == "hemisphere")
## 
## REML criterion at convergence: -16982.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.7921 -0.6605  0.0035  0.6753  3.2170 
## 
## Random effects:
##  Groups            Name        Variance Std.Dev.
##  Sample:Diagnostic (Intercept) 0.000000 0.00000 
##  Sample            (Intercept) 0.000000 0.00000 
##  Residual                      0.006801 0.08247 
## Number of obs: 7912, groups:  Sample:Diagnostic, 13; Sample, 11
## 
## Fixed effects:
##                    Estimate Std. Error  t value
## (Intercept)       1.053e+01  2.255e-03 4669.930
## Age              -3.074e-03  4.397e-05  -69.929
## DiagnosticAD     -1.938e-01  2.356e-02   -8.226
## Age:DiagnosticAD  1.702e-03  3.136e-04    5.426
## 
## Correlation of Fixed Effects:
##             (Intr) Age    DgnsAD
## Age         -0.896              
## DiagnostcAD -0.096  0.086       
## Ag:DgnstcAD  0.126 -0.140 -0.992
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

```r
coefs <- data.frame(coef(summary(m.1)))
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs
```

```
##                      Estimate   Std..Error     t.value          p.z
## (Intercept)      10.528806334 2.254596e-03 4669.930202 0.000000e+00
## Age              -0.003074501 4.396583e-05  -69.929319 0.000000e+00
## DiagnosticAD     -0.193800114 2.355978e-02   -8.225888 2.220446e-16
## Age:DiagnosticAD  0.001701526 3.135693e-04    5.426315 5.752943e-08
```

Mean value:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
as_tibble(mean)[1,2]
```

```
## # A tibble: 1 x 1
##   lsmean
##    <dbl>
## 1   10.4
```

Mean values:

```r
mean <- lsmeans(m.1, ~ Diagnostic, var="Age")
pairs_means <- pairs(mean) 
mean <- as_tibble(mean)

pairs_means %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> 0.1086367 </td>
   <td style="text-align:right;"> 0.0082284 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 13.20271 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

```r
paste("Difference in pct ", as_tibble(pairs_means)[2]*100/mean$lsmean[1])
```

```
## [1] "Difference in pct  1.04710882243449"
```

Trends:

```r
m.lstk <- lstrends(m.1, ~ Diagnostic, var="Age")
pairs_trends <- pairs(m.lstk) 
m.lstk <- as_tibble(m.lstk)

pairs_trends %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL - AD </td>
   <td style="text-align:right;"> -0.0017015 </td>
   <td style="text-align:right;"> 0.0003136 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -5.426315 </td>
   <td style="text-align:right;"> 1e-07 </td>
  </tr>
</tbody>
</table>

```r
lst <- m.lstk %>% mutate(
    TX = paste(signif(Age.trend, 2), "±", signif(SE, 2)))

tableI <- full_join(mean, lst, by = c("Diagnostic")) %>%
  mutate(percent = Age.trend*100/abs(lsmean), Variable = "I_shiftc")
```

Figures:


```r
fig_TX_I_age_c <- ggplot(data = as_tibble(pairs), aes(
    x = contrast,
    y = estimate,
    ymin = estimate - SE,
    ymax = estimate + SE)) +
    geom_hline(yintercept = 0,
               linetype = "11",
               colour = "grey60") +
    geom_pointrange( position = position_dodge(width = 0.2)) +
   geom_text(aes(label = str_c("p adj ", signif(p.value, digits = 2))), nudge_x = 0.3, nudge_y = 0.0003) +
    coord_flip() + 
    labs(y =  expression('Differences in slopes ('*Delta*'I/'*Delta*'Age)'), x = "Contrast") +
  theme_pubr() + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10), text = element_text(size = 10))

e <- effect("Age:Diagnostic", m.1)
e <- as.data.frame(e)

fig_I_age_c <- ggplot(e, aes(x = Age, y = fit,ymin=lower, ymax=upper, shape = Diagnostic)) + geom_pointrange() + geom_line() + theme_pubr() + labs(x = "Age [years]", y = "I", shape = "Diagnostic") + 
  theme(axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    text = element_text(size = 10), legend.position = "none")
```

Natural fluctuation error:

```r
signif(sigma.hat(m.1)$sigma$data,2)
```

```
## [1] 0.082
```

Type B error:

```r
signif(sigma.hat(m.1)$sigma$Sample,2)
```

```
## (Intercept) 
##           0
```

### rates



**TABLE 3**


```r
rate %>%
  kable(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Variable </th>
   <th style="text-align:left;"> CTL </th>
   <th style="text-align:left;"> AD </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> T_shiftc </td>
   <td style="text-align:left;"> -0.0044 ± 5.6e-05 ( -0.18 ) </td>
   <td style="text-align:left;"> -0.0014 ± 4e-04 ( -0.062 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AT_shiftc </td>
   <td style="text-align:left;"> -240 ± 5.1 ( -0.24 ) </td>
   <td style="text-align:left;"> -92 ± 36 ( -0.096 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AE_shiftc </td>
   <td style="text-align:left;"> -41 ± 1.5 ( -0.11 ) </td>
   <td style="text-align:left;"> -39 ± 11 ( -0.1 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> K_shiftc </td>
   <td style="text-align:left;"> -0.00086 ± 8.8e-06 ( -0.16 ) </td>
   <td style="text-align:left;"> 3.4e-05 ± 6.2e-05 ( 0.0059 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S_shiftc </td>
   <td style="text-align:left;"> 0.0016 ± 6.5e-05 ( 0.018 ) </td>
   <td style="text-align:left;"> 0.00027 ± 0.00046 ( 0.0029 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> I_shiftc </td>
   <td style="text-align:left;"> -0.0031 ± 4.4e-05 ( -0.03 ) </td>
   <td style="text-align:left;"> -0.0014 ± 0.00031 ( -0.013 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> localGI_shiftc </td>
   <td style="text-align:left;"> -0.0035 ± 5.1e-05 ( -0.13 ) </td>
   <td style="text-align:left;"> 3e-04 ± 0.00036 ( 0.012 ) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GM_shiftc </td>
   <td style="text-align:left;"> -1100 ± 14 ( -0.42 ) </td>
   <td style="text-align:left;"> -330 ± 100 ( -0.15 ) </td>
  </tr>
</tbody>
</table>

## Correlation between K, S, and I with Age for the Healthy Control group

### K 


```r
figS6a_alt <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"),
       aes(x = Age , y = K_shiftc)) +
  geom_point(aes(color = Sample, fill = Sample, alpha = 0.4)) +
  geom_smooth(color = "black", method = "lm") +
  # stat_cor() +
  theme_pubr() +
  guides(alpha = "none")+ 
  labs(x = "Age [years]", y ="K (after harmonization)") +
  scale_x_continuous(limits = c(0,100))

cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$K_shiftc, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$K_shiftc and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -107.66, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.8024991 -0.7849177
## sample estimates:
##        cor 
## -0.7938742
```

### S


```r
figS6b_alt <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"),
       aes(x = Age, y = S_shiftc)) +
  geom_point(aes(color = Sample, fill = Sample, alpha = 0.4)) +
  geom_smooth(color = "black", method = "lm") +
  # stat_cor() +
  theme_pubr() +
  guides(alpha = "none", color = FALSE,fill = FALSE)+ 
  labs(x = "Age [years]", y ="S (after harmonization)") +
  theme(legend.position= "none") +
  scale_x_continuous(limits = c(0,100))

cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$S_shiftc, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$S_shiftc and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = 25.606, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2747209 0.3180738
## sample estimates:
##       cor 
## 0.2965501
```

### I


```r
figS6c_alt <- ggplot(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL"),
       aes(x = Age, y = I_shiftc)) +
  geom_point(aes(color = Sample, fill = Sample, alpha = 0.4)) + 
  geom_smooth(color = "black", method = "lm") +
  # stat_cor() +
  theme_pubr() +
  guides(alpha = "none", color = FALSE,fill = FALSE)+ 
  labs(x = "Age [years]", y ="I (after harmonization)") +
  theme(legend.position= "none") +
  scale_x_continuous(limits = c(0,100))

cor.test(filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$I_shiftc, filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age)
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$I_shiftc and filter(dados_datasetscomp, ROI == "hemisphere", Diagnostic == "CTL")$Age
## t = -73.86, df = 6800, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.6801680 -0.6537884
## sample estimates:
##        cor 
## -0.6671873
```


```r
figS6_alt <- ggarrange(figS6a_alt, ggarrange(figS6b_alt, figS6c_alt,labels = c("B","C"), nrow = 1, ncol = 2, font.label = list(size = 11)),labels = c("A"), nrow = 2, ncol = 1, font.label = list(size = 11), common.legend = TRUE, legend = "top")

ggsave("figS6_alt.pdf", plot = figS6_alt, width = 18, height = 22, units = "cm", device = "pdf")
```

**FIGURE 1_alternativa**

```r
figS6_alt
```

<div class="figure">
<img src="mainresults_review_files/figure-html/figure1alt-1.png" alt="\label{fig:figure1}" width="681.6" />
<p class="caption">\label{fig:figure1}</p>
</div>

### Figures


```r
fig_age_c <- ggarrange(fig_K_age_c, ggarrange(fig_S_age_c, fig_I_age_c,labels = c("B","C"), nrow = 1, ncol = 2, font.label = list(size = 11)),labels = c("A"), nrow = 2, ncol = 1, font.label = list(size = 11), common.legend = TRUE, legend = "top")

ggsave("fig_age_c.pdf", plot = fig_age_c, width = 18, height = 18, units = "cm", device = "pdf")
```




**FIGURE4**


```r
fig_age_c
```

<div class="figure">
<img src="mainresults_review_files/figure-html/figure4-1.png" alt="\label{fig:figure4}Fitted values extracted from the linear regression model after data harmonization. Bars represent the 95 confidence interval. Complete summary of linear models in Supplementary Material. (A) Concerning K, Alzheimer's Disease (AD) has a shallow slope, meaning small changes with Age. Compared to healthy controls, AD values of K are almost constant and similar to older subjects (AD slope p-value&lt;0.0001 and CTL slope p-value&lt;0.0001; pairwise comparison estimate -0.000903, p&lt;0.0001). (B) For S, the AD and CTL patterns have similar intercepts, but statistically different slopes (AD slope p-value=0.004 and CTL slope p-value&lt;0.0001; pairwise comparison estimate 0.0013, p=0.004). (C) For I, that reflects brain volume, CTL has a decreasing volume with aging, while AD has a smaller slope (AD slope p-value&lt;0.0001 and CTL slope p-value&lt;0.0001; pairwise comparison estimate -0.00168, p&lt;0.0001)." width="681.6" />
<p class="caption">\label{fig:figure4}Fitted values extracted from the linear regression model after data harmonization. Bars represent the 95 confidence interval. Complete summary of linear models in Supplementary Material. (A) Concerning K, Alzheimer's Disease (AD) has a shallow slope, meaning small changes with Age. Compared to healthy controls, AD values of K are almost constant and similar to older subjects (AD slope p-value<0.0001 and CTL slope p-value<0.0001; pairwise comparison estimate -0.000903, p<0.0001). (B) For S, the AD and CTL patterns have similar intercepts, but statistically different slopes (AD slope p-value=0.004 and CTL slope p-value<0.0001; pairwise comparison estimate 0.0013, p=0.004). (C) For I, that reflects brain volume, CTL has a decreasing volume with aging, while AD has a smaller slope (AD slope p-value<0.0001 and CTL slope p-value<0.0001; pairwise comparison estimate -0.00168, p<0.0001).</p>
</div>

## Lobes


```r
dados_datasetscomp_lobes <- filter(dados_datasetscomp, ROI != "hemisphere", Diagnostic == "CTL" | Diagnostic == "AD", Sample != "IDOR-CCD-Control")
dados_datasetscomp_lobes$Diagnostic <- factor(dados_datasetscomp_lobes$Diagnostic, levels = c("CTL", "MCI","AD"))
dados_datasetscomp_lobes$Sample <- as.factor(dados_datasetscomp_lobes$Sample)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Sample </th>
   <th style="text-align:left;"> Diagnostic </th>
   <th style="text-align:right;"> N </th>
   <th style="text-align:left;"> age </th>
   <th style="text-align:left;"> age_range </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ADNI </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 868 </td>
   <td style="text-align:left;"> 75 ± 6.5 </td>
   <td style="text-align:left;"> 56 ;  96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ADNI </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:right;"> 542 </td>
   <td style="text-align:left;"> 75 ± 8 </td>
   <td style="text-align:left;"> 56 ;  92 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHEAD </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:left;"> 42 ± 19 </td>
   <td style="text-align:left;"> 24 ;  76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AOMICPIOP1 </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 208 </td>
   <td style="text-align:left;"> 22 ± 1.8 </td>
   <td style="text-align:left;"> 18 ;  26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AOMICPIOP2 </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 224 </td>
   <td style="text-align:left;"> 22 ± 1.8 </td>
   <td style="text-align:left;"> 18 ;  26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HCPr900 </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 881 </td>
   <td style="text-align:left;"> 29 ± 3.6 </td>
   <td style="text-align:left;"> 24 ;  37 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IDOR </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 77 </td>
   <td style="text-align:left;"> 66 ± 8.4 </td>
   <td style="text-align:left;"> 43 ;  80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IDOR </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:left;"> 77 ± 6.1 </td>
   <td style="text-align:left;"> 63 ;  86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IXI-Guys </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 314 </td>
   <td style="text-align:left;"> 51 ± 16 </td>
   <td style="text-align:left;"> 20 ;  86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IXI-HH </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 181 </td>
   <td style="text-align:left;"> 47 ± 17 </td>
   <td style="text-align:left;"> 20 ;  82 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IXI-IOP </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:left;"> 42 ± 17 </td>
   <td style="text-align:left;"> 20 ;  86 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NKI </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:left;"> 34 ± 19 </td>
   <td style="text-align:left;"> 4 ;  85 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> OASIS </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 312 </td>
   <td style="text-align:left;"> 45 ± 24 </td>
   <td style="text-align:left;"> 18 ;  94 </td>
  </tr>
</tbody>
</table>

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Diagnostic </th>
   <th style="text-align:right;"> N </th>
   <th style="text-align:left;"> age </th>
   <th style="text-align:left;"> age_range </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:right;"> 3095 </td>
   <td style="text-align:left;"> 46 ± 23 </td>
   <td style="text-align:left;"> 4 ;  96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:right;"> 555 </td>
   <td style="text-align:left;"> 75 ± 8 </td>
   <td style="text-align:left;"> 56 ;  92 </td>
  </tr>
</tbody>
</table>

### K \~ Age

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: K_shiftc ~ Age * Diagnostic * ROI + (1 | Sample:Diagnostic:ROI)
##    Data: dados_datasetscomp_lobes
## 
## REML criterion at convergence: -139166.2
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -6.3529 -0.6169  0.0203  0.6414  5.3896 
## 
## Random effects:
##  Groups                Name        Variance  Std.Dev.
##  Sample:Diagnostic:ROI (Intercept) 0.0000000 0.00000 
##  Residual                          0.0004934 0.02221 
## Number of obs: 29188, groups:  Sample:Diagnostic:ROI, 48
## 
## Fixed effects:
##                         Estimate Std. Error  t value
## (Intercept)           -5.062e-01  6.403e-04 -790.494
## Age                   -9.198e-04  1.247e-05  -73.752
## DiagnosticAD          -9.008e-02  6.343e-03  -14.201
## ROIO                   4.222e-02  9.055e-04   46.621
## ROIP                   4.204e-02  9.055e-04   46.429
## ROIT                   2.639e-02  9.055e-04   29.142
## Age:DiagnosticAD       9.744e-04  8.447e-05   11.536
## Age:ROIO               2.355e-04  1.764e-05   13.352
## Age:ROIP              -2.060e-04  1.764e-05  -11.680
## Age:ROIT              -3.828e-05  1.764e-05   -2.170
## DiagnosticAD:ROIO      3.362e-02  8.970e-03    3.748
## DiagnosticAD:ROIP     -4.009e-02  8.970e-03   -4.469
## DiagnosticAD:ROIT      6.636e-03  8.970e-03    0.740
## Age:DiagnosticAD:ROIO -4.070e-04  1.195e-04   -3.407
## Age:DiagnosticAD:ROIP  4.587e-04  1.195e-04    3.840
## Age:DiagnosticAD:ROIT -2.046e-04  1.195e-04   -1.713
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-128-1.png" width="768" /><img src="mainresults_review_files/figure-html/unnamed-chunk-128-2.png" width="768" />

Trends:

```r
pairs %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> Contrast.1 </th>
   <th style="text-align:left;"> Contrast.2 </th>
   <th style="text-align:left;"> Diagnostic.1 </th>
   <th style="text-align:left;"> ROI.1 </th>
   <th style="text-align:left;"> Diagnostic.2 </th>
   <th style="text-align:left;"> ROI.2 </th>
   <th style="text-align:left;"> contrast_ROI </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0009744 </td>
   <td style="text-align:right;"> 8.45e-05 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -11.535618 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> CTL F </td>
   <td style="text-align:left;"> AD F </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> F-F </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0005674 </td>
   <td style="text-align:right;"> 8.45e-05 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -6.716977 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> CTL O </td>
   <td style="text-align:left;"> AD O </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> O-O </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0014331 </td>
   <td style="text-align:right;"> 8.45e-05 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -16.966067 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> CTL P </td>
   <td style="text-align:left;"> AD P </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> P-P </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0007698 </td>
   <td style="text-align:right;"> 8.45e-05 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -9.113069 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> CTL T </td>
   <td style="text-align:left;"> AD T </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> T-T </td>
  </tr>
</tbody>
</table>

### S \~ Age


```
## Linear mixed model fit by REML ['lmerMod']
## Formula: S_shiftc ~ Age * Diagnostic * ROI + (1 | Sample:Diagnostic:ROI)
##    Data: dados_datasetscomp_lobes
## 
## REML criterion at convergence: -29380.8
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.2192 -0.6718 -0.0240  0.6368  5.8473 
## 
## Random effects:
##  Groups                Name        Variance Std.Dev.
##  Sample:Diagnostic:ROI (Intercept) 0.00000  0.0000  
##  Residual                          0.02126  0.1458  
## Number of obs: 29188, groups:  Sample:Diagnostic:ROI, 48
## 
## Fixed effects:
##                         Estimate Std. Error  t value
## (Intercept)            8.731e+00  4.203e-03 2076.999
## Age                    2.681e-03  8.187e-05   32.752
## DiagnosticAD           2.785e-01  4.164e-02    6.689
## ROIO                  -2.149e-01  5.945e-03  -36.156
## ROIP                   6.513e-01  5.945e-03  109.566
## ROIT                  -4.981e-02  5.945e-03   -8.379
## Age:DiagnosticAD      -2.415e-03  5.545e-04   -4.355
## Age:ROIO              -1.212e-03  1.158e-04  -10.464
## Age:ROIP              -1.198e-04  1.158e-04   -1.035
## Age:ROIT              -5.236e-04  1.158e-04   -4.522
## DiagnosticAD:ROIO     -1.154e-01  5.889e-02   -1.960
## DiagnosticAD:ROIP      2.092e-01  5.889e-02    3.552
## DiagnosticAD:ROIT     -2.767e-02  5.889e-02   -0.470
## Age:DiagnosticAD:ROIO  9.338e-04  7.842e-04    1.191
## Age:DiagnosticAD:ROIP -2.442e-03  7.842e-04   -3.114
## Age:DiagnosticAD:ROIT  9.578e-04  7.842e-04    1.221
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-130-1.png" width="768" /><img src="mainresults_review_files/figure-html/unnamed-chunk-130-2.png" width="768" />

Trends:

```r
pairs %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> Contrast.1 </th>
   <th style="text-align:left;"> Contrast.2 </th>
   <th style="text-align:left;"> Diagnostic.1 </th>
   <th style="text-align:left;"> ROI.1 </th>
   <th style="text-align:left;"> Diagnostic.2 </th>
   <th style="text-align:left;"> ROI.2 </th>
   <th style="text-align:left;"> contrast_ROI </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> 0.0024148 </td>
   <td style="text-align:right;"> 0.0005545 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 4.354852 </td>
   <td style="text-align:right;"> 0.0003550 </td>
   <td style="text-align:left;"> CTL F </td>
   <td style="text-align:left;"> AD F </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> F-F </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> 0.0014810 </td>
   <td style="text-align:right;"> 0.0005545 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 2.670780 </td>
   <td style="text-align:right;"> 0.1315406 </td>
   <td style="text-align:left;"> CTL O </td>
   <td style="text-align:left;"> AD O </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> O-O </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> 0.0048567 </td>
   <td style="text-align:right;"> 0.0005545 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 8.758395 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:left;"> CTL P </td>
   <td style="text-align:left;"> AD P </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> P-P </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> 0.0014570 </td>
   <td style="text-align:right;"> 0.0005545 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 2.627543 </td>
   <td style="text-align:right;"> 0.1459513 </td>
   <td style="text-align:left;"> CTL T </td>
   <td style="text-align:left;"> AD T </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> T-T </td>
  </tr>
</tbody>
</table>

### I \~ Age


```
## Linear mixed model fit by REML ['lmerMod']
## Formula: I_shiftc ~ Age * Diagnostic * ROI + (1 | Sample:Diagnostic:ROI)
##    Data: dados_datasetscomp_lobes
## 
## REML criterion at convergence: -47911.9
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -4.5924 -0.6391  0.0025  0.6538  4.0292 
## 
## Random effects:
##  Groups                Name        Variance Std.Dev.
##  Sample:Diagnostic:ROI (Intercept) 0.00000  0.0000  
##  Residual                          0.01127  0.1061  
## Number of obs: 29188, groups:  Sample:Diagnostic:ROI, 48
## 
## Fixed effects:
##                         Estimate Std. Error  t value
## (Intercept)            1.029e+01  3.060e-03 3363.538
## Age                   -3.634e-03  5.959e-05  -60.975
## DiagnosticAD          -2.590e-01  3.031e-02   -8.546
## ROIO                  -7.612e-01  4.327e-03 -175.911
## ROIP                   3.175e-01  4.327e-03   73.373
## ROIT                   7.542e-02  4.327e-03   17.429
## Age:DiagnosticAD       2.767e-03  4.036e-04    6.855
## Age:ROIO               7.472e-04  8.428e-05    8.866
## Age:ROIP               5.915e-04  8.428e-05    7.018
## Age:ROIT               6.682e-04  8.428e-05    7.928
## DiagnosticAD:ROIO      1.045e-01  4.286e-02    2.437
## DiagnosticAD:ROIP     -2.042e-03  4.286e-02   -0.048
## DiagnosticAD:ROIT      8.302e-03  4.286e-02    0.194
## Age:DiagnosticAD:ROIO -1.215e-03  5.708e-04   -2.129
## Age:DiagnosticAD:ROIP  5.997e-05  5.708e-04    0.105
## Age:DiagnosticAD:ROIT -7.663e-04  5.708e-04   -1.342
## optimizer (nloptwrap) convergence code: 0 (OK)
## boundary (singular) fit: see help('isSingular')
```

<img src="mainresults_review_files/figure-html/unnamed-chunk-132-1.png" width="768" /><img src="mainresults_review_files/figure-html/unnamed-chunk-132-2.png" width="768" />

Trends:

```r
pairs %>%
  kable() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> contrast </th>
   <th style="text-align:right;"> estimate </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> z.ratio </th>
   <th style="text-align:right;"> p.value </th>
   <th style="text-align:left;"> Contrast.1 </th>
   <th style="text-align:left;"> Contrast.2 </th>
   <th style="text-align:left;"> Diagnostic.1 </th>
   <th style="text-align:left;"> ROI.1 </th>
   <th style="text-align:left;"> Diagnostic.2 </th>
   <th style="text-align:left;"> ROI.2 </th>
   <th style="text-align:left;"> contrast_ROI </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0027668 </td>
   <td style="text-align:right;"> 0.0004036 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -6.854835 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:left;"> CTL F </td>
   <td style="text-align:left;"> AD F </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> F </td>
   <td style="text-align:left;"> F-F </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0015513 </td>
   <td style="text-align:right;"> 0.0004036 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -3.843424 </td>
   <td style="text-align:right;"> 0.0030554 </td>
   <td style="text-align:left;"> CTL O </td>
   <td style="text-align:left;"> AD O </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> O </td>
   <td style="text-align:left;"> O-O </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0028267 </td>
   <td style="text-align:right;"> 0.0004036 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -7.003407 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:left;"> CTL P </td>
   <td style="text-align:left;"> AD P </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> P </td>
   <td style="text-align:left;"> P-P </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CTL-AD </td>
   <td style="text-align:right;"> -0.0020005 </td>
   <td style="text-align:right;"> 0.0004036 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> -4.956275 </td>
   <td style="text-align:right;"> 0.0000198 </td>
   <td style="text-align:left;"> CTL T </td>
   <td style="text-align:left;"> AD T </td>
   <td style="text-align:left;"> CTL </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> AD </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> T-T </td>
  </tr>
</tbody>
</table>


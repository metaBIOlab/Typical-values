# PREPARACAO ----
#source("analises/preparing-yujiang-datasets.R")
#source("analises/preparing-yujiang-datasets-agedecay_25.R")
source("analises/preparing-Motahouzel-dataset.R")
source("Wang2016_OASIS.R")
source("Wang2016_HCP.R")
source("Wang2019_ADNI_SEMCSF.R")
source("Wang2019_HCP2.R")
source("Wang2019_IXI2.R")
source("Wang2019_NKI.R")
source("IDOR.R")
source("AOMIC.R")
# source("ACC.R")
# source("ACC_caltechucsf.R")
source("AHEAD.R")
# source("ZK.R")
# source("Micro.R")
source("AOMIC-PIOP2.R")
# source("HCPINFANT.R")
# source("ds000030.R")

# junta as amostras para comparacao ----
dados_datasetscomp <-
  full_join(dados_idor, dados_AOMIC) %>%
  full_join(dados_AHEAD) %>%
  full_join(dados_AOMICPIOP2)
  # %>% full_join(dados_ds000030) # %>% full_join(dados_HPC)
# %>% full_join(dados_ACC) %>% full_join(dados_caltechucsf) %>% full_join(dados_zk) %>% full_join(dados_micro)

HCPr900$SUBJ <- as.character(HCPr900$SUBJ)
HCPr900$Age <- as.double(HCPr900$Age)

dados_datasetscomp <- dados_datasetscomp %>%
  mutate(
    hemi = Hemisphere,
    AvgThickness = AvgCortThickness,
    TotalArea = PialArea,
    ExposedArea = SmoothPialArea
  ) 

source("analises/variables.R")

Wang <- full_join(ADNI, HCPr900) %>%
  full_join(NKI) %>%
  full_join(OASIS) %>%
  full_join(IXI) %>% full_join(HCP500r)

source("analises/variables2.R")

dados_datasetscomp <- full_join(dados_datasetscomp, Wang)

dados_datasetscomp$ROI <- as.factor(dados_datasetscomp$ROI)
dados_datasetscomp$Diagnostic <- as.factor(dados_datasetscomp$Diagnostic)
dados_datasetscomp$Sample <- as.factor(dados_datasetscomp$Sample)
dados_datasetscomp$Gender <- as.factor(dados_datasetscomp$Gender)

dados_datasetscomp <- dados_datasetscomp %>%
  dplyr::select(
    -c(Lobe,
      SubjectID,
      Session,
      birthday,
      acq_date,
      Lobe,
      AvgCortThickness,
      PialArea,
      WhiteArea,
      SmoothPialArea,
      GreymatterVol,
      BMI,
      handedness,
      ESC,
      education_category,
      religious_now,
      raven_score,
      NEO_N,
      NEO_E,
      NEO_O,
      NEO_A,
      NEO_C,
      Hemisphere
    )
  )  %>%
  filter(
    !is.na(TotalArea),
    ExposedArea != 0,!is.na(localGI),!is.infinite(AvgThickness),!is.na(Diagnostic), AvgThickness !=0,
    localGI != 0,!is.infinite(AvgThickness),!is.na(AvgThickness)
    ) %>%
  droplevels() %>%
  unique()

# dados_datasetscomp$Diagnostic <- revalue(dados_datasetscomp$Diagnostic, c("CONTROLE"="CTL", "Control" = "CTL", "ALZ"="AD", "CCL" = "MCI"))
# dados_datasetscomp$Diagnostic <- revalue(dados_datasetscomp$Diagnostic, c("CONTROLE"="CTL", "ALZ"="AD", "CCL" = "MCI", "AgCC" = "CCD"))

dados_datasetscomp$Diagnostic <- as.character(dados_datasetscomp$Diagnostic)

dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "CONTROLE"] <- "CTL"
dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "Control"] <- "CTL"
dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "CONTROL"] <- "CTL"
dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "CCL"] <- "MCI"
dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "ALZ"] <- "AD"
dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "AgCC"] <- "CCD"
dados_datasetscomp$Diagnostic[dados_datasetscomp$Diagnostic == "MICRO"] <- "MICRO-ZIKA"

dados_datasetscomp$Diagnostic <- factor(dados_datasetscomp$Diagnostic)

dados_datasetscomp$Gender[dados_datasetscomp$Gender == "F"] <- "FEM"
dados_datasetscomp$Gender[dados_datasetscomp$Gender == "M"] <- "MASC"
dados_datasetscomp$Gender[dados_datasetscomp$Gender == "f"] <- "FEM"
dados_datasetscomp$Gender[dados_datasetscomp$Gender == "m"] <- "MASC"
#dados_datasetscomp$Gender[dados_datasetscomp$Gender == "u"] <- "U"

dados_datasetscomp <- filter(dados_datasetscomp, !is.nan(Age), !is.na(Age), Gender != "u") %>%
  droplevels()

# decaimento da idade ----
# source("analises/decaimento_idade_25.R")
# source("analises/deaging.R")

source("harm_deageg_lobes.R")

dados_datasetscomp$machine[dados_datasetscomp$Sample == "HCP900r"] <- "Siemens-Skyra-modified-3T"
dados_datasetscomp$FieldStrenght[dados_datasetscomp$Sample == "HCP900r"] <- 3

dados_datasetscomp$machine[dados_datasetscomp$Sample == "NKI"] <- "Siemens-Magnetom-3T"
dados_datasetscomp$FieldStrenght[dados_datasetscomp$Sample == "NKI"] <- 3

dados_datasetscomp$machine[dados_datasetscomp$Sample == "OASIS"] <- "Siemens-Vision-1.5T"
dados_datasetscomp$FieldStrenght[dados_datasetscomp$Sample == "OASIS"] <- 1.5

dados_datasetscomp$machine[dados_datasetscomp$Sample == "ADNI"] <- "multiple-3T"
dados_datasetscomp$FieldStrenght[dados_datasetscomp$Sample == "ADNI"] <- 3

dados_datasetscomp$machine[dados_datasetscomp$Sample == "IXI"] <- "Philips-Gyroscan Intera-1.5T"
dados_datasetscomp$FieldStrenght[dados_datasetscomp$Sample == "IXI"] <- 1.5

# ADNI,
# HCPr900,
# NKI,
# OASIS,
# AHEAD,
# AOMIC,
# AOMICPIOP2,


rm(
  HCPr900_SmoothArea,
  HCPr900_SmoothArea_CHwoB,
  HCPr900_subj,
  HCPr900_T,
  HCPr900_TotalArea,
  all_sessions,
  data_caltechucsf,
  data_lobes_caltechucsf,
  dataAHEAD,
  dataAHEADlobes,
  dataAOMIC,
  dataAOMIClobes,
  part_ses,
  participants,
  tabela_sujeitos,
  tabela_sujeitos_caltechucsf,
  UCSF_CALTECH_IDOR_AgCCDatabase,
  data,
  data_lobes,
  decay_AvgThickness,
  decay_logExposedArea,
  decay_logTotalArea,
  IXI_SmoothArea,
  IXI_SmoothArea_CHwoB,
  IXI_subj,
  IXI_T,
  IXI_TotalArea,
  NKI_SmoothArea,
  NKI_SmoothArea_CHwoB,
  NKI_subj,
  NKI_T,
  NKI_TotalArea,
  OASIS_lh,
  OASIS_rh,
  ADNI_SmoothArea,
  ADNI_SmoothArea_CHwoB,
  ADNI_subj,
  ADNI_T,
  ADNI_TotalArea
)

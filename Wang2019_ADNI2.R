# ADNI ####
## Hemisphere #####
ADNI <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI.csv")
ADNI <-
  ADNI %>% mutate(Sample = "ADNI", Diagnostic = Diagnosis, ROI = "hemisphere") 

colnames(ADNI)[which(names(ADNI) == "SubjID")] <- "SUBJ"

ADNI_subj <- ADNI %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic)

ADNI_lh <- ADNI %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic, ROI, contains("_1")) %>% mutate(hemi = "L")
ADNI_rh <- ADNI %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic, ROI, contains("_2")) %>% mutate(hemi = "R")

colnames(ADNI_lh)[which(names(ADNI_lh) == "PialArea_1")] <-
  "TotalArea"

colnames(ADNI_lh)[which(names(ADNI_lh) == "SmoothPialArea_1")] <-
  "ExposedArea"

colnames(ADNI_lh)[which(names(ADNI_lh) == "AvgCortThickness_1")] <-
  "AvgThickness"

colnames(ADNI_lh)[which(names(ADNI_lh) == "WhiteArea_1")] <-
  "WhiteArea"

colnames(ADNI_lh)[which(names(ADNI_lh) == "PialFullArea_1")] <-
  "PialFullArea"

colnames(ADNI_lh)[which(names(ADNI_lh) == "WhiteFullArea_1")] <-
  "WhiteFullArea"

colnames(ADNI_lh)[which(names(ADNI_lh) == "SmoothPialFullArea_1")] <-
  "SmoothPialFullArea"

colnames(ADNI_lh)[which(names(ADNI_lh) == "ConvexHullFullArea_1")] <-
  "ConvexHullFullArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "PialArea_2")] <-
  "TotalArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "SmoothPialArea_2")] <-
  "ExposedArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "AvgCortThickness_2")] <-
  "AvgThickness"

colnames(ADNI_rh)[which(names(ADNI_rh) == "WhiteArea_2")] <-
  "WhiteArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "PialFullArea_2")] <-
  "PialFullArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "WhiteFullArea_2")] <-
  "WhiteFullArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "SmoothPialFullArea_2")] <-
  "SmoothPialFullArea"

colnames(ADNI_rh)[which(names(ADNI_rh) == "ConvexHullFullArea_2")] <-
  "ConvexHullFullArea"

ADNI <- rbind(ADNI_lh, ADNI_rh)
ADNI$SUBJ <- as.character(ADNI$SUBJ)

## Lobes #####
### Total Area #####
ADNI_TotalArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_TotalArea_L.csv", 
                            col_names = FALSE)
ADNI_TotalArea_L <- cbind(ADNI_subj, ADNI_TotalArea_L)

ADNI_TotalArea_L <- ADNI_TotalArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

ADNI_TotalArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_TotalArea_R.csv", 
                            col_names = FALSE)
ADNI_TotalArea_R <- cbind(ADNI_subj, ADNI_TotalArea_R)

ADNI_TotalArea_R <- ADNI_TotalArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

ADNI_TotalArea <- rbind(ADNI_TotalArea_L, ADNI_TotalArea_R)
ADNI_TotalArea$SUBJ <- as.character(ADNI_TotalArea$SUBJ)

### Avg Thickness #####
ADNI_T_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_T_L.csv", 
                    col_names = FALSE)
ADNI_T_L <- cbind(ADNI_subj, ADNI_T_L)

ADNI_T_L <- ADNI_T_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

ADNI_T_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_T_R.csv", 
                    col_names = FALSE)
ADNI_T_R <- cbind(ADNI_subj, ADNI_T_R)

ADNI_T_R <- ADNI_T_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

ADNI_T <- rbind(ADNI_T_L, ADNI_T_R)
ADNI_T$SUBJ <- as.character(ADNI_T$SUBJ)

### Exposed Area #####
ADNI_SmoothArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_L.csv", 
                             col_names = FALSE)
ADNI_SmoothArea_L <- cbind(ADNI_subj, ADNI_SmoothArea_L)

ADNI_SmoothArea_L <- ADNI_SmoothArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea")

ADNI_SmoothArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_R.csv", 
                             col_names = FALSE)
ADNI_SmoothArea_R <- cbind(ADNI_subj, ADNI_SmoothArea_R)

ADNI_SmoothArea_R <- ADNI_SmoothArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea")

ADNI_SmoothArea <- rbind(ADNI_SmoothArea_L, ADNI_SmoothArea_R)
ADNI_SmoothArea$SUBJ <- as.character(ADNI_SmoothArea$SUBJ)

### Exposed Area CHwoB#####
ADNI_SmoothArea_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_CHwoB_L.csv", 
                                   col_names = FALSE)
ADNI_SmoothArea_CHwoB_L <- cbind(ADNI_subj, ADNI_SmoothArea_CHwoB_L)

ADNI_SmoothArea_CHwoB_L <- ADNI_SmoothArea_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea_CHwoB")

ADNI_SmoothArea_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_CHwoB_R.csv", 
                                   col_names = FALSE)
ADNI_SmoothArea_CHwoB_R <- cbind(ADNI_subj, ADNI_SmoothArea_CHwoB_R)

ADNI_SmoothArea_CHwoB_R <- ADNI_SmoothArea_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea_CHwoB")

ADNI_SmoothArea_CHwoB <- rbind(ADNI_SmoothArea_CHwoB_L, ADNI_SmoothArea_CHwoB_R)
ADNI_SmoothArea_CHwoB$SUBJ <- as.character(ADNI_SmoothArea_CHwoB$SUBJ)

## ADNI #####
ADNI <- full_join(ADNI_T, ADNI_SmoothArea) %>%
  full_join(ADNI_TotalArea) %>%
  full_join(ADNI_SmoothArea_CHwoB) %>%
  mutate(c = ExposedArea_CHwoB/ExposedArea) %>%
  full_join(ADNI)

rm(ADNI_SmoothArea_CHwoB_R,
   ADNI_SmoothArea_CHwoB_L,
   ADNI_SmoothArea_R,
   ADNI_SmoothArea_L,
   ADNI_TotalArea_L,
   ADNI_TotalArea_R,
   ADNI_T_L,
   ADNI_T_R,
   ADNI_lh,
   ADNI_rh)

ADNI$SUBJ <- as.character(ADNI$SUBJ)
ADNI$Age <- as.double(ADNI$Age)

ADNIMERGE <- read_csv(
  "~/idor/Gyrification/data/ADNIMERGE.csv",
  col_types = cols(
    PIB_bl = col_character(),
    ABETA = col_character(),
    ABETA_bl = col_character(),
    TAU = col_character(),
    PTAU = col_character(),
    PTAU_bl = col_character(),
    TAU_bl = col_character()
  )
) %>% mutate(SUBJ = str_c(
  PTID,
  "_",
  str_sub(EXAMDATE, 1, 4),
  str_sub(EXAMDATE, 6, 7),
  str_sub(EXAMDATE, 9, 10)
)) %>%
  filter(!is.na(TAU), !is.na(ABETA), !is.nan(TAU), !is.nan(ABETA)) %>%
  mutate(TAU = as.double(TAU),
         ABETA = ifelse(ABETA == ">1700", as.double(1700), as.double(ABETA)))

ADNI <- left_join(ADNI, ADNIMERGE)

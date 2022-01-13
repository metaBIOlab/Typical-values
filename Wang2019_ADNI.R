# ADNI ####
## Hemisphere #####
ADNI <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI.csv")
ADNI <-
  ADNI %>% mutate(Sample = "ADNI", Diagnostic = Diagnosis, ROI = "hemisphere") 

ADNI_subj <- ADNI %>%
  dplyr::select(SubjID, Sample, Gender, Age)

ADNI_lh <- ADNI %>%
  dplyr::select(SubjID, Sample, Gender, Age, Diagnostic, ROI, contains("_1")) %>% mutate(hemi = "L")
ADNI_rh <- ADNI %>%
  dplyr::select(SubjID, Sample, Gender, Age, Diagnostic, ROI, contains("_2")) %>% mutate(hemi = "R")

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

colnames(ADNI)[which(names(ADNI) == "SubjID")] <- "SUBJ"

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

### Exposed Area #####
ADNI_SmoothArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_L.csv", 
                             col_names = FALSE)
ADNI_SmoothArea_L <- cbind(ADNI_subj, ADNI_SmoothArea_L)

ADNI_SmoothArea_L <- ADNI_SmoothArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea")

ADNI_SmoothArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_R.csv", 
                             col_names = FALSE)
ADNI_SmoothArea_R <- cbind(ADNI_subj, ADNI_SmoothArea_R)

ADNI_SmoothArea_R <- ADNI_SmoothArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea")

ADNI_SmoothArea <- rbind(ADNI_SmoothArea_L, ADNI_SmoothArea_R)

ADNI_SmoothArea_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_CHwoB_L.csv", 
                                   col_names = FALSE)
ADNI_SmoothArea_CHwoB_L <- cbind(ADNI_subj, ADNI_SmoothArea_CHwoB_L)

ADNI_SmoothArea_CHwoB_L <- ADNI_SmoothArea_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea_CHwoB")


ADNI_SmoothArea_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/ADNI_SmoothArea_CHwoB_R.csv", 
                                   col_names = FALSE)
ADNI_SmoothArea_CHwoB_R <- cbind(ADNI_subj, ADNI_SmoothArea_CHwoB_R)

ADNI_SmoothArea_CHwoB_R <- ADNI_SmoothArea_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea_CHwoB")

ADNI_SmoothArea_CHwoB <- rbind(ADNI_SmoothArea_CHwoB_L, ADNI_SmoothArea_CHwoB_R)

## ADNI #####
ADNI <- full_join(ADNI, ADNI_SmoothArea) %>%
  full_join(ADNI_T) %>%
  full_join(ADNI_TotalArea) %>%
  filter(ROI != "X5", ROI != "X6") %>%
  dplyr::select(-c(SubjID))

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

ADNI$ROI[ADNI$ROI == "X1"] <- "F"
ADNI$ROI[ADNI$ROI == "X2"] <- "P"
ADNI$ROI[ADNI$ROI == "X3"] <- "T"
ADNI$ROI[ADNI$ROI == "X4"] <- "O"
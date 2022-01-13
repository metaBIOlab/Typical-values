# NKI ####
## Hemisphere #####
NKI <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI.csv")
NKI <-
  NKI %>% mutate(Sample = "NKI", Diagnostic = "CTL", ROI = "hemisphere") 

colnames(NKI)[which(names(NKI) == "SubjID")] <- "SUBJ"

NKI_subj <- NKI %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic)

NKI_lh <- NKI %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic, ROI, contains("_1")) %>% mutate(hemi = "L")
NKI_rh <- NKI %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic, ROI, contains("_2")) %>% mutate(hemi = "R")

colnames(NKI_lh)[which(names(NKI_lh) == "PialArea_1")] <-
  "TotalArea"

colnames(NKI_lh)[which(names(NKI_lh) == "SmoothPialArea_1")] <-
  "ExposedArea"

colnames(NKI_lh)[which(names(NKI_lh) == "AvgCortThickness_1")] <-
  "AvgThickness"

colnames(NKI_lh)[which(names(NKI_lh) == "WhiteArea_1")] <-
  "WhiteArea"

colnames(NKI_lh)[which(names(NKI_lh) == "PialFullArea_1")] <-
  "PialFullArea"

colnames(NKI_lh)[which(names(NKI_lh) == "WhiteFullArea_1")] <-
  "WhiteFullArea"

colnames(NKI_lh)[which(names(NKI_lh) == "SmoothPialFullArea_1")] <-
  "SmoothPialFullArea"

colnames(NKI_lh)[which(names(NKI_lh) == "ConvexHullFullArea_1")] <-
  "ConvexHullFullArea"

colnames(NKI_rh)[which(names(NKI_rh) == "PialArea_2")] <-
  "TotalArea"

colnames(NKI_rh)[which(names(NKI_rh) == "SmoothPialArea_2")] <-
  "ExposedArea"

colnames(NKI_rh)[which(names(NKI_rh) == "AvgCortThickness_2")] <-
  "AvgThickness"

colnames(NKI_rh)[which(names(NKI_rh) == "WhiteArea_2")] <-
  "WhiteArea"

colnames(NKI_rh)[which(names(NKI_rh) == "PialFullArea_2")] <-
  "PialFullArea"

colnames(NKI_rh)[which(names(NKI_rh) == "WhiteFullArea_2")] <-
  "WhiteFullArea"

colnames(NKI_rh)[which(names(NKI_rh) == "SmoothPialFullArea_2")] <-
  "SmoothPialFullArea"

colnames(NKI_rh)[which(names(NKI_rh) == "ConvexHullFullArea_2")] <-
  "ConvexHullFullArea"

NKI <- rbind(NKI_lh, NKI_rh)
NKI$SUBJ <- as.character(NKI$SUBJ)

## Lobes #####
### Total Area #####
NKI_TotalArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_TotalArea_L.csv", 
                                col_names = FALSE)
NKI_TotalArea_L <- cbind(NKI_subj, NKI_TotalArea_L)

NKI_TotalArea_L <- NKI_TotalArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

NKI_TotalArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_TotalArea_R.csv", 
                                col_names = FALSE)
NKI_TotalArea_R <- cbind(NKI_subj, NKI_TotalArea_R)

NKI_TotalArea_R <- NKI_TotalArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

NKI_TotalArea <- rbind(NKI_TotalArea_L, NKI_TotalArea_R)
NKI_TotalArea$SUBJ <- as.character(NKI_TotalArea$SUBJ)

### Avg Thickness #####
NKI_T_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_T_L.csv", 
                        col_names = FALSE)
NKI_T_L <- cbind(NKI_subj, NKI_T_L)

NKI_T_L <- NKI_T_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

NKI_T_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_T_R.csv", 
                        col_names = FALSE)
NKI_T_R <- cbind(NKI_subj, NKI_T_R)

NKI_T_R <- NKI_T_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

NKI_T <- rbind(NKI_T_L, NKI_T_R)
NKI_T$SUBJ <- as.character(NKI_T$SUBJ)

### Exposed Area #####
NKI_SmoothArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_SmoothArea_L.csv", 
                                 col_names = FALSE)
NKI_SmoothArea_L <- cbind(NKI_subj, NKI_SmoothArea_L)

NKI_SmoothArea_L <- NKI_SmoothArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea")

NKI_SmoothArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_SmoothArea_R.csv", 
                                 col_names = FALSE)
NKI_SmoothArea_R <- cbind(NKI_subj, NKI_SmoothArea_R)

NKI_SmoothArea_R <- NKI_SmoothArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea")

NKI_SmoothArea <- rbind(NKI_SmoothArea_L, NKI_SmoothArea_R)
NKI_SmoothArea$SUBJ <- as.character(NKI_SmoothArea$SUBJ)

### Exposed Area CHwoB#####
NKI_SmoothArea_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_SmoothArea_CHwoB_L.csv", 
                                       col_names = FALSE)
NKI_SmoothArea_CHwoB_L <- cbind(NKI_subj, NKI_SmoothArea_CHwoB_L)

NKI_SmoothArea_CHwoB_L <- NKI_SmoothArea_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea_CHwoB")

NKI_SmoothArea_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_SmoothArea_CHwoB_R.csv", 
                                       col_names = FALSE)
NKI_SmoothArea_CHwoB_R <- cbind(NKI_subj, NKI_SmoothArea_CHwoB_R)

NKI_SmoothArea_CHwoB_R <- NKI_SmoothArea_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea_CHwoB")

NKI_SmoothArea_CHwoB <- rbind(NKI_SmoothArea_CHwoB_L, NKI_SmoothArea_CHwoB_R)
NKI_SmoothArea_CHwoB$SUBJ <- as.character(NKI_SmoothArea_CHwoB$SUBJ)

### Gaussian Curvature #####
NKI_K_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_K_CHwoB_L.csv", 
                                   col_names = FALSE)
NKI_K_CHwoB_L <- cbind(NKI_subj, NKI_K_CHwoB_L)

NKI_K_CHwoB_L <- NKI_K_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "K_CHwoB")

NKI_K_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/NKI_K_CHwoB_R.csv", 
                                   col_names = FALSE)
NKI_K_CHwoB_R <- cbind(NKI_subj, NKI_K_CHwoB_R)

NKI_K_CHwoB_R <- NKI_K_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "K_CHwoB")

NKI_K_CHwoB <- rbind(NKI_K_CHwoB_L, NKI_K_CHwoB_R)
NKI_K_CHwoB$SUBJ <- as.character(NKI_K_CHwoB$SUBJ)

## NKI #####
NKI <- full_join(NKI_T, NKI_SmoothArea) %>%
  full_join(NKI_TotalArea) %>%
  full_join(NKI_K_CHwoB) %>%
  full_join(NKI_SmoothArea_CHwoB) %>%
  mutate(c = 4*pi/K_CHwoB) %>%
  full_join(NKI)

rm(NKI_SmoothArea_CHwoB_R,
   NKI_SmoothArea_CHwoB_L,
   NKI_SmoothArea_R,
   NKI_SmoothArea_L,
   NKI_TotalArea_L,
   NKI_TotalArea_R,
   NKI_T_L,
   NKI_T_R,
   NKI_lh,
   NKI_rh)

NKI$SUBJ <- as.character(NKI$SUBJ)
NKI$Age <- as.double(NKI$Age)

# HCPr900 ####
## Hemisphere #####
HCPr900 <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900.csv")
HCPr900 <-
  HCPr900 %>% mutate(Sample = "HCPr900", Diagnostic = "CTL", ROI = "hemisphere",
                     Age = ifelse(AgeCat == "26-30", 28,
                                  ifelse(AgeCat == "31-35", 33,
                                         ifelse(AgeCat == "22-25", 23.5,
                                                ifelse(AgeCat == "36+", 37, "")))))

colnames(HCPr900)[which(names(HCPr900) == "SubjID")] <- "SUBJ"

HCPr900_subj <- HCPr900 %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic)

HCPr900_lh <- HCPr900 %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic, ROI, contains("_1")) %>% mutate(hemi = "L")
HCPr900_rh <- HCPr900 %>%
  dplyr::select(SUBJ, Sample, Gender, Age, Diagnostic, ROI, contains("_2")) %>% mutate(hemi = "R")

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "PialArea_1")] <-
  "TotalArea"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "SmoothPialArea_1")] <-
  "ExposedArea"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "AvgCortThickness_1")] <-
  "AvgThickness"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "WhiteArea_1")] <-
  "WhiteArea"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "PialFullArea_1")] <-
  "PialFullArea"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "WhiteFullArea_1")] <-
  "WhiteFullArea"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "SmoothPialFullArea_1")] <-
  "SmoothPialFullArea"

colnames(HCPr900_lh)[which(names(HCPr900_lh) == "ConvexHullFullArea_1")] <-
  "ConvexHullFullArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "PialArea_2")] <-
  "TotalArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "SmoothPialArea_2")] <-
  "ExposedArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "AvgCortThickness_2")] <-
  "AvgThickness"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "WhiteArea_2")] <-
  "WhiteArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "PialFullArea_2")] <-
  "PialFullArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "WhiteFullArea_2")] <-
  "WhiteFullArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "SmoothPialFullArea_2")] <-
  "SmoothPialFullArea"

colnames(HCPr900_rh)[which(names(HCPr900_rh) == "ConvexHullFullArea_2")] <-
  "ConvexHullFullArea"

HCPr900 <- rbind(HCPr900_lh, HCPr900_rh)
HCPr900$SUBJ <- as.character(HCPr900$SUBJ)

## Lobes #####
### Total Area #####
HCPr900_TotalArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_TotalArea_L.csv", 
                            col_names = FALSE)
HCPr900_TotalArea_L <- cbind(HCPr900_subj, HCPr900_TotalArea_L)

HCPr900_TotalArea_L <- HCPr900_TotalArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

HCPr900_TotalArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_TotalArea_R.csv", 
                            col_names = FALSE)
HCPr900_TotalArea_R <- cbind(HCPr900_subj, HCPr900_TotalArea_R)

HCPr900_TotalArea_R <- HCPr900_TotalArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

HCPr900_TotalArea <- rbind(HCPr900_TotalArea_L, HCPr900_TotalArea_R)
HCPr900_TotalArea$SUBJ <- as.character(HCPr900_TotalArea$SUBJ)

### Avg Thickness #####
HCPr900_T_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_T_L.csv", 
                    col_names = FALSE)
HCPr900_T_L <- cbind(HCPr900_subj, HCPr900_T_L)

HCPr900_T_L <- HCPr900_T_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

HCPr900_T_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_T_R.csv", 
                    col_names = FALSE)
HCPr900_T_R <- cbind(HCPr900_subj, HCPr900_T_R)

HCPr900_T_R <- HCPr900_T_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

HCPr900_T <- rbind(HCPr900_T_L, HCPr900_T_R)
HCPr900_T$SUBJ <- as.character(HCPr900_T$SUBJ)

### Exposed Area #####
HCPr900_SmoothArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_SmoothArea_L.csv", 
                             col_names = FALSE)
HCPr900_SmoothArea_L <- cbind(HCPr900_subj, HCPr900_SmoothArea_L)

HCPr900_SmoothArea_L <- HCPr900_SmoothArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea")

HCPr900_SmoothArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_SmoothArea_R.csv", 
                             col_names = FALSE)
HCPr900_SmoothArea_R <- cbind(HCPr900_subj, HCPr900_SmoothArea_R)

HCPr900_SmoothArea_R <- HCPr900_SmoothArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea")

HCPr900_SmoothArea <- rbind(HCPr900_SmoothArea_L, HCPr900_SmoothArea_R)
HCPr900_SmoothArea$SUBJ <- as.character(HCPr900_SmoothArea$SUBJ)

### Exposed Area CHwoB#####
HCPr900_SmoothArea_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_SmoothArea_CHwoB_L.csv", 
                                   col_names = FALSE)
HCPr900_SmoothArea_CHwoB_L <- cbind(HCPr900_subj, HCPr900_SmoothArea_CHwoB_L)

HCPr900_SmoothArea_CHwoB_L <- HCPr900_SmoothArea_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea_CHwoB")

HCPr900_SmoothArea_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_SmoothArea_CHwoB_R.csv", 
                                   col_names = FALSE)
HCPr900_SmoothArea_CHwoB_R <- cbind(HCPr900_subj, HCPr900_SmoothArea_CHwoB_R)

HCPr900_SmoothArea_CHwoB_R <- HCPr900_SmoothArea_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "ExposedArea_CHwoB")

HCPr900_SmoothArea_CHwoB <- rbind(HCPr900_SmoothArea_CHwoB_L, HCPr900_SmoothArea_CHwoB_R)
HCPr900_SmoothArea_CHwoB$SUBJ <- as.character(HCPr900_SmoothArea_CHwoB$SUBJ)

### Gaussian Curvature #####
HCPr900_K_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_K_CHwoB_L.csv", 
                                       col_names = FALSE)
HCPr900_K_CHwoB_L <- cbind(HCPr900_subj, HCPr900_K_CHwoB_L)

HCPr900_K_CHwoB_L <- HCPr900_K_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "K_CHwoB")

HCPr900_K_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/HCPr900_K_CHwoB_R.csv", 
                                       col_names = FALSE)
HCPr900_K_CHwoB_R <- cbind(HCPr900_subj, HCPr900_K_CHwoB_R)

HCPr900_K_CHwoB_R <- HCPr900_K_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "K_CHwoB")

HCPr900_K_CHwoB <- rbind(HCPr900_K_CHwoB_L, HCPr900_K_CHwoB_R)
HCPr900_K_CHwoB$SUBJ <- as.character(HCPr900_K_CHwoB$SUBJ)

## HCPr900 #####
HCPr900 <- full_join(HCPr900_T, HCPr900_SmoothArea) %>%
  full_join(HCPr900_TotalArea) %>%
  full_join(HCPr900_SmoothArea_CHwoB) %>%
  full_join(HCPr900_K_CHwoB) %>%
  mutate(c = 4*pi/K_CHwoB) %>%
  full_join(HCPr900)

rm(HCPr900_SmoothArea_CHwoB_R,
   HCPr900_SmoothArea_CHwoB_L,
   HCPr900_SmoothArea_R,
   HCPr900_SmoothArea_L,
   HCPr900_TotalArea_L,
   HCPr900_TotalArea_R,
   HCPr900_T_L,
   HCPr900_T_R,
   HCPr900_lh,
   HCPr900_rh)

HCPr900$SUBJ <- as.character(HCPr900$SUBJ)
HCPr900$Age <- as.double(HCPr900$Age)

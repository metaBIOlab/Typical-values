# IXI ####
## Hemisphere #####
IXI <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI.csv")
IXI <-
  IXI %>% mutate(Sample = "IXI", Diagnostic = "CTL", ROI = "hemisphere") 

IXI_subj <- IXI %>%
  dplyr::select(SubjID, Sample, Gender, Age)

IXI_lh <- IXI %>%
  dplyr::select(SubjID, Sample, Gender, Age, Diagnostic, ROI, contains("_1")) %>% mutate(hemi = "L")
IXI_rh <- IXI %>%
  dplyr::select(SubjID, Sample, Gender, Age, Diagnostic, ROI, contains("_2")) %>% mutate(hemi = "R")

colnames(IXI_lh)[which(names(IXI_lh) == "PialArea_1")] <-
  "TotalArea"

colnames(IXI_lh)[which(names(IXI_lh) == "SmoothPialArea_1")] <-
  "ExposedArea"

colnames(IXI_lh)[which(names(IXI_lh) == "AvgCortThickness_1")] <-
  "AvgThickness"

colnames(IXI_lh)[which(names(IXI_lh) == "WhiteArea_1")] <-
  "WhiteArea"

colnames(IXI_lh)[which(names(IXI_lh) == "PialFullArea_1")] <-
  "PialFullArea"

colnames(IXI_lh)[which(names(IXI_lh) == "WhiteFullArea_1")] <-
  "WhiteFullArea"

colnames(IXI_lh)[which(names(IXI_lh) == "SmoothPialFullArea_1")] <-
  "SmoothPialFullArea"

colnames(IXI_lh)[which(names(IXI_lh) == "ConvexHullFullArea_1")] <-
  "ConvexHullFullArea"

colnames(IXI_rh)[which(names(IXI_rh) == "PialArea_2")] <-
  "TotalArea"

colnames(IXI_rh)[which(names(IXI_rh) == "SmoothPialArea_2")] <-
  "ExposedArea"

colnames(IXI_rh)[which(names(IXI_rh) == "AvgCortThickness_2")] <-
  "AvgThickness"

colnames(IXI_rh)[which(names(IXI_rh) == "WhiteArea_2")] <-
  "WhiteArea"

colnames(IXI_rh)[which(names(IXI_rh) == "PialFullArea_2")] <-
  "PialFullArea"

colnames(IXI_rh)[which(names(IXI_rh) == "WhiteFullArea_2")] <-
  "WhiteFullArea"

colnames(IXI_rh)[which(names(IXI_rh) == "SmoothPialFullArea_2")] <-
  "SmoothPialFullArea"

colnames(IXI_rh)[which(names(IXI_rh) == "ConvexHullFullArea_2")] <-
  "ConvexHullFullArea"

IXI <- rbind(IXI_lh, IXI_rh)

colnames(IXI)[which(names(IXI) == "SubjID")] <- "SUBJ"

## Lobes #####
### Total Area #####
IXI_TotalArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_TotalArea_L.csv", 
                            col_names = FALSE)
IXI_TotalArea_L <- cbind(IXI_subj, IXI_TotalArea_L)

IXI_TotalArea_L <- IXI_TotalArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

IXI_TotalArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_TotalArea_R.csv", 
                            col_names = FALSE)
IXI_TotalArea_R <- cbind(IXI_subj, IXI_TotalArea_R)

IXI_TotalArea_R <- IXI_TotalArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "TotalArea")

IXI_TotalArea <- rbind(IXI_TotalArea_L, IXI_TotalArea_R)

### Avg Thickness #####
IXI_T_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_T_L.csv", 
                    col_names = FALSE)
IXI_T_L <- cbind(IXI_subj, IXI_T_L)

IXI_T_L <- IXI_T_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

IXI_T_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_T_R.csv", 
                    col_names = FALSE)
IXI_T_R <- cbind(IXI_subj, IXI_T_R)

IXI_T_R <- IXI_T_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "AvgThickness")

IXI_T <- rbind(IXI_T_L, IXI_T_R)

### Exposed Area #####
IXI_SmoothArea_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_SmoothArea_L.csv", 
                             col_names = FALSE)
IXI_SmoothArea_L <- cbind(IXI_subj, IXI_SmoothArea_L)

IXI_SmoothArea_L <- IXI_SmoothArea_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea")

IXI_SmoothArea_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_SmoothArea_R.csv", 
                             col_names = FALSE)
IXI_SmoothArea_R <- cbind(IXI_subj, IXI_SmoothArea_R)

IXI_SmoothArea_R <- IXI_SmoothArea_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea")

IXI_SmoothArea <- rbind(IXI_SmoothArea_L, IXI_SmoothArea_R)

IXI_SmoothArea_CHwoB_L <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_SmoothArea_CHwoB_L.csv", 
                                   col_names = FALSE)
IXI_SmoothArea_CHwoB_L <- cbind(IXI_subj, IXI_SmoothArea_CHwoB_L)

IXI_SmoothArea_CHwoB_L <- IXI_SmoothArea_CHwoB_L  %>%
  mutate(hemi = "L") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea_CHwoB")

IXI_SmoothArea_CHwoB_R <- read_csv("~/idor/Gyrification/data/Wang2016e2019/IXI_SmoothArea_CHwoB_R.csv", 
                                   col_names = FALSE)
IXI_SmoothArea_CHwoB_R <- cbind(IXI_subj, IXI_SmoothArea_CHwoB_R)

IXI_SmoothArea_CHwoB_R <- IXI_SmoothArea_CHwoB_R %>%
  mutate(hemi = "R") %>%
  pivot_longer(X1:X6, names_to = "ROI", values_to = "SmoothArea_CHwoB")

IXI_SmoothArea_CHwoB <- rbind(IXI_SmoothArea_CHwoB_L, IXI_SmoothArea_CHwoB_R)

## IXI #####
IXI <- full_join(IXI, IXI_SmoothArea) %>%
  full_join(IXI_T) %>%
  full_join(IXI_TotalArea) %>%
  filter(ROI != "X5", ROI != "X6") %>%
  dplyr::select(-c(SubjID))

rm(IXI_SmoothArea_CHwoB_R,
   IXI_SmoothArea_CHwoB_L,
   IXI_SmoothArea_R,
   IXI_SmoothArea_L,
   IXI_TotalArea_L,
   IXI_TotalArea_R,
   IXI_T_L,
   IXI_T_R,
   IXI_lh,
   IXI_rh)

IXI$SUBJ <- as.character(IXI$SUBJ)
IXI$Age <- as.double(IXI$Age)
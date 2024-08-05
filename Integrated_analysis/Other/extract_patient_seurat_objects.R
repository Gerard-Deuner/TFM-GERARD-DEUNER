# get patient-specific seurat objects from the merged seurat object
data.dir <- "/Users/gerarddeunercos/Documents/SERPENTINE/data/inputdata/"
merged.sobj <- readRDS(paste0(data.dir, "SERPENTINE_raw_18_01_2024.rds"))

# patient 01
p01.sobj <- merged.sobj %>%
  subset(patient == "01") %>%
  saveRDS(file = paste0(data.dir, "Patient_01_SCR_C02_harmony_integrated_TCR_01-12-2023.rds"))

# patient 02
p02.sobj <- merged.sobj %>%
  subset(patient == "02") %>%
  saveRDS(file = paste0(data.dir, "Patient_02_SCR_C02_harmony_integrated_TCR_01-12-2023.2.rds"))

# patient 03
p03.sobj <- merged.sobj %>%
  subset(patient == "03") %>%
  saveRDS(file = paste0(data.dir, "/Patient_03_SCR_C02_harmony_integrated_TCR_01-12-2023.rds"))

# patient 08
p08.sobj <- merged.sobj %>%
  subset(patient == "08") %>%
  saveRDS(file = paste0(data.dir, "/Patient_08_SCR_C02_harmony_integrated_TCR_01-12-2023.rds"))

# patient 10
p10.sobj <- merged.sobj %>%
  subset(patient == "10") %>%
  saveRDS(file = paste0(data.dir, "/Patient_10_SCR_C02_harmony_integrated_TCR_01-12-2023.rds"))

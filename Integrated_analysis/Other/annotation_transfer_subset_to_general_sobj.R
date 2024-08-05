####################################################################
# TRANSFER ANNOTATION LABELS FROM SUBSETS TO GENERAL SEURAT OBJECT #
####################################################################

# load packages
library(dplyr)
library(qs)
library(Seurat)

# load env vars
source("/Users/gerarddeunercos/Documents/cluster1/gdeuner/config.R")

# patients vector
patients <- c("Patient_01", "Patient_02", "Patient_03", "Patient_08", "Patient_10")

# subsets vector
subsets <- c("TNK", "myeloid")

# iterate over patients
for (patient in patients){
  
  # load patient preprocessed seurat object
  print(paste("Loading", patient, "Seurat Object..."))
  data <- qread(file.path(sp.data.dir, "outputdata", paste0(patient, "_SCR_CO2_harmony_integrated_pp_annotated_TCR_22-01-24")))
  print(paste(patient, "Seurat Object Loaded!"))
  
  # add metadata column were the annotation will be stored
  data$Annotation_2.0 <- rep(NA, nrow(data@meta.data))
  
  # iterate over subsets
  for (subset in subsets){
    
    print(paste("Loading", patient, subset, "Subset Seurat Object"))
    
    # load subset seurat object
    if (subset == "TNK") {
      
      data.s <- qread(file.path(sp.data.dir, "outputdata", paste0(patient, "_SCR_CO2_", subset, "_subset_pp_annotated_no_clean_30-01-24")))
      
    } else {
      
      data.s <- qread(file.path(sp.data.dir, "outputdata", paste0(patient, "_SCR_CO2_", subset, "_subset_pp_annotated_no_clean_05-02-24")))
      
    }
    
    
    # transfer annotation 
    common.rows <- intersect(rownames(data.s@meta.data), rownames(data@meta.data))
    
    if (subset == "TNK") {
      
      tnk.labels <- data.s@meta.data %>% 
        dplyr::select(Annotation_2.0) 
    
    } else {
      
      mye.labels <- data.s@meta.data %>% 
        dplyr::select(Annotation_2.0) 
    
    }

    
    # add a metadata column indicating the subset (main cell type)
    data.s$cell_type <- rep(subset, nrow(data.s@meta.data))
    subset_column <- data.s@meta.data %>%
      select(cell_type)
    data <- AddMetaData(data, subset_column)
    
    
  }
  
  # add column in main object
  labels <- rbind(tnk.labels, mye.labels)
  data <- AddMetaData(data, labels)
  
  # save Annotation1.0 as Annotation2.0 for the other cell types 
  data$Annotation_2.0 <- ifelse(is.na(data$Annotation_2.0), as.character(data$Annotation_1.0), as.character(data$Annotation_2.0))
  # same for main cell type
  data$cell_type <- ifelse(is.na(data$cell_type), as.character(data$Annotation_1.0), as.character(data$cell_type))
  
  
  # save seurat object
  print("Saving Newly Annotated Seurat Object..")
  qsave(data, file.path(sp.data.dir, "outputdata", paste0(patient, "_SCR_CO2_harmony_integrated_pp_annotated_2.0_TCR_13-02-24")))
  print("Seurat Object Saved!")
  
  # remove seurat object of memory
  rm(data)
  
}

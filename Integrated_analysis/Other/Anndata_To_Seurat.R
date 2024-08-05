#####################
# Anndata To Seurat #
#####################

# load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)

# import path variables 
source("/Users/gerarddeunercos/Documents/cluster1/gdeuner/config.R")

# convert to h5seurat
Convert("/Users/gerarddeunercos/Downloads/Combined_SCR_CO2_TNK_annotated_13-04-24.h5ad", dest = "h5ad", overwrite = TRUE)


##############

library(sceasy)
library(reticulate)
use_condaenv('/Users/gerarddeunercos/Documents/cluster1/gdeuner/envs/')
loompy <- reticulate::import('loompy')

devtools::install_github('satijalab/seurat-data')

sceasy::convertFormat(file.path(sp.data.dir, "outputdata", "combined", "Combined_SCR_CO2_TCR_full-integrated_annot_22-03-24.h5ad"), from="anndata", to="seurat",
                      outFile='file.path(sp.data.dir, "outputdata", "combined", "Combined_SCR_CO2_TCR_full-integrated_annot_22-03-24.rds')

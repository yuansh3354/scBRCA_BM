rm(list = ls())
gc(reset = T)
# packages
library(cowplot)
library(Seurat)
library(tidyr)
library(viridis)
library(stringr)
library(scales)
library(magrittr)
library(reshape2)
library(numbat)
library(dplyr)
library(tibble)
options(stringsAsFactors = F)
plot.path='plts/'
source('script/Functions.R')
########################## My Functions ##########################
numbat_analysis <- function(sce.obj,sID) {
  input.file=paste0('result/OUT_Numbat/',sID,'/',sID,'_allele_counts.tsv.gz')
  out.dir=paste0('result/OUT_Numbat/',sID,'/')
  
  sce.obj = subset(sce.obj, subset = SampleID==sID)
  # Extract count matrix from Seurat object
  count_mat <- sce.obj@assays$RNA@counts
  colnames(count_mat) <- gsub("_\\d+",'',colnames(count_mat))
  # Import allele counts data and select overlapping cells
  df_allele <- read.table(input.file, header = TRUE, sep = "\t", quote = "", comment.char = "")
  ids <- intersect(df_allele$cell, colnames(count_mat))
  df_allele <- df_allele[which(df_allele$cell %in% c(ids)),] 
  count_mat = count_mat[,ids]
  
  use.cells = sce.obj@meta.data[which(sce.obj$celltype =='MTC'),] %>% rownames()
  use.cells = gsub("_\\d+",'',use.cells)
  use.cells = intersect(use.cells,ids)
  df_allele = df_allele[which(df_allele$cell %in% c(use.cells)),] 
  count_mat = count_mat[,use.cells]
  # Run numBat
  out <- run_numbat(count_mat, ref_hca, df_allele,
                    genome = "hg38", t = 1e-5, ncores = 32,
                    plot = T, out_dir = out.dir)
  
  return(out)
}
get_barcode_cell_annotation = function(sce.obj,sID){
  #ids=gsub(pattern = "[-]|230269A_", replacement = "", x = sID)
  sce.obj = subset(sce.obj, subset = SampleID==sID)
  barcode_cell_annotation = sce.obj@meta.data[,'Cell_type',drop=F]
  barcode_cell_annotation = tibble::add_column(barcode_cell_annotation, 
                                               Index = gsub('_\\d+','',rownames(barcode_cell_annotation)),
                                               .before = 1)
  colnames(barcode_cell_annotation) = c('Index','Cell_type')
  print(dim(barcode_cell_annotation))
  folder_path <- paste0("result/OUT_Numbat/", sID)
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
    # cat(sprintf("Folder '%s' was created.\n", folder_path))
  }
  write.table(barcode_cell_annotation, 
              file =paste0(folder_path,'/',sID,'.cell_barcode_annotations.tsv'),
              sep = "\t", row.names = FALSE,quote = F)
  return(barcode_cell_annotation)
}
Combined.SeuratObj = readRDS('Combined.SeuratObj')
unique(Combined.SeuratObj$SampleID)
MTC.meta = read.csv('MTC.SubCell.csv',row.names = 1)
MTC.meta = MTC.meta[,c('cell_Subtype'),drop=F]
Combined.SeuratObj = AddMetaData(Combined.SeuratObj, metadata = MTC.meta)
use.cells = rownames(Combined.SeuratObj@meta.data[which(Combined.SeuratObj$celltype %in%c('Oligodendrocytes','MTC','Microglia','T Cell','Macrophage')),])
Combined.SeuratObj = subset(Combined.SeuratObj,cells=use.cells)
SampleIDs = c('230269A_JFP01_cellranger710','230269A_JFP02_cellranger710','GSM4259357_cellranger710',
              'LHP-307-001-10XSC3_cellranger710','LXZ-307-001-10XSC3_cellranger710','MZH-307-001-10XSC3_cellranger710',
              'ZMY-307-001-10XSC3_cellranger710','ZMY-307-002-10XSC3_cellranger710','ZY-307-001-10XSC3_cellranger710',
              'ZY-307-002-10XSC3_cellranger710','GSM4555888_cellranger710','GSM4555889_cellranger710',
              'GSM4555891_cellranger710','GSM6123277_cellranger710','CTT-307-CTC-10XSC3_cellranger710',
              'LHP-307-CTC-10XSC3_cellranger710','LXZ-307-CTC-10XSC3_cellranger710')
Combined.SeuratObj$Cell_type = ifelse(Combined.SeuratObj$celltype == 'MTC',
                                      Combined.SeuratObj$cell_Subtype,as.character(Combined.SeuratObj$celltype))

for(sample in SampleIDs){
  print(sample)
  out = get_barcode_cell_annotation(Combined.SeuratObj,sample)
}

Combined.SeuratObj = readRDS('Combined.SeuratObj')
SampleIDs = c(
  "LHP-307-CTC-10XSC3_cellranger710",
  "LHP-307-001-10XSC3_cellranger710",
  "LXZ-307-CTC-10XSC3_cellranger710",
  "LXZ-307-001-10XSC3_cellranger710",
  "230269A_JFP01_cellranger710",
  "230269A_JFP02_cellranger710",
  "ZY-307-001-10XSC3_cellranger710",
  "ZY-307-002-10XSC3_cellranger710",
  "CTT-307-CTC-10XSC3_cellranger710","GSM4259357_cellranger710",
              "GSM4555888_cellranger710","GSM4555889_cellranger710",
              "GSM4555891_cellranger710","GSM6123277_cellranger710",
              
              "MZH-307-001-10XSC3_cellranger710",
              "ZMY-307-001-10XSC3_cellranger710",
              
              "ZMY-307-002-10XSC3_cellranger710")

for(sample in SampleIDs){
  print(sample)
  out = numbat_analysis(Combined.SeuratObj, sample)
}
out.dir=paste0('result/OUT_Numbat_1/','MTC','/')

sce.obj = subset(Combined.SeuratObj, subset = celltype=='MTC')
# Extract count matrix from Seurat object
count_mat <- sce.obj@assays$RNA@counts


colnames(count_mat) <- paste0(sce.obj$SampleID,'_',sce.obj$CellID)
# Import allele counts data and select overlapping cells
dfs_allele = NULL
for (sID in SampleIDs) {
  input.file=paste0('result/OUT_Numbat/',sID,'/',sID,'_allele_counts.tsv.gz')
  df_allele <- read.table(input.file, header = TRUE, sep = "\t", quote = "", comment.char = "")
  df_allele$cell = paste0(sID,'_',df_allele$cell)
  dfs_allele = rbind(dfs_allele,df_allele)
}
ids <- intersect(dfs_allele$cell, colnames(count_mat)) %>% unique()

dfs_allele <- dfs_allele[which(dfs_allele$cell %in% c(ids)),] 
count_mat = count_mat[,ids]
dfs_allele=as.data.frame(dfs_allele)
# Run numBat
out <- run_numbat(count_mat, ref_hca, dfs_allele,
                  genome = "hg38", t = 1e-5, ncores = 32,
                  plot = T, out_dir = out.dir)

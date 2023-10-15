rm(list = ls())
gc(reset = T)
source('script/Functions.R')
options(stringsAsFactors = F)
ids.SampleType = c('Blood','Brain','CSF')

sce = readRDS('Combined.SeuratObj')
sce$Sample.Name = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',sce$SampleID)

plot.path='plts/'
CKmix = c('KRT7','KRT8','KRT18','KRT19')
CellMarkers = list(
  'MTC' = c('EPCAM','SCGB2A2','CD24','KRT7','KRT8','KRT18','KRT19'), # 
  'Neutrophils' = c('CSF3R','S100A8','S100A9'),
  'Macrophage' = c('AIF1','LYZ','CD163'),
  'Microglia' = c('TMEM119','CSF1R','ITGAM'), # Microglia
  "T Cell" = c("CD3D",'CD3E','CD2'),
  'vSMC' = c('RGS5','TAGLN','ACTA2'),
  'Endothelial' = c('CLDN5','PECAM1','RAMP2'), # Endothelial
  'Astrocytes' = c('GFAP','SLC1A2','AQP4'),
  'Oligodendrocytes'= c('MBP','MOBP','MOG') # Oligodendrocytes
)
meta = sce@meta.data
if(F){
  df.order = as.data.frame(table(sce$SampleID,sce$orig.ident))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  grouped_data <- meta %>%
    group_by(orig.ident, SampleID, celltype) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data$SampleID = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',grouped_data$SampleID)
  grouped_data$SampleID = factor(grouped_data$SampleID,levels = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',as.character(df.order$Var1)))
  x = ggplot(grouped_data, aes(x = SampleID, y = percentage, fill = celltype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell.type.cls) +
    labs(title = "Cell Type Distribution by Sample",
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal()
  
  df = table(meta$celltype,meta$orig.ident) %>% as.data.frame()
  df$Var2 = factor(df$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df$Var1 = factor(df$Var1, levels = c("Oligodendrocytes","Astrocytes","Endothelial","vSMC","Microglia","Neutrophils","T Cell","Macrophage","MTC"))
  
  # df = df[order(df$Freq),]
  # df$Var1 = factor(df$Var1, levels = df$Var1)
  y = ggplot(df,aes(Freq,Var1,fill=Var1)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cell.type.cls) +
    labs(x = 'Number of cells', 
         fill = "Cell Type") +scale_x_continuous(expand = c(0.01,0))+
    facet_wrap(~Var2,nrow = 1)+
    theme_minimal()
  y
  plt = ggarrange(x, y, ncol = 2, labels = c("a", "b"))
}
if(F){
  sce$UMAP1 = sce@reductions$umap@cell.embeddings[,1]
  sce$UMAP2 = sce@reductions$umap@cell.embeddings[,2]
  x = RidgePlot(sce,features = c('UMAP1','UMAP2'),sort=T,stack=T)+
    ggridges::stat_density_ridges(rel_min_height = 0.01,alpha=0.5)
  topptx(x,'plts/UMAP-Rideg.pptx',width = 12)
  
}

if(F){
  DimPlot(sce,raster=FALSE,label = F,group.by = 'Sample.Name') + 
    scale_color_manual(values = sample.cls)  + NoLegend()+
    theme(
      axis.title = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()) + ggtitle('')
  a = DimPlot(sce,label =F,group.by = 'Sample.Name')+theme(legend.position = 'bottom') +
    scale_color_manual(values = sample.cls)
}

if(F){
  ids =c("Primary.ER","Primary.PR","Primary.HER2","BM.ER","BM.PR","BM.HER2")
  for(i in ids){
    DimPlot(sce,raster=FALSE,label = F,group.by = i) + 
      scale_color_npg(alpha = 0.7) + NoLegend()+
      theme(
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        axis.line = element_blank()) + ggtitle('')
    a = DimPlot(sce,label =F,group.by = i)+
      scale_color_npg(alpha = 0.7)
  }
}
Combined.SeuratObj = readRDS('Combined.SeuratObj')
Combined.SeuratObj$Sample.Name = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',Combined.SeuratObj$SampleID)

sce = subset(Combined.SeuratObj,subset = orig.ident == "Adjacent Tumor")
if(T){
  sce = SCTransform(sce, method = "glmGamPoi", 
                    vars.to.regress = c("pMT"),seed.use = 777,
                    ncells=dim(sce)[2]*0.2,variable.features.n = 3500)
  sce = Seurat::RunPCA(sce)
  sce = FindNeighbors(sce, reduction = "pca", dims = 1:30)
  sce = FindClusters(sce,resolution = 1)
  sce = RunUMAP(sce, reduction = "pca", dims = 1:30,learning.rate = 0.5)
}
DimPlot(sce, group.by = 'celltype')  + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = cell.type.cls) + NoLegend()+
  theme(
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')
ggsave(paste0('plts/AT_UMAP.pdf'),height = 15,width =15,units = 'cm')

sce = subset(Combined.SeuratObj,subset = orig.ident == "Blood")
if(T){
  sce = SCTransform(sce, method = "glmGamPoi", 
                    vars.to.regress = c("pMT"),seed.use = 777,
                    ncells=dim(sce)[2]*0.2,variable.features.n = 3500)
  sce = Seurat::RunPCA(sce)
  sce = FindNeighbors(sce, reduction = "pca", dims = 1:30)
  sce = FindClusters(sce,resolution = 1)
  sce = RunUMAP(sce, reduction = "pca", dims = 1:30,learning.rate = 0.5)
}
DimPlot(sce, group.by = 'celltype')  + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = cell.type.cls) + NoLegend()+
  theme(
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')
ggsave(paste0('plts/BLOOD_UMAP.pdf'),height = 15,width =15,units = 'cm')

sce = subset(Combined.SeuratObj,subset = orig.ident == "CSF")
if(T){
  sce = SCTransform(sce, method = "glmGamPoi", 
                    vars.to.regress = c("pMT"),seed.use = 777,
                    ncells=dim(sce)[2]*0.2,variable.features.n = 3500)
  sce = Seurat::RunPCA(sce)
  sce = FindNeighbors(sce, reduction = "pca", dims = 1:30)
  sce = FindClusters(sce,resolution = 1)
  sce = RunUMAP(sce, reduction = "pca", dims = 1:30,learning.rate = 0.5)
}
DimPlot(sce, group.by = 'celltype')  + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = cell.type.cls) + NoLegend()+
  theme(
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')
ggsave(paste0('plts/CSF_UMAP.pdf'),height = 15,width =15,units = 'cm')

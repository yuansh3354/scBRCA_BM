
rm(list = ls())
gc(reset = T)
# packages
pkgs <- c('scRNAtoolVis','ClusterGVis','scrabble','monocle3','tidyverse')
lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
source('script/Functions.R')
options(stringsAsFactors = F)
# Combined.SeuratObj = readRDS('Combined.SeuratObj')
# Combined.SeuratObj$Sample.Name = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',Combined.SeuratObj$SampleID)
# meta = Combined.SeuratObj@meta.data
plot.path='plts/'
palette = c("Grays","Light Grays","Blues2","Blues3","Purples2","Purples3","Reds2","Reds3","Greens2")
immport.genes = read.csv('GeneList.csv')


if(F){
  sce = subset(Combined.SeuratObj,subset = celltype %in% c('Microglia','Macrophage'))
  sce = CreateSeuratObject(counts = sce@assays$SCT@counts)
  sce = AddMetaData(sce, metadata = Combined.SeuratObj@meta.data)
  sce = NormalizeData(object = sce,normalization.method =  "LogNormalize",  scale.factor = 1e4)
  sce = FindVariableFeatures(object = sce,selection.method = "vst", nfeatures = 2000)
  sce = ScaleData(object = sce)
  sce = Seurat::RunPCA(object = sce, do.print = FALSE)
  sce <- RunUMAP(sce, reduction = "pca", dims = 1:5)
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:5)
  sce <- FindClusters(sce,resolution = 0.3)
  DimPlot(sce,group.by = 'orig.ident')
  DimPlot(sce) 
  Markers = c("SPP1","SELL","PRKCB","P2RY12","NLRP3","LST1","GPX1","EGR1","C1QA","APOE")
  DotPlot(sce,features = (Markers),cols = c("#FFFFFF", "blue"),col.min = 0.5,cluster.idents = T) + RotatedAxis()+
    theme_bw()
  sce$cell_Subtype = NA
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(12),'PRKCB+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(5,9),'SELL+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(2,10),'LST1+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(1),'APOE+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(8),'GPX1+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(0),'EGR1+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(3),'P2RY12+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(4),'NLRP3+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(11),'C1QA+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(6,7),'SPP1+',sce$cell_Subtype)
  sce$cell_Subtype.1 = paste0(sce$celltype,' ',sce$cell_Subtype)
  table(sce$cell_Subtype.1)
  ids = names(which(table(sce$cell_Subtype.1) > 50))
  sce = sce[,which(sce$cell_Subtype.1 %in% ids)]
  ids = gsub('\\+','',unique(sce$cell_Subtype))
  ids = ids[order(ids,decreasing = T)]
  x1 = jjDotPlot(object = sce,
                gene = ids,
                xtree = F,
                ytree = F,
                id = 'cell_Subtype',
                rescale = T,
                rescale.min = -5,
                rescale.max = 5,
                point.geom = F,
                tile.geom = T)
  x1
  topptx(x1,'plts/Mph-Marker.pptx',height = 10,width = 10)
  saveRDS(sce,file = 'MPH.STAD')
}
sce.std = readRDS('MPH.STAD')
ids = gsub('\\+','',unique(sce.std$cell_Subtype))
ids = ids[order(ids,decreasing = T)]
sce = sce.std
if(T){
  df.order = as.data.frame(table(sce$SampleID,sce$orig.ident))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  meta = sce@meta.data
  grouped_data <- meta %>%
    group_by(orig.ident, SampleID, cell_Subtype.1) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data$SampleID = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',grouped_data$SampleID)
  grouped_data$SampleID = factor(grouped_data$SampleID,levels = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',as.character(df.order$Var1)))
  x1 = ggplot(grouped_data, aes(x = SampleID, y = percentage, fill = cell_Subtype.1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Mg.colors) +
    labs(title = "Cell Type Distribution by Sample",
         x = NULL, #
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_classic()
  x1
  grouped_data = grouped_data %>%
    group_by(SampleID) %>%
    filter(percentage == max(percentage))
  y1 = ggplot(grouped_data, aes(SampleID,percentage,fill=cell_Subtype.1)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Mg.colors)+
    theme_classic()
}
if(T){
  df.order = as.data.frame(table(sce$SampleID,sce$SampleType))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  meta = sce@meta.data
  grouped_data <- meta %>%
    group_by(SampleType,cell_Subtype.1) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  x2 = ggplot(grouped_data, aes(x =  SampleType, y = percentage, fill = cell_Subtype.1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Mg.colors) +
    labs(title = "Cell Type Distribution by Position",
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() 
  x2
}
if(T){
  sce$BC.Type = ifelse(sce$BM.HER2 == 'Pos','HER2',
                       ifelse((sce$BM.ER == 'Pos')|(sce$BM.PR == 'Pos'),'HR+',
                              ifelse(sce$BM.HER2 == 'Neg' & sce$BM.ER == 'Neg' & sce$BM.PR == 'Neg','TNBC',NA)))
  df.order = as.data.frame(table(sce$SampleID,sce$BC.Type))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  meta = sce@meta.data
  grouped_data <- meta %>%
    group_by(BC.Type,cell_Subtype.1) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data = na.omit(grouped_data)
  x3 = ggplot(grouped_data, aes(x =  BC.Type, y = percentage, fill = cell_Subtype.1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Mg.colors) +
    labs(title = "Cell Type Distribution by Sample",
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal()
  x3
  grouped_data = grouped_data %>%
    group_by(BC.Type) %>%
    filter(percentage == max(percentage))
  y3 = ggplot(grouped_data, aes(BC.Type,percentage,fill=cell_Subtype.1)) + 
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Mg.colors)+
    theme_classic() +
    theme(axis.text.x =element_text(angle = 45, hjust = 1),
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))+
    scale_y_continuous(limits = c(0,1),expand = c(0,0.01))
}
if(T){
  sce$BC.Type = ifelse(sce$BM.HER2 == 'Pos','HER2',
                       ifelse((sce$BM.ER == 'Pos')|(sce$BM.PR == 'Pos'),'HR+',
                              ifelse(sce$BM.HER2 == 'Neg' & sce$BM.ER == 'Neg' & sce$BM.PR == 'Neg','TNBC',NA)))
  df.order = as.data.frame(table(sce$SampleID,sce$Gender))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  meta = sce@meta.data
  grouped_data <- meta %>%
    group_by(Gender,cell_Subtype.1) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data = na.omit(grouped_data)
  x4 = ggplot(grouped_data, aes(x =  Gender, y = percentage, fill = cell_Subtype.1)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = Mg.colors) +
    labs(title = "Cell Type Distribution by Sample",
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() 
  
}
plt = ggarrange(x1,x2,x3,x4, ncol = 1, labels = c("a", "b",'c','d'))
plt
topptx(plt,file=paste0(plot.path,"Mph_CellSubtype_Freq_bar-in-datatype.pptx"),width = 10,height = 20)
meta = sce@meta.data
write.csv(meta,'Mph.SubCell.csv')
sce.std = readRDS('MPH.STAD')
sce.std$BC.Type = ifelse(sce.std$BM.HER2 == 'Pos','HER2',
                         ifelse((sce.std$BM.ER == 'Pos')|(sce.std$BM.PR == 'Pos'),'HR+',
                                ifelse(sce.std$BM.HER2 == 'Neg' & sce.std$BM.ER == 'Neg' & 
                                         sce.std$BM.PR == 'Neg','TNBC',NA)))
ids = gsub('\\+','',unique(sce.std$cell_Subtype))
ids = ids[order(ids,decreasing = T)]
if(F){
  sce = sce.std[c(immport.genes$Symbol,ids),]
  sce = NormalizeData(object = sce,normalization.method =  "LogNormalize",  scale.factor = 1e4)
  sce = FindVariableFeatures(object = sce,selection.method = "vst")
  sce = ScaleData(object = sce)
  sce = Seurat::RunPCA(object = sce, do.print = FALSE)
  sce <- RunUMAP(sce, reduction = "pca", dims = 1:5)
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:5)
  sce <- FindClusters(sce,resolution = 0.3)
  DimPlot(sce,group.by = 'cell_Subtype.1') + scale_color_manual(values = Mg.colors)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(face = "plain")) +
    labs(x=NULL,y=NULL,title=NULL) + NoLegend()+
    theme(
      axis.title = element_blank(),  
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()) + ggtitle('')

  DimPlot(sce,group.by = 'celltype') + scale_color_simpsons(alpha = 0.7)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(face = "plain")) +
    labs(x=NULL,y=NULL,title=NULL) + NoLegend()+
    theme(
      axis.title = element_blank(), 
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()) + ggtitle('')
  dt = sce@meta.data
  for (variable in ids) {
    dt[[variable]] = sce@assays$RNA@data[variable,]
    sce@meta.data[[variable]] = sce@assays$RNA@data[variable,]
  }
  library(ggpmisc)
  
  ggplot(dt, aes(APOE, SPP1)) +
    geom_point() +
    scale_color_manual(values = Mg.colors) +
    geom_smooth(method = "lm", se = T) +    
    stat_cor(method = "pearson",  label.x.npc = 0.02, label.y.npc = 0.98, vjust = 1, hjust = 0, )+
    theme(legend.position = "none")
  

  
  
  df <- dt %>%
    group_by(cell_Subtype.1) %>%
    summarise(
      "SPP1" =median(SPP1),
      "SELL"=median(SELL),
      "PRKCB"=median(PRKCB),
      "P2RY12"=median(P2RY12),
      "NLRP3"=median(NLRP3),
      "LST1"=median(LST1),
      "GPX1"=median(GPX1),
      "EGR1"=median(EGR1),
      "C1QA"=median(C1QA),
      "APOE"=median(APOE)) 
  library(ggradar)
  x1=ggradar(df,
             grid.min = min(df[,2:11]), 
             grid.mid = mymean(df[,2:11]), 
             grid.max = max(df[,2:11]), 
             group.colours = Mg.colors,
             group.point.size = 2,
             group.line.width = 1, 
             background.circle.transparency = 0, 
             legend.position = 'right', 
             legend.text.size = 12, 
             fill = TRUE, 
             fill.alpha = 0.3 
  ) + facet_wrap(~cell_Subtype.1)
  x1
  Idents(sce) = sce$cell_Subtype.1
  ids.sce = sce
  degs = FindAllMarkers(sce)
  degs$pct = degs$pct.1 - degs$pct.2
  degs$label = ifelse(degs$pct>0.3 & degs$avg_log2FC>1.5,'Positive',ifelse(degs$pct< -0.3 & degs$avg_log2FC< -2,'Negative','Other'))
  x = ggplot(degs, aes(pct, avg_log2FC)) + 
    geom_point(aes(color = label), position = position_jitter(width = 0.05, height = 0.05),size=0.2) + 
    facet_wrap(~ cluster) +
    theme_bw() +
    geom_smooth(method = 'lm') +
    scale_color_manual(values = c("#698EC3", "black", "#c10534"))
  x
  
  result = degs %>% 
    group_by(cluster) %>% 
    summarise(Count = n())
  result = result[order(result$Count),]
  result$cluster = factor(result$cluster,levels = result$cluster)
  
  x = ggplot(result, aes(x = cluster, y = Count,color=cluster,fill=cluster)) +
    geom_segment(aes(x = cluster, xend = cluster, y = 0, yend = Count)) +
    geom_point(size = 5, alpha = 0.7, shape = 21, stroke = 1) +
    theme_minimal_vgrid() +
    coord_flip() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_fill_manual(values = Mg.colors) +  
    scale_color_manual(values = Mg.colors)  +
    expand_limits(x = c(0, 0))+
    labs(x=NULL,y=NULL)
}

if(F){
  sce = sce.std
  Idents(sce) = sce$cell_Subtype.1
  degs.std = readRDS('degs.Microglia_degs.Macrophage')
  lapply(seq_along(unique(degs.std$cluster)), function(x){
    tmp <- degs.std |> 
      dplyr::filter(cluster == unique(degs.std$cluster)[x]) |> 
      dplyr::arrange(desc(avg_log2FC))
    
    # gene list
    glist <- tmp$avg_log2FC
    names(glist) <- tmp$gene
    
    # enrichment
    ego <- clusterProfiler::gseGO(geneList = glist,
                                  ont = "BP",
                                  OrgDb = org.Hs.eg.db,
                                  keyType = "SYMBOL",
                                  pvalueCutoff = 0.1)
    
    # to data.frame
    if (min(ego@result$NES[order(ego@result$NES,decreasing = T)][1:5])>0) {
      df <- data.frame(ego) |> 
        dplyr::arrange(pvalue,desc(NES))
    }else {
      df <- data.frame(ego) |> 
        dplyr::arrange(pvalue,NES)
    }
    
    # plot
    lapply(1:5, function(x){
      GseaVis::gseaNb(object = ego,
                      subPlot = 1,
                      geneSetID = df$ID[x],
                      addPval = T,
                      pvalX = 0.9,
                      pvalY = 0.8)
    }) -> plist
    
    # combine
    pplist <- cowplot::plot_grid(plotlist = plist,nrow = 1,align = 'hv')
    
    return(pplist)
  }) -> gglist
  names(gglist) <- paste("C",1:12,sep = "")
  
  for(i in 1:12){
    assign(paste0('x',i),gglist[[i]])
  }
  
  
  plt = ggarrange(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12, ncol = 1, labels = c("a", "b",'c','d','e','f','g','h','i','j','k','l'))
  plt
  
}

if(F){
  sce = subset(sce.std, subset = cell_Subtype.1 %in% c("Microglia SPP1+","Microglia EGR1+","Microglia P2RY12+",'Microglia NLRP3+','Macrophage NLRP3+'))
  counts <- as.data.frame(sce@assays$RNA@data) %>% 
    rownames_to_column(var = "Gene")
  
  meta <- data.frame(Cell=rownames(sce@meta.data), 
                     cell_type=sce@meta.data$cell_Subtype.1)
  
  fwrite(counts, "result/CellphoneDB/counts.txt", row.names=F, sep='\t')
  fwrite(meta, "result/CellphoneDB/meta.txt", row.names=F, sep='\t')
  pvals <- read.delim('result/CellphoneDB/Output/pvalues.txt', check.names = FALSE)
  means <- read.delim('result/CellphoneDB/Output/means.txt', check.names = FALSE)
  decon <- read.delim('result/CellphoneDB/Output/deconvoluted.txt', check.names = FALSE)
  library(ktplots)
  sce = subset(sce.std, subset = cell_Subtype.1 %in% c("Microglia SPP1+","Microglia EGR1+","Microglia P2RY12+",'Microglia NLRP3+','Macrophage NLRP3+'))
  sce$celltype = sce$cell_Subtype.1
  sce <- SingleCellExperiment(list(counts = sce@assays$RNA@counts), colData = sce@meta.data, rowData = rownames(sce))
  desiredInteractions = list(
    c("Macrophage NLRP3+", "Microglia P2RY12+"),
    c("Macrophage NLRP3+", "Microglia EGR1+"),
    c("Macrophage NLRP3+", "Microglia NLRP3+"),
    c("Macrophage NLRP3+", "Microglia SPP1+"))
  
  plot_cpdb2(
    scdata = sce,
    cell_type1 = "Macrophage NLRP3+",
    cell_type2 = ".",
    celltype_key = "celltype", # column name where the cell ids are located in the metadata
    means = means,
    pvals = pvals,
    deconvoluted = decon, # new options from here on specific to plot_cpdb2
    desiredInteractions = desiredInteractions
  )
  
  cellphoneDB_Dotplot(pvals.data = pvals,means.data = means,
                      target.cells_1 =grep('Micr',unique(sce.std$cell_Subtype.1),value = T),
                      target.cells_2 =c("Macrophage APOE+","Macrophage NLRP3+","Macrophage LST1+","Macrophage SPP1+") ,
                      gene_a = grep('^HLA|^COL',rownames(sce.std),value = T))
}

library(CellChat)
sce = sce.std
sce$cell_type = sce$cell_Subtype.1
cellchat <-createCellChat(object=sce,group.by ="cell_type")
groupSize <-as.numeric(table(cellchat@idents))  
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB<-CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)cellchat <- projectData(cellchat,PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 300)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")
cellchat <- computeCommunProbPathway(cellchat)
for(x in seq_along(unique(df.net$pathway_name))){
  plt = netVisual_chord_cell(cellchat, sources.use =c( 'Macrophage SPP1+',
   "Macrophage NLRP3+","Macrophage LST1+","Macrophage PRKCB+","Macrophage SELL+","Macrophage APOE+"),
                       targets.use = c("Microglia SPP1+","Microglia EGR1+","Microglia P2RY12+"),
                       signaling = unique(df.net$pathway_name)[x])
  plt
  topptx(plt,paste0('plts/cellchat/',unique(df.net$pathway_name)[x],'.pptx'),width = 7,height = 7)
}

pplist <- cowplot::plot_grid(plotlist = plist,nrow = 6,align = 'hv')
topptx(pplist,'plts/cellchat.pptx',width = 40,height = 40)

# table(sce$celltype, sce$orig.ident)
sce.mph = subset(sce.std, subset = celltype == 'Macrophage')
sce.mph$type = ifelse(sce.mph$orig.ident %in% c('Adjacent Tumor','Tumor'),'Resident','Circulating')
Idents(sce.mph) = sce.mph$type
degs.mph = FindMarkers(sce.mph,ident.1 = 'Resident',ident.2 = 'Circulating')
degs.mph = degs.mph[order(degs.mph$avg_log2FC,decreasing = T),]
glist <- degs.mph$avg_log2FC
names(glist) <- rownames(degs.mph)
ego.mph = clusterProfiler::gseGO(geneList = glist,
                       ont = "BP",
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       pvalueCutoff = 0.1)
ego.df = ego.mph@result
ego.df = ego.df[order(ego.df$NES,decreasing = T),]
df = rbind(head(ego.df,n=10),tail(ego.df,n=10))
df = df[order(df$NES),]
df$Description = factor(df$Description, levels = df$Description)
x1 = ggplot(df,aes(NES,Description)) + 
  geom_bar(stat = 'identity') +
  labs(y=NULL,x='Normalized Enrichment Score')+
  scale_x_continuous(limits = c(-4.51,2.6),expand = c(0,0.01))+
  theme_bw()
x1

degs = FindAllMarkers(sce.mph,verbose = T,only.pos = T)
ids.degs = degs[which(abs(degs$pct.1 - degs$pct.2)>0.3),]
dim(ids.degs)
table(ids.degs$cluster)
Resident = ids.degs[which(ids.degs$cluster == 'Resident'),]
Resident = Resident[order(Resident$avg_log2FC,decreasing = T),]
Resident$gene[1:10]
Circulating = ids.degs[which(ids.degs$cluster == 'Circulating'),]
Circulating = Circulating[order(Circulating$avg_log2FC,decreasing = T),]
Circulating$gene[11:20]

Markers = c('CD274','LAG3','IDO1','FOXP3','CD276','TIGIT','CD96','CDKN2D',"SOCS1",'TGFB1')
x1 = jjDotPlot(object = sce.mph,
          gene = Markers,
          xtree = F,
          ytree = F,
          id = 'type',
          rescale = T,
          rescale.min = -5,
          rescale.max = 5,
          point.geom = F,
          tile.geom = T)
x1
DimPlot(sce.mph)+ NoLegend()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')


df = sce.mph@meta.data
result = df %>% group_by(type,cell_Subtype.1) %>% 
  summarise(Count = n())
df <- result %>%
  group_by(type) %>%
  mutate(Percent = Count / sum(Count) * 100)

pie_plot <- ggplot(df, aes(x = "", y = Percent, fill = cell_Subtype.1)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  labs(x = NULL, y = NULL, fill = "cell_Subtype.1") +
  theme_void() +
  theme(legend.position = "right")+
  scale_fill_manual(values = Mg.colors)

plot1 <- pie_plot +
  facet_wrap(~ type, nrow = 1) +
  geom_text(aes(label = ifelse(Percent > 15, paste0(round(Percent, 1), "%"), "")),
            position = position_stack(vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5))
plot2 <- pie_plot +
  facet_wrap(~ type, nrow = 1) +
  geom_text(aes(label = ifelse(Percent > 15, paste0(round(Percent, 1), "%"), "")),
            position = position_stack(vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5))

print(plot2)
Idents(sce.mph) = sce.mph$orig.ident
degs1 = FindMarkers(sce.mph, ident.1 = 'Blood',ident.2 = 'Tumor',logfc.threshold = 0)
degs2 = FindMarkers(sce.mph, ident.1 = 'CSF',ident.2 = 'Tumor',logfc.threshold = 0)
ids = intersect(rownames(degs1),rownames(degs2))
degs1 = degs1[ids,]
degs2 = degs2[ids,]
degs = cbind(degs1,degs2)
colnames(degs) = c(paste0('V',seq(1:10)))
degs$label = ifelse(degs$V2< -0.25 & degs$V7< -0.25,'Down in Circulating',
                    ifelse(degs$V2>0.25 & degs$V7>0.25,'Up in Circulating','#Ng'))
degs$gene = rownames(degs)
degs = degs[!grepl("^RP[SL][[:digit:]]",rownames(degs)),]
top_right <- head(degs[order(degs$V2 + degs$V7, decreasing = TRUE), "target"], 10)
bottom_left <- head(degs[order(degs$V2 + degs$V7, decreasing = FALSE), "target"], 10)

top_right_pos <- degs[degs$target %in% top_right, c("V2", "V7")]
bottom_left_pos <- degs[degs$target %in% bottom_left, c("V2", "V7")]

x = ggplot(degs, aes(V2, V7)) +
  theme_bw() +
  geom_point(aes(color = label)) +
  scale_color_manual(values = c('Down in Circulating' = '#4798CD',
                                'Up in Circulating' = '#a80c11',
                                'Ng' = '7f7f7f')) +
  geom_smooth(method = "lm", span = 2) +
  stat_cor(method = "spearman", label.x = -6, label.y = 4) +
  labs(x = 'Differential Blood vs MTC (log2FC)',
       y = 'Differential CSF vs MTC (log2FC)') +
  xlim(c(-7, 7)) +
  ylim(c(-4, 4)) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "solid", color = "gray") +
  geom_hline(yintercept = c(-0.25, 0.25), linetype = "solid", color = "gray") +
  theme(legend.position = c(0.75, 0.95),
        legend.justification = c(0.5, 1),
        legend.direction = "horizontal")+
  geom_text_repel(data = bottom_left_pos, aes(label = rownames(bottom_left_pos)), vjust = -1) +
  geom_text_repel(data = top_right_pos, aes(label = rownames(top_right_pos)), vjust = 1)
sce.std = readRDS('MPH.STAD')
if(T){
  sce = subset(sce.std, subset = celltype == 'Macrophage')
  sce = subset(sce, subset = Sample.Name %in%  c("230269A_JFP01","230269A_JFP02","CTT-307-CTC","GSM4259357","GSM4555888",
                                                 "GSM4555889","GSM4555891","GSM6123277","LHP-307-001","LXZ-307-001","LXZ-307-CTC",
                                                 "MZH-307-001","ZMY-307-001","ZY-307-001","ZY-307-002","LHP-307-CTC","ZMY-307-002"))
  
  ids = names(which(table(sce$SampleID) > 50))
  sce = sce[,which(sce$SampleID %in% ids)]
  sce = CreateSeuratObject(counts = sce@assays$RNA@counts,min.cells = 500)
  sce = AddMetaData(sce, metadata = sce.std@meta.data)
  sce = NormalizeData(object = sce,normalization.method =  "LogNormalize",  scale.factor = 1e4)
  sce = FindVariableFeatures(object = sce,selection.method = "vst")
  sce = ScaleData(object = sce)
  sce = Seurat::RunPCA(object = sce, do.print = FALSE)
  sce <- RunUMAP(sce, reduction = "pca", dims = 1:5)
  DimPlot(sce,group.by = 'orig.ident')
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:5)
  sce <- FindClusters(sce,resolution = 0.3)
}

expression_matrix = sce@assays$RNA@data
cell_metadata = data.frame(sce@meta.data)
gene_annotation = data.frame(expression_matrix[,1])
gene_annotation[,1] = row.names(gene_annotation)
colnames(gene_annotation)=c("gene_short_name")
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


cds <- preprocess_cds(cds, num_dim = 10)
cds <- align_cds(cds, alignment_group = "orig.ident")

cds <- reduce_dimension(cds,cores=16)
cds.ids <- reduce_dimension(cds,cores=16)

cds <- cluster_cells(cds,resolution = 0.0000001)
cds <- learn_graph(cds)

myselect <- function(cds,select.classify,my_select){
  cell_ids <- which(colData(cds)[,select.classify] == my_select)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
  }

cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'orig.ident',my_select = "Blood"))

plot_cells(cds, 
           color_cells_by = "pseudotime",cell_size = 0.7,
           label_branch_points=FALSE,show_trajectory_graph=F)+theme_classic()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + NoLegend()
  plot_cells(cds,
             color_cells_by = "cell_Subtype.1",
             label_cell_groups=FALSE,cell_size = 0.7,
             label_leaves=FALSE,show_trajectory_graph=F,
             label_branch_points=FALSE)+theme_classic() +scale_color_manual(values = Mg.colors)+
    theme(
      axis.title = element_blank(),  
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()) + NoLegend()

Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=16)
cds = readRDS('cds.mph.RDS')
Track_genes = readRDS('Track_genes.mph.RDS')
ggplot(Track_genes,aes(morans_I))


ids = colnames(sce) %>% 
  gsub('_\\d+','',.)
ids = paste0(ids,'_',sce$SampleID)
sce = Seurat::RenameCells(sce,new.names = ids)
pseudo.df = as.data.frame(pseudotime(cds)) %>% set_colnames('pseudotime')
pseudo.df$pseudotime.raw = pseudo.df$pseudotime
pseudo.df$pseudotime = round(pseudo.df$pseudotime)
# sce = sce.std
sce = AddMetaData(sce.std, metadata = pseudo.df)
df = sce@meta.data
df_grouped <- df %>%
  group_by(pseudotime, cell_Subtype.1) %>%
  summarise(Count = n()) %>%
  ungroup()

df_max_percent <- df_grouped %>%
  group_by(pseudotime) %>%
  mutate(Percent = Count / sum(Count)) %>% na.omit()
x = ggplot(df_max_percent, aes(pseudotime,Percent,color=cell_Subtype.1 )) + 
  geom_line() + 
  scale_color_manual(values = Mg.colors)
x
df_max_percent.ids <- df_grouped %>%
  group_by(pseudotime) %>%
  mutate(Percent = Count / sum(Count)) %>%
  filter(Percent == max(Percent)) %>%
  ungroup() %>% na.omit()
gene.list = Track_genes[which(Track_genes$p_value ==0),'gene_short_name']
write.csv(df_max_percent.ids,'df_max_percent.ids.csv')
celltype='Mph'

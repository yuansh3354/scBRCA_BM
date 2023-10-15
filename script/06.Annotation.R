srm(list = ls())
gc(reset = T)
# packages
pkgs <- c('scRNAtoolVis','ClusterGVis','scrabble')
lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
source('script/Functions.R')
options(stringsAsFactors = F)
Combined.SeuratObj = readRDS('Combined.SeuratObj')
Combined.SeuratObj$Sample.Name = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',Combined.SeuratObj$SampleID)
meta = Combined.SeuratObj@meta.data
plot.path='plts/'
palette = c("Grays","Light Grays","Blues2","Blues3","Purples2","Purples3","Reds2","Reds3","Greens2")

sce = subset(Combined.SeuratObj,idents = c('MTC'))
sce = CreateSeuratObject(counts = sce@assays$RNA@counts)
sce = AddMetaData(sce, metadata = Combined.SeuratObj@meta.data)
if(T){
  sce <- SCTransform(sce, method = "glmGamPoi", vars.to.regress = "nFeature_RNA", verbose = FALSE)
  sce = Seurat::RunPCA(object = sce, do.print = FALSE)
  sce <- RunUMAP(sce, reduction = "pca", dims = 1:30,learning.rate = 0.5)
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:30)
  sce <- FindClusters(sce,resolution = 2)
  DimPlot(sce,group.by = 'Sample.Name')  + labs(title = NULL)+
    scale_color_manual(values = sample.cls) + NoLegend()+
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()) + ggtitle('')
  ggsave('plts/MTC-SampleID.pdf',height = 15,width = 15,units = 'cm')
  Markers =  list(  
    "HER2 Family" = c('ERBB2'),
    "AR Positive" =c('AR'),
    "Hormone Receptor" = c('PRLR','ESR1'), 
    #'NOTCH Related' = c('NOTCH1','NOTCH2','NOTCH3','JAG1'),
    "BM Related" = c('AXL','ITGB1','ST6GALNAC5'),#'YTHDF3','GAS6','FN1',
    "Neuro Related" = c('NLGN3','SNCA','SYT11','GPHN')
  )


  DotPlot(sce,features = (Markers),cols = c("#FFFFFF", "blue"),col.min = 1,cluster.idents=T) + RotatedAxis()
  meta = sce@meta.data
  # [,which(is.na(sce$cell_Subtype))]
  sce$cell_Subtype = 'MTC blank'
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(47,55,51,56,12,16),'Neuro Related',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(39,61,23,33,19,11,17),'AR+HR+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(0,21,22,25,54,58,49,44,29,1,20),'HR+',sce$cell_Subtype)
  sce$cell_Subtype = ifelse(Idents(sce) %in% c(11,15,40,60,62,42,59,37,28,34,53,8,24,41,13,32,4,2,6,38,9,14,52,50),'HER2+',sce$cell_Subtype)
  
  DimPlot(sce,group.by = 'cell_Subtype') + scale_color_manual(values = MTC.cls) + NoLegend()+
    theme(
      axis.title = element_blank(),  
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()) + ggtitle('')
  ggsave('plts/MTC-sub.pdf',width = 8,height = 8)
  sce$cell_Subtype = factor(sce$cell_Subtype, 
                            levels = c('MTC blank','Neuro Related',
                                       'HR+','AR+HR+','HER2+' ))
  x = jjDotPlot(object = sce,
                gene = unlist(Markers),
                xtree = F,
                ytree = F,
                id = 'cell_Subtype',
                rescale = T,
                rescale.min = -5,
                rescale.max = 5,
                point.geom = F,
                tile.geom = T)
  x
  topptx(x,'plts/MTC-Marker.pptx')
  # saveRDS(sce,file = 'MTC.STAD')
}
sce.std = sce
sce = sce.std
MTC = c('EPCAM','CD44','CD24','KRT7','KRT8','KRT18','KRT19')
a = DotPlot(sce,features = MTC,group.by = 'orig.ident')
df = a$data
x = ggplot(df,aes(id,pct.exp,alpha=log(avg.exp))) + theme_minimal() +
  labs(title = NULL,
       y = "Percentage Expression (%)",alpha = "Relative Expr")+
  geom_bar(stat = 'identity',fill='#1f77b4') + facet_wrap(~features.plot,ncol = 1)
x

if(T){
  df.order = as.data.frame(table(sce$SampleID,sce$orig.ident))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  meta = sce@meta.data
  grouped_data <- meta %>%
    group_by(orig.ident, SampleID, cell_Subtype) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data$SampleID = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',grouped_data$SampleID)
  grouped_data$SampleID = factor(grouped_data$SampleID,levels = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',as.character(df.order$Var1)))
  x1 = ggplot(grouped_data, aes(x = SampleID, y = percentage, fill = cell_Subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cls) +
    labs(title = "Cell Type Distribution by Sample",
         x = NULL, # 
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x =element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_blank(),  
          axis.title.x = element_blank(),  
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  x1
}
if(T){
  df.order = as.data.frame(table(sce$SampleID,sce$SampleType))
  df.order$Freq = ifelse(df.order$Freq == 0, NA, df.order$Freq)
  df.order = na.omit(df.order)
  df.order$Var2 = factor(df.order$Var2,levels = c('Tumor','Adjacent Tumor', 'Blood','CSF'))
  df.order = df.order[order(df.order$Var2),]
  meta = sce@meta.data
  grouped_data <- meta %>%
    group_by(SampleType,cell_Subtype) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  x2 = ggplot(grouped_data, aes(x =  SampleType, y = percentage, fill = cell_Subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cls) +
    labs(title = "Cell Type Distribution by Position",
         x = NULL, 
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x =element_text(angle = 45, hjust = 1), # 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(),
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
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
    group_by(BC.Type,cell_Subtype) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data = na.omit(grouped_data)
  x3 = ggplot(grouped_data, aes(x =  BC.Type, y = percentage, fill = cell_Subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cls) +
    labs(title = "Cell Type Distribution by Sample",
         x = NULL, # 
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x =element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  x3
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
    group_by(Gender,cell_Subtype) %>%
    summarise(count = n()) %>%
    mutate(total_count = sum(count), percentage = count / total_count)
  grouped_data = na.omit(grouped_data)
  x4 = ggplot(grouped_data, aes(x =  Gender, y = percentage, fill = cell_Subtype)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cls) +
    labs(title = "Cell Type Distribution by Sample",
         x = NULL, 
         y = "Percentage (%)",
         fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x =element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_blank(), 
          axis.title.x = element_blank(), 
          plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))
  x4
}
DimPlot(sce,group.by = 'cell_Subtype') + scale_color_manual(values = MTC.cls)
ggsave('plts/MTC-sub.pdf',width = 10,height = 8)
plt = ggarrange(x1,x2,x3,x4, ncol = 1, labels = c("a", "b",'c','d'))
plt
topptx(plt,file=paste0(plot.path,"MTC_CellSubtype_Freq_bar-in-datatype.pptx"),width = 10,height = 10)
meta = sce@meta.data
write.csv(meta,'MTC.SubCell.csv')
plt = ggarrange(x1,x2,x3,x4, 
                ncol = 1, 
                labels = c("a", "b",'c','d'))
plt

meta = sce@meta.data
grouped_data <- meta %>%
  group_by(SampleType,cell_Subtype) %>%
  summarise(count = n()) %>%
  mutate(total_count = sum(count), 
         percentage = count / total_count,
         mean_percentage = mean(percentage),
         sd_percentage = sd(percentage))

error_df <- grouped_data %>%
  group_by(SampleType, cell_Subtype) %>%
  summarise(mean_percentage = mean(percentage),
            sd_percentage = sd(percentage),
            n = n())

plot <- ggplot(grouped_data, aes(x = SampleType, 
                                 y = percentage, fill = cell_Subtype)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal_hgrid() + 
  facet_wrap(~cell_Subtype,ncol=1) + 
  scale_fill_manual(values = MTC.cls)
plot
t.test(grouped_data[which(grouped_data$orig.ident %in%
                            c('Adjacent Tumor','Tumor')),'percentage'],
       grouped_data[which(grouped_data$orig.ident %in% c('Blood')),
                    'percentage'])
sce.std = readRDS('MTC.STAD')
# sce.int = readRDS('MTC.integrated')
sce.std$BC.Type = ifelse(sce.std$BM.HER2 == 'Pos','HER2',
                     ifelse((sce.std$BM.ER == 'Pos')|(sce.std$BM.PR == 'Pos'),'HR+',
                            ifelse(sce.std$BM.HER2 == 'Neg' & sce.std$BM.ER == 'Neg' & 
                                     sce.std$BM.PR == 'Neg','TNBC',NA)))

Idents(sce.std) = sce.std$cell_Subtype
# degs.std = FindAllMarkers(sce.std)
# Idents(sce.int) = sce.int$cell_Subtype
# degs.int = FindAllMarkers(sce.int)
# saveRDS(degs.std,file = 'degs.MTC.std')
degs.std = readRDS('degs.MTC.std')
degs.std.markers <- degs.std %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_log2FC)

x = split(degs.std$gene,degs.std$cluster)

x=jjVolcano(diffData = degs.std,
          log2FC.cutoff = 1,
          size = 3.5,
          topGeneN = 10,
          tile.col = MTC.cls)
topptx(x,filename = 'plts/jjVolcano.pptx')

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
  df <- data.frame(ego) |> 
    dplyr::arrange(pvalue,desc(NES))
  # df = df[which(df$NES > 0),]
  
  # plot
  lapply(c(2,3,5,7,20), function(x){
    GseaVis::gseaNb(object = ego,
                    subPlot = 2,
                    geneSetID = df$ID[x],
                    addPval = T,
                    pvalX = 0.9,
                    pvalY = 0.8)
  }) -> plist
  
  # combine
  pplist <- cowplot::plot_grid(plotlist = plist,nrow = 1,align = 'hv')
  
  return(pplist)
}) -> gglist

names(gglist) <- paste("C",1:5,sep = "")
sce.std = NormalizeData(object = sce.std,assay='RNA',normalization.method =  "LogNormalize",  scale.factor = 1e6)
sce.std = FindVariableFeatures(object = sce.std,assay='RNA',selection.method = "vst", nfeatures = 3000)
sce.std = ScaleData(object = sce.std,assay='RNA')

st.data <- prepareDataFromscRNA(object = sce.std,
                               diffData = degs.std.markers,
                               showAverage = TRUE)
x = visCluster(object = st.data,
           plot.type = "both",
           line.side = "left",
           column_names_rot = 45,
           markGenes = degs.std.markers$gene)
x
table(degs.std[which(degs.std$p_val_adj<0.05 & abs(degs.std$avg_log2FC) >1),]$cluster)
degs.std.selected = degs.std[which(degs.std$p_val_adj<0.05),]

x1 = gglist[[1]]
x2 = gglist[[2]]
x3 = gglist[[3]]
x4 = gglist[[4]]
x5 = gglist[[5]]
plt = ggarrange(x1,x2,x3,x4,x5, ncol = 1, labels = c("a", "b",'c','d','e'))
plt
get_module = function(degs,cell.name){
  module = degs[which(degs$cluster %in% cell.name),]
  module = module[order(module$avg_log2FC,decreasing = T),]
  module = module[1:100,]
  module = module$gene
  return(module)
}
module.blank = get_module(degs.std.selected,'MTC blank')
module.neuro = get_module(degs.std.selected,'Neuro Related')
module.hr = get_module(degs.std.selected,c('HR+','AR+HR+'))
module.her2 = get_module(degs.std.selected,'HER2+')

modules = list(
  module.blank=module.blank,
  module.neuro=module.neuro,
  module.hr=module.hr,
  module.her2=module.her2
)

sce.std = AddModuleScore(sce.std,features = modules)
df = scrabble::hierarchy(sce.std@meta.data[,c("Cluster1","Cluster2","Cluster3","Cluster4")])
df$cell_Subtype = sce.std$cell_Subtype
groups = rownames(df)
names(groups) = df$cell_Subtype %>% as.character()

mynames = c("MTC blank" ,'Neuro Related' ,"HR+",    "HER2+" )
scrabble::plot_hierarchy(df,legend.pos = "top",
                         quadrant.names = mynames,
                         groups=groups,group.cols = MTC.cls,
                         legend = F)

scrabble::plot_hierarchy(df,
                         quadrant.names =mynames,
                         groups=groups[groups %in% rownames(df[which(df$cell_Subtype == 'MTC blank'),])],
                         group.cols = MTC.cls,legend = F)scrabble::plot_hierarchy(df,
                         quadrant.names = mynames,
                         groups=groups[groups %in% rownames(df[which(df$cell_Subtype == 'HR+'),])],
                         group.cols = MTC.cls,legend = F)
scrabble::plot_hierarchy(df,
                         quadrant.names =  mynames,
                         groups=groups[groups %in% rownames(df[which(df$cell_Subtype == 'AR+HR+'),])],
                         group.cols = MTC.cls,legend = F)
scrabble::plot_hierarchy(df,
                         quadrant.names =  mynames,
                         groups=groups[groups %in% rownames(df[which(df$cell_Subtype == 'HER2+'),])],
                         group.cols = MTC.cls,legend = F)
scrabble::plot_hierarchy(df,
                         quadrant.names =  mynames,
                         groups=groups[groups %in% rownames(df[which(df$cell_Subtype == 'Neuro Related'),])],
                         group.cols = MTC.cls,legend = F)

mmc = read.csv('MMC.csv')
modules = list(
 p1 = mmc[,1],
 p2 = mmc[,2],
 p3 = mmc[,3],
 p4 = mmc[,4],
 p5 = mmc[,5],
 p6 = mmc[,6],
 p7 = mmc[,7],
)
mymean = function(x){mean(apply(x,1,mean))}
sce.std = AddModuleScore(sce.std,features = modules,ctrl = 30)
dt = sce.std@meta.data
df <- dt %>%
  group_by(orig.ident, cell_Subtype) %>%
  summarise("DNA Replication" = mean(Cluster1),
            "G2/M" = mean(Cluster2),
            "Stress Response" = mean(Cluster3),
            "Developmental Growth" = mean(Cluster4),
            "pre-mRNA Matuation / core Spliceosome" = mean(Cluster5),
            "ECM deposition/ EMT" = mean(Cluster6),
            "Inflammatory Response" = mean(Cluster7)
            ) 
library(ggradar)
df = df[,-5]
ids = which(df$orig.ident == 'Adjacent Tumor')
x1=ggradar(df[ids,2:9],
           grid.min = min(df[ids,3:9]), 
           grid.mid = mymean(df[ids,3:9]), 
           grid.max = max(df[ids,3:9]), 
           group.colours = MTC.cls,
           group.point.size = 4,
           group.line.width = 1, 
           background.circle.colour = 'grey', 
           background.circle.transparency = 0, 
           legend.position = 'right', 
           legend.text.size = 12, 
           fill = TRUE, 
           fill.alpha = 0.3 
)
ids = which(df$orig.ident == 'Tumor')
x2=ggradar(df[ids,2:9],
           grid.min = min(df[ids,3:9]), 
           grid.mid = mymean(df[ids,3:9]), 
           grid.max = max(df[ids,3:9]), 
           group.colours = MTC.cls,
           group.point.size = 4,
           group.line.width = 1, 
           background.circle.colour = 'grey', 
           background.circle.transparency = 0, 
           legend.position = 'right', 
           legend.text.size = 12, 
           fill = TRUE, 
           fill.alpha = 0.3 
)
ids = which(df$orig.ident == 'CSF')
x3=ggradar(df[ids,2:9],
           grid.min = min(df[ids,3:9]), 
           grid.mid = mymean(df[ids,3:9]), 
           grid.max = max(df[ids,3:9]), 
           group.colours = MTC.cls,
           group.point.size = 4,
           group.line.width = 1,
           background.circle.colour = 'grey',
           background.circle.transparency = 0, 
           legend.text.size = 12, 
           fill = TRUE, 
           fill.alpha = 0.3  
)
ids = which(df$orig.ident == 'Blood')
x4=ggradar(df[ids,2:9],
           grid.min = min(df[ids,3:9]),
           grid.mid = mymean(df[ids,3:9]), 
           grid.max = max(df[ids,3:9]), 
           group.colours = MTC.cls,
           group.point.size = 4,
           group.line.width = 1, 
           background.circle.colour = 'grey', 
           background.circle.transparency = 0,
           legend.position = 'right', 
           legend.text.size = 12,            
           fill = TRUE, 
           fill.alpha = 0.3 
)
plt = ggarrange(x1,x2,x3,x4, ncol = 4, labels = c("a", "b",'c','d'))
plt

df <- dt %>%
  group_by(CellID,orig.ident) %>%
  summarise("DNA Replication" = mean(Cluster1),
            "G2/M" = mean(Cluster2),
            "Stress Response" = mean(Cluster3),
            "Developmental Growth" = mean(Cluster4),
            "pre-mRNA Matuation / core Spliceosome" = mean(Cluster5),
            "ECM deposition/ EMT" = mean(Cluster6),
            "Inflammatory Response" = mean(Cluster7)
  ) %>% as.data.frame()

library(Mfuzz)
dt <- new("ExpressionSet",exprs = df)




df = reshape2::melt(df,
          id.var=c('CellID','orig.ident'))
df$orig.ident = factor(df$orig.ident, levels = c("Tumor","Blood","CSF","Adjacent Tumor"))
ggplot(df,aes(orig.ident,value)) + geom_boxplot(color ='#C4425E' )+
  facet_wrap(~variable,nrow = 4) +
  theme_minimal() +
  theme(axis.text.x =element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"))+
  geom_signif(comparisons = list(c("Tumor", "Blood")), y_position = 1.2, tips_length = 0.02)


 
sce.std = readRDS('MTC.STAD')
# sce.int = readRDS('MTC.integrated')
sce.std$BC.Type = ifelse(sce.std$BM.HER2 == 'Pos','HER2',
                         ifelse((sce.std$BM.ER == 'Pos')|(sce.std$BM.PR == 'Pos'),'HR+',
                                ifelse(sce.std$BM.HER2 == 'Neg' & sce.std$BM.ER == 'Neg' & 
                                         sce.std$BM.PR == 'Neg','TNBC',NA)))
DimPlot(sce,group.by = 'cell_Subtype')

sce = CreateSeuratObject(counts = sce.std@assays$SCT@counts,
                         min.cells = 55901*0.85)

genes = rownames(sce) %>% as.data.frame()
genes = genes[!grepl('^MT',genes[,1]),,drop=F]
genes = genes[!grepl('^RP[SL][[:digit:]]',genes[,1]),,drop=F]
# gtf <- rtracklayer::import('refdata-gex-GRCh38-2020-A/genes/genes.gtf')
# gtf<-as.data.frame(gtf)
gtf = gtf[which(gtf$gene_name %in% genes[,1]),]
df = gtf[,c('seqnames','gene_name','transcript_type')]
df = unique(df)
df = df[which(df$transcript_type == 'protein_coding'),]
df$seqnames = as.character(df$seqnames)
ids = table(df$seqnames) %>% as.data.frame()
ids = ids[order(ids$Freq,decreasing = T),]
ids
ggplot(df,aes(seqnames)) + geom_bar()











rm(list = ls())
gc(reset = T)
pkgs = c('scRNAtoolVis','rtracklayer')
lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
source('script/Functions.R')
options(stringsAsFactors = F)
file_path <- "refdata-gex-GRCh38-2020-A/genes/lncRNA_data.csv"
lncRNA_data = read.csv(file_path)
lncRNA = lncRNA_data$gene_name
########################## Step.01 ########################## 
seurat.info = 'Step.01.read.seurat.log.txt'
if(file.exists(seurat.info)){file.remove(seurat.info)}
sink(seurat.info, append = TRUE)
for (ids in ids.SampleType) {
  assign(ids,dir(ids,pattern = '_cellranger710'))
  files = get(ids)
  for(file in files){
    file.path = paste(ids,file,'outs/filtered_feature_bc_matrix',sep = '/')
    temp = CreateSeuratObject(Read10X(file.path),min.cells = 100,min.features = 200)
    temp$SampleID = file
    temp$SampleType = ids
    temp$CellID = rownames(temp@meta.data)
    temp = eliminate_genes(temp)
    print(file)
    print(temp)
    cat("\n")
    assign(file,temp)
  }
}
sink()
rm(temp)

ls(pattern = "_cellranger710$")
sce.merge.obj = merge(`230269A_JFP01_cellranger710`,
                      y=c(
                        `230269A_JFP02_cellranger710`,`CTT-307-CTC-10XSC3_cellranger710`,
                        `GSM4259357_cellranger710`,`GSM4555888_cellranger710`,`GSM4555889_cellranger710`,
                        `GSM4555891_cellranger710`,`GSM5645891_cellranger710`,`GSM5645892_cellranger710`,
                        `GSM5645893_cellranger710`,`GSM6123277_cellranger710`,`GSM7475324_cellranger710`,
                        `GSM7475325_cellranger710`,`GSM7475326_cellranger710`,`LHP-307-001-10XSC3_cellranger710`,
                        `LXZ-307-001-10XSC3_cellranger710`,`LXZ-307-CTC-10XSC3_cellranger710`,
                        `MZH-307-001-10XSC3_cellranger710`,`ZMY-307-001-10XSC3_cellranger710`,
                        `ZY-307-001-10XSC3_cellranger710`,`ZY-307-002-10XSC3_cellranger710`,
                        `LHP-307-CTC-10XSC3_cellranger710`,`ZMY-307-002-10XSC3_cellranger710`
                        )
                      )
sce.merge.obj$`Gender` = ifelse(grepl("(ZY-307)", sce.merge.obj$SampleID),'M',
                                      ifelse(grepl("(230269A_)|(CTT-307)|(LHP-307)|(LXZ-307)|(MZH-307)|(ZMY-307)|(GSM5645)|(GSM74753)",
                                                   sce.merge.obj$SampleID),'F',NA))
sce.merge.obj$orig.ident = sce.merge.obj$SampleType
sce.merge.obj$orig.ident = ifelse(sce.merge.obj$SampleID %in% c("230269A_JFP02_cellranger710",
                                                                  "ZY-307-002-10XSC3_cellranger710",
                                                                  "ZMY-307-002-10XSC3_cellranger710"),'Adjacent Tumor',
                                  ifelse(sce.merge.obj$SampleType == 'Brain','Tumor',sce.merge.obj$SampleType))

sce.merge.obj$`Age` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),65,
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),54,
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),52,
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),35,
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),41,
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),53,
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),33,
                      ifelse(grepl("GSM5645891", sce.merge.obj$SampleID),75,
                      ifelse(grepl("GSM5645892", sce.merge.obj$SampleID),38,
                      ifelse(grepl("GSM5645893", sce.merge.obj$SampleID),72,
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),41,
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),55,
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),56,NA)))))))))))))

sce.merge.obj$`Time to BM (months)` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),0,
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),60,
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),30,
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),17,
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),55,
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),50,
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),48,
                      ifelse(grepl("GSM5645891", sce.merge.obj$SampleID),31,
                      ifelse(grepl("GSM5645892", sce.merge.obj$SampleID),26,
                      ifelse(grepl("GSM5645893", sce.merge.obj$SampleID),52,NA))))))))))

sce.merge.obj$`Number of BM` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'>1',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'>1',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),"1",
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'1',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'1',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),'>1',
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'>1',NA)))))))

sce.merge.obj$`Primary ER` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),'Pos',NA))))))))))

sce.merge.obj$`Primary PR` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),'Pos',NA))))))))))

sce.merge.obj$`Primary HER2` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),'Pos',NA))))))))))

sce.merge.obj$`BM ER` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),NA,
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),NA,
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),'Pos',NA))))))))))

sce.merge.obj$`BM PR` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),NA,
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),NA,
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),'Neg',NA))))))))))

sce.merge.obj$`BM HER2` = ifelse(grepl("JFP01_", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("ZY-307", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("LXZ-307", sce.merge.obj$SampleID),NA,
                      ifelse(grepl("ZMY-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("MZH-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("LHP-307", sce.merge.obj$SampleID),NA,
                      ifelse(grepl("CTT-307", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("GSM7475324", sce.merge.obj$SampleID),'Pos',
                      ifelse(grepl("GSM7475325", sce.merge.obj$SampleID),'Neg',
                      ifelse(grepl("GSM7475326", sce.merge.obj$SampleID),'Pos',NA))))))))))

meta = sce.merge.obj@meta.data[,-c(1:3)]
sce = CreateSeuratObject(counts = sce.merge.obj@assays$RNA@counts,
                         min.features = 300,
                         min.cells = dim(sce.merge.obj)[2]*0.001)
sce = AddMetaData(sce,metadata = meta)

sce[["pMT"]] = PercentageFeatureSet(object = sce, pattern = "^MT-")
sce[["pHB"]] = PercentageFeatureSet(sce,  pattern = "^HBA|^HBB")
sce[["pRP"]] = PercentageFeatureSet(sce, pattern = "^RP[SL][[:digit:]]")
sce
nFeature_lower = 300
nFeature_upper = 10000
nCount_lower = 1000
nCount_upper = 60000
pMT_lower = 0
pMT_upper = 20
pHB_lower = 0
pHB_upper = 5


sce = subset(sce, subset =  nFeature_RNA < nFeature_upper & nCount_RNA < nCount_upper & pMT < pMT_upper)
table(sce$SampleID)
table(sce$SampleType)
sce$orig.ident = sce$SampleType
sce$orig.ident = ifelse(sce$SampleID %in% c("230269A_JFP02_cellranger710",
                                                                "ZY-307-002-10XSC3_cellranger710",
                                                                "ZMY-307-002-10XSC3_cellranger710"),'Adjacent Tumor',
                                  ifelse(sce$SampleType == 'Brain','Tumor',sce$SampleType))
if(T){
  # sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 100000)
  # sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3500)
  # sce <- ScaleData(sce)
  sce = SCTransform(sce, method = "glmGamPoi", 
                    vars.to.regress = c("pMT"),seed.use = 777,
                    ncells=dim(sce)[2]*0.2,variable.features.n = 3500)
  sce = Seurat::RunPCA(sce)
  sce = FindNeighbors(sce, reduction = "pca", dims = 1:30)
  sce = FindClusters(sce,resolution = 1)
  sce = RunUMAP(sce, reduction = "pca", dims = 1:30,learning.rate = 0.5)
}
DimPlot(sce,group.by = 'SampleID',raster=FALSE)
DimPlot(sce, group.by = 'SampleType',raster=FALSE)
CKmix = c('KRT7','KRT8','KRT18','KRT19')
CellMarkers = list(
  'Epi' = c(CKmix,'EPCAM','CDH1','CD24'), # Epithelial
  'Neu' = c('CSF3R','S100A8','S100A9'),
  "T" = c("CD3D",'CD3E','CD2'),
  'AS' = c('GFAP','SLC1A2','AQP4'), # Astrocytes
  'OD'= c('MBP','MOBP','MOG'), # Oligodendrocytes
  'Maph' = c('AIF1','LYZ','CD163','PTPRC'),
  'MG' = c('TMEM119','CSF1R','ITGAM'),# Microglia
  'vSMC' = c('RGS5','TAGLN','ACTA2','ISLR','CTHRC1'), # vascular smooth muscle cells
  'End' = c('CLDN5','PECAM1','RAMP2'), # Endothelial
  'Others' = c("NR2F2", 'FLT4', "VWF")
  )
DotPlot(sce,features = (CellMarkers),cols = c("#FFFFFF", "blue"),col.min = 0) + RotatedAxis() + 
  NoLegend()
sce$celltype = NA
sce$celltype = ifelse(Idents(sce) %in% c(59),'Astrocytes',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(4,24,31,40,52,57),'Oligodendrocytes',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(0,2,3,6,8,11,17,20,22,23,32,33,34,35,37,38,39,41,42,45,47,48,50,51,53,54,55,58,62),'MTC',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(1,12,13,19,25,28,29,44,49),'Macrophage',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(9,15,18,60),'Microglia',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(7,30,46),'vSMC',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(36,56),'Endothelial',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(5,10,14,21,27,43,61),'T Cell',sce$celltype)
sce$celltype = ifelse(Idents(sce) %in% c(16,26,60,62),'Neutrophils',sce$celltype)
sce$celltype = factor(sce$celltype, levels = c("Oligodendrocytes","Astrocytes","Endothelial","vSMC","T Cell","Microglia","Macrophage","Neutrophils","MTC"))

DimPlot(sce,raster = F,label = T)
DimPlot(sce,group.by ='celltype',raster = F,label = T)

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

x=jjDotPlot(object = sce,
          gene = unlist(CellMarkers),
          xtree = F,
          ytree = F,
          id = 'celltype',
          rescale = T,
          rescale.min = -5,
          rescale.max = 5,
          point.geom = F,
          tile.geom = T)
x
cls = c("#224C5E", "#90353b", "#55752f", "#e37e00", "#6e8e84", 
        "#c10534", "#938dd2", "#cac27e", "#a0522d", "#7b92a8", 
        "#2d6d66", "#9c8847", "#bfa19c", "#ffd200", "#1f77b4", 
        "#4dbbd5", "#00a087", "#3c5488", "#f39b7f", "#8491b4", 
        "#91d1c2", "#dc0000", "#7e6148", "#e64b35", "#698ec3",
        "#fb9a99", "#fdbf6f", "#9467bd", "#FAEA71", "#6a3d9a")
cls = c('#FAEA71', '#E68649', '#78C679', '#BD82BD', '#E0181A',
        '#B1DB68', '#A6CFE3', '#1F75AC', '#D6BBDC', '#8AA5B0',
        '#FFD2A5')
cls = c('#FDD1A1', '#8AA5AF', '#FCCDE5', '#E41A1B', '#A6CEE3',
        '#1F78B4', '#B3DE69', '#FDED6F', '#BC80BD', '#A65629',
        '#FD8D3B', '#66C2A5', '#CCEBC5')
cls = c('#FEEA9A', '#88D0A6', '#FD94A5', '#6989BE', '#9467bd')
cls = c('#224C5E', '#98B597', '#EDC877', '#E1934F', '#DF7E66',
        '#B75348', '#6D2F20')
show_col(cls)
cell.cls = colorRampPalette(brewer.pal(8, "Set1"))(9)
show_col(cell.cls)


Idents(sce) = sce$celltype
DimPlot(sce,raster=FALSE,label = F) + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = cell.type.cls) + NoLegend()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')
a = DimPlot(sce,label =F)+theme(legend.position = 'bottom') +
  scale_color_manual(values = cell.type.cls)
saveRDS(sce,file = 'Combined.SeuratObj')

sce = Combined.SeuratObj
DimPlot(sce,raster=FALSE,label = F,group.by='orig.ident') + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = orig.cls)+ NoLegend()+
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')

a = DimPlot(sce,raster=FALSE,label = F,group.by='orig.ident') + 
  theme(legend.position = 'bottom') +
  scale_color_manual(values = orig.cls)


DimPlot(sce,raster=FALSE,label = F,group.by='Gender') + 
  theme(legend.position = 'bottom') +NoLegend()+
  scale_color_lancet(alpha = 0.7)+ 
  theme(
    axis.title = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')

sce$GEO = ifelse(grepl('^GSM',sce$Sample.Name),'Other Study','This Study')
DimPlot(sce,raster=FALSE,label = F,group.by='GEO') + 
  theme(legend.position = 'bottom') +
  scale_color_simpsons(alpha = 0.7)+ NoLegend()+
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.line = element_blank()) + ggtitle('')

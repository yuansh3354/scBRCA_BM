
rm(list = ls())
gc(reset = T)
# #packages
# pkgs <- c('scRNAtoolVis','ClusterGVis','scrabble','monocle3','tidyverse')
# lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
source('script/Functions.R')
options(stringsAsFactors = F)
Combined.SeuratObj = readRDS('Combined.SeuratObj')
sce.std = readRDS('MTC.STAD')
Combined.SeuratObj$Sample.Name = gsub('(_cellranger710)|(-10XSC3_cellranger710)','',Combined.SeuratObj$SampleID)
Combined.SeuratObj$cell_Subtype = ifelse(Combined.SeuratObj$celltype == 'MTC',
                                         as.character(sce.std$cell_Subtype), as.character(Combined.SeuratObj$celltype))
meta = Combined.SeuratObj@meta.data
plot.path='plts/'
file.path = 'result/inferCNV/'
#dir.create(file.path)
# reticulate::use_condaenv("scPytorch",conda = '/home/yuansh/anaconda3/bin/conda')
# reticulate::py_config()
# reticulate::conda_binary(conda='/home/yuansh/anaconda3/bin/conda')

library(infercnv)
####################################################
if(F){
  sce = Combined.SeuratObj
  sce = subset(sce, subset = Sample.Name %in% c("GSM4259357","GSM4555888","GSM4555889","GSM4555891","GSM5645891",
                                                "GSM5645892","GSM5645893","GSM6123277","GSM7475324","GSM7475325",
                                                "GSM7475326"))
  expFile=paste0(file.path,'Other_expFile.txt')
  groupFiles=paste0(file.path,'Other_groupFiles.txt')
  geneFile=paste0(file.path,'Other_geneFile.txt')
  unique(sce$cell_Subtype)
  if(sum(file.exists(c(expFile,groupFiles,geneFile))) != 3){
    dfcount = as.data.frame(sce@assays$SCT@counts)
    groupinfo= data.frame(cellId = colnames(dfcount))
    groupinfo$cellType = sce$cell_Subtype
    geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
    geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
    geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
    dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
    dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
    write.table(dfcount ,file = expFile,sep = '\t',quote = F)
    write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
    write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
    other.infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                              annotations_file=groupFiles,
                                              delim="\t",
                                              gene_order_file= geneFile,
                                              ref_group_names=c("Oligodendrocytes","vSMC",
                                                                "Microglia","Macrophage","T Cell",
                                                                "Endothelial","Neutrophils")) # 如果有正常细胞的话，把正常细胞的分组填进去
    # Run infer CNV
    saveRDS(other.infercnv_obj,file = 'result/inferCNV/Other-infercnv_obj')
  }

  
  sce = Combined.SeuratObj
  sce = subset(sce, subset = Sample.Name %in% c("230269A_JFP01","230269A_JFP02","CTT-307-CTC","LHP-307-001",
                                                "LXZ-307-001","LXZ-307-CTC","MZH-307-001","ZMY-307-001","ZY-307-001",
                                                "ZY-307-002","LHP-307-CTC","ZMY-307-002"))
  expFile=paste0(file.path,'My_expFile.txt')
  groupFiles=paste0(file.path,'My_groupFiles.txt')
  geneFile=paste0(file.path,'My_geneFile.txt')
  
  if(sum(file.exists(c(expFile,groupFiles,geneFile))) != 3){
    dfcount = as.data.frame(sce@assays$SCT@counts)
    groupinfo= data.frame(cellId = colnames(dfcount))
    groupinfo$cellType = sce$cell_Subtype
    geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
    geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
    geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
    dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
    dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
    write.table(dfcount ,file = expFile,sep = '\t',quote = F)
    write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
    write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
    my.infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                           annotations_file=groupFiles,
                                           delim="\t",
                                           gene_order_file= geneFile,
                                           ref_group_names=c("Oligodendrocytes","vSMC","Astrocytes",
                                                             "Microglia","Macrophage","T Cell",
                                                             "Endothelial","Neutrophils")) 
    # Run infer CNV
    saveRDS(my.infercnv_obj,file = 'result/inferCNV/My-infercnv_obj')
  }
  
  sce = Combined.SeuratObj
  expFile=paste0(file.path,'all_expFile.txt')
  groupFiles=paste0(file.path,'all_groupFiles.txt')
  geneFile=paste0(file.path,'all_geneFile.txt')
  
  if(sum(file.exists(c(expFile,groupFiles,geneFile))) != 3){
    dfcount = as.data.frame(sce@assays$SCT@counts)
    groupinfo= data.frame(cellId = colnames(dfcount))
    groupinfo$cellType = sce$cell_Subtype
    geneInfor=annoGene(rownames(dfcount),"SYMBOL",'human')
    geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
    geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
    dfcount =dfcount [rownames(dfcount ) %in% geneInfor[,1],]
    dfcount =dfcount [match( geneInfor[,1], rownames(dfcount) ),] 
    write.table(dfcount ,file = expFile,sep = '\t',quote = F)
    write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
    write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
    my.infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                           annotations_file=groupFiles,
                                           delim="\t",
                                           gene_order_file= geneFile,
                                           ref_group_names=c("Oligodendrocytes","vSMC","Astrocytes",
                                                             "Microglia","Macrophage","T Cell",
                                                             "Endothelial","Neutrophils")) 
    # Run infer CNV
    saveRDS(my.infercnv_obj,file = 'result/inferCNV/all-infercnv_obj')
  }
}



##### check 
#  Import inferCNV dendrogram
if(F){
  infercnv.obs = read.table('result/inferCNV/My/infercnv.observations.txt',header = T)
  infercnv.ref = read.table('result/inferCNV/My/infercnv.references.txt',header = T)
  df = Combined.SeuratObj@meta.data
  colnames(infercnv.obs) = gsub('\\.','-',colnames(infercnv.obs))
  colnames(infercnv.ref) = gsub('\\.','-',colnames(infercnv.ref))
  gene.anno = read.table('result/inferCNV/My_geneFile.txt')
  infercnv.combine = cbind(infercnv.ref,infercnv.obs) %>% as.data.frame()
  ids = intersect(gene.anno$V1,rownames(infercnv.combine))
  gene.anno = gene.anno[which(gene.anno$V1 %in% ids),]
  chr.list = split(gene.anno$V1,gene.anno$V2)
  chr.score = lapply(seq_along(unique(names(chr.list))), function(x){
    temp = infercnv.combine[chr.list[[x]],] %>% apply(.,2,median)
    return(temp)
  })
  names(chr.score) = names(chr.list)
  df = as.data.frame(chr.score)
  use.cells = intersect(rownames(df),rownames(meta))
  chr.score.df = cbind(meta[use.cells,],df)
  saveRDS(chr.score.df,'my.infercnv.chr.score.df')
  saveRDS(infercnv.obs,'my.infercnv.obs')
  saveRDS(infercnv.ref,'my.infercnv.ref')
}
if(F){
  infercnv.obs = read.table('result/inferCNV/Other/infercnv.observations.txt',header = T)
  infercnv.ref = read.table('result/inferCNV/Other/infercnv.references.txt',header = T)
  df = Combined.SeuratObj@meta.data
  colnames(infercnv.obs) = gsub('\\.','-',colnames(infercnv.obs))
  colnames(infercnv.ref) = gsub('\\.','-',colnames(infercnv.ref))
  gene.anno = read.table('result/inferCNV/Other_geneFile.txt')
  infercnv.combine = cbind(infercnv.ref,infercnv.obs) %>% as.data.frame()
  ids = intersect(gene.anno$V1,rownames(infercnv.combine))
  gene.anno = gene.anno[which(gene.anno$V1 %in% ids),]
  chr.list = split(gene.anno$V1,gene.anno$V2)
  chr.score = lapply(seq_along(unique(names(chr.list))), function(x){
    temp = infercnv.combine[chr.list[[x]],] %>% apply(.,2,median)
    return(temp)
  })
  names(chr.score) = names(chr.list)
  df = as.data.frame(chr.score)
  use.cells = intersect(rownames(df),rownames(meta))
  chr.score.df = cbind(meta[use.cells,],df)
  saveRDS(chr.score.df,'other.infercnv.chr.score.df')
  saveRDS(infercnv.obs,'other.infercnv.obs')
  saveRDS(infercnv.ref,'other.infercnv.ref')
}

other = readRDS('other.infercnv.chr.score.df')
my = readRDS('my.infercnv.chr.score.df')
chr.score.df = rbind(other,my)
chr.score.df$cell_Subtype  = ifelse(chr.score.df$celltype !='MTC','Normal',chr.score.df$cell_Subtype )
# ids = mean(apply(chr.score.df[which(chr.score.df$celltype!='MTC'),27:48],1,mean))
# chr.score.df = chr.score.df[which(chr.score.df$celltype == 'MTC'),]
df.plot = chr.score.df[,26:48]
df.plot = melt(df.plot)
df.plot$value = df.plot$value - 1
# df.plot$cell_Subtype = factor(df.plot$cell_Subtype,levels = cell.sub.level)
df.plot$variable = factor(df.plot$variable,levels = paste0('chr',1:22))
# df.plot$value = ifelse(df.plot$value > -0.025 & df.plot$value<0.025,0,df.plot$value)
df.plot = na.omit(df.plot)
library(ggridges)
x1 = ggplot(df.plot, aes(x = value, y = cell_Subtype, fill =cell_Subtype)) +
  stat_density_ridges(
    jittered_points=F, scale = 1, rel_min_height = 0.01,
    geom = "density_ridges_gradient", calc_ecdf = TRUE
  ) +
  theme_cleveland()+
  scale_fill_simpsons(alpha = 0.8) +
  labs(x = "Feature", y = "Cell Type") +
  facet_wrap(~variable,scales = 'free') +
  theme(panel.background = element_rect(fill = "white"))
x1

seurat.info = 'infercnv.chr.pvalue.log.txt'
if(file.exists(seurat.info)){file.remove(seurat.info)}
sink(seurat.info, append = TRUE)
for(chr in c('chr8','chr13','chr14','chr16','chr21')){
  for(cell in MTC.levels){
    print(paste0(chr,": ",cell))
    print(t.test(chr.score.df[which(chr.score.df$cell_Subtype =='Normal'),chr],chr.score.df[which(chr.score.df$cell_Subtype ==cell),chr]))
    cat("\n")
  }
}
sink()

df.plot = chr.score.df
df.plot = df.plot[which(df.plot$celltype == 'MTC'),]
df.plot = df.plot[,c(1,26:48)]
df.plot = melt(df.plot,id.vars = c("orig.ident","cell_Subtype"))
df.plot$value = df.plot$value - 1
df.plot$variable = factor(df.plot$variable,levels = paste0('chr',1:22))
df.plot = na.omit(df.plot)
df.plot$value = ifelse(df.plot$value > 0.05,'Gain',ifelse(df.plot$value< -0.05, 'Loss','None'))
df.plot = df.plot[which(df.plot$value !='None'),]
x4 = ggplot(df.plot, aes(x = value,fill =cell_Subtype)) +
  geom_bar(position = 'dodge')+
  theme_cleveland()+
  scale_fill_manual(values = MTC.cls)+
  labs(x = NULL, y = 'Number of Cells') +
  facet_wrap(~variable) +
  theme(panel.background = element_rect(fill = "white"),legend.position = 'None',
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_line(size = 0.5),)
x4





df.plot = chr.score.df[,27:48]
df.plot= df.plot -1
df.plot[df.plot>-0.025 & df.plot<0.025] =0
# df.plot[df.plot>0] = 1
# df.plot[df.plot<0] = -1
df.plot = df.plot[,c(paste0('chr',1:22))]
df.plot = df.plot[,c(5,9,7,17,8,10,13,14,16,20,21)]
df = cor(df.plot)
x2 = pheatmap::pheatmap(df,cluster_rows = F,cluster_cols = F,
                   breaks =seq(-1, 1, length.out = 101),
                   display_numbers = TRUE,number_color = "black"
                   )
df.plot = chr.score.df[,27:48]
df.plot= df.plot -1
df.plot[df.plot>-0.025 & df.plot<0.025] =0
# df.plot[df.plot>0] = 1
# df.plot[df.plot<0] = -1
df.plot = df.plot[,c(paste0('chr',1:22))]
df.plot = df.plot[,c(5,9,7,17,8,10,13,14,16,20)]
df = cor(df.plot)
x3 = pheatmap::pheatmap(df,cluster_rows = F,cluster_cols = F,
                        breaks =seq(-1, 1, length.out = 101),
                        display_numbers = TRUE,number_color = "black"
)

topptx(x1,'plts/inferCNV.pptx',height = 25,width = 40)
topptx(x2,'plts/inferCNV.cormut.pptx')
topptx(x3,'plts/inferCNV.cormut.cor.pptx',height = 15,width = 15)
topptx(x4,'plts/inferCNV-bar.pptx',height = 25,width = 40)





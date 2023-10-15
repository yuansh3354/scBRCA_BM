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
####################################################
sce.std = readRDS('MTC.STAD')
sce = sce.std
use.cells = sce@meta.data[which(sce$Sample.Name %in% c("230269A_JFP01","230269A_JFP02","CTT-307-CTC","GSM4259357","GSM4555888",
                                                       "GSM4555889","GSM4555891","GSM6123277","LHP-307-001","LXZ-307-001","LXZ-307-CTC",
                                                       "MZH-307-001","ZMY-307-001","ZY-307-001","ZY-307-002","LHP-307-CTC","ZMY-307-002")),] %>% rownames()
sce = subset(sce,cells=use.cells)
ids = names(which(table(sce$SampleID) > 50))
sce = sce[,which(sce$SampleID %in% ids)]
sce = CreateSeuratObject(counts = sce@assays$SCT@counts)
sce = AddMetaData(sce, metadata = sce.std@meta.data)
if(T){
  sce = NormalizeData(object = sce,normalization.method =  "LogNormalize",  scale.factor = 1e4)
  sce = FindVariableFeatures(object = sce,selection.method = "vst", nfeatures = 2000)
  sce = ScaleData(object = sce)
  sce = Seurat::RunPCA(object = sce, do.print = FALSE)
  sce <- RunUMAP(sce, reduction = "pca", dims = 1:5, n.components = 3)
  sce <- FindNeighbors(sce, reduction = "pca", dims = 1:5)
  sce <- FindClusters(sce)
  DimPlot(sce,group.by = 'orig.ident')
}
if(F){
  expression_matrix = sce@assays$RNA@data
  cell_metadata = data.frame(sce@meta.data)
  gene_annotation = data.frame(expression_matrix[,1])
  gene_annotation[,1] = row.names(gene_annotation)
  colnames(gene_annotation)=c("gene_short_name")
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  cds <- preprocess_cds(cds, num_dim = num_dim,norm_method = c("none"))
  cds <- align_cds(cds, alignment_group = "orig.ident")
  
  cds <- reduce_dimension(cds,cores=16)
  
  cds <- cluster_cells(cds,resolution = 0.0000001)
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(sce, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  cds <- learn_graph(cds)
  
  myselect <- function(cds,select.classify,my_select){
    cell_ids <- which(colData(cds)[,select.classify] == my_select)
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    }
  
  cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'orig.ident',my_select = "Blood"))

  plot_cells(cds, 
             color_cells_by = "pseudotime",cell_size = 1,
             label_branch_points=FALSE,show_trajectory_graph=F)+theme_classic()+
    theme(

      axis.ticks = element_blank(),
      axis.line = element_blank()) + NoLegend()
  x = plot_cells(cds, 
             color_cells_by = "pseudotime",cell_size = 0,
             label_branch_points=FALSE,show_trajectory_graph=F)+theme_classic()+
    theme(
      axis.ticks = element_blank(),
      axis.line = element_blank())
    plot_cells(cds,
               color_cells_by = "orig.ident",
               label_cell_groups=FALSE,cell_size = 1,
               label_leaves=FALSE,show_trajectory_graph=F,
               label_branch_points=FALSE)+theme_classic() +scale_color_manual(values = orig.cls)+
    theme(
      axis.title = element_blank(),  
      axis.text = element_blank(), 
      axis.ticks = element_blank(),
      axis.line = element_blank()) + NoLegend()
  
  Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=16)
  saveRDS(cds,file = 'cds.RDS')
  saveRDS(Track_genes,file='Track_genes.RDS')
}
cds = readRDS('cds.RDS')
Track_genes = readRDS('Track_genes.RDS')
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
df = colData(cds)
pseudo.df = as.data.frame(pseudotime(cds)) %>% set_colnames('pseudotime')
pseudo.df$pseudotime.raw = pseudo.df$pseudotime
pseudo.df$pseudotime = round(pseudo.df$pseudotime)
# sce = sce.std
sce = AddMetaData(sce, metadata = pseudo.df)

df = sce@meta.data
df$cell_Subtype = factor(df$cell_Subtype, 
                         levels = c('MTC blank','Neuro Related',
                                    'HR+','AR+HR+','HER2+' ))
#which(is.na(df$pseudotime.raw) | is.infinite(df$pseudotime.raw))
df = df[,c('pseudotime.raw','orig.ident','cell_Subtype')]
df = na.omit(df)
sce$orig.ident = factor(sce$orig.ident,levels = c('Blood','CSF','Tumor','Adjacent Tumor' ))
x1 = RidgePlot(sce, features = 'pseudotime',group.by = 'orig.ident')+
  scale_fill_manual(values = orig.cls)

df = sce@meta.data
df.plot = df %>%
  group_by(pseudotime,cell_Subtype) %>%
  summarise(Count = n()) %>% 
  mutate(percentage = Count/sum(Count)) 
df.plot[which(df.plot$pseudotime == 1),'percentage'] %>% sum
df.plot$cell_Subtype = factor(df.plot$cell_Subtype, 
                              levels = c('MTC blank','Neuro Related',
                                         'HR+','AR+HR+','HER2+' ))
x2 = ggplot(df.plot, 
       aes(pseudotime,percentage))+ 
  facet_wrap(~cell_Subtype,nrow = 1)+
  theme_bw()+
  geom_point(size = 0.5)+
  geom_smooth(method = "loess", span =2)+
  stat_cor(method = "spearman",  label.x = 0,
           label.y = 0.9)

sce$Sample.Name = factor(sce$Sample.Name, levels = sample.levels)
df.plot = df %>%
  group_by(Sample.Name,cell_Subtype) %>%
  summarise(Count = n())%>% 
  mutate(percentage = Count/sum(Count)) 

x4 = ggplot(df.plot, aes(Sample.Name, percentage,fill=cell_Subtype)) + facet_wrap(~cell_Subtype,ncol = 1)+
  geom_bar(stat = 'identity') +
  theme(axis.text.x =element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), 
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")) +   scale_fill_manual(values =MTC.cls)

library(ClusterGVis)
genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.35))
set.seed(9463)
random_numbers <- sample(1:36022, 15000, replace = FALSE)
mat <- pre_pseudotime_matrix(cds_obj = cds[,random_numbers],
                             gene_list = genes)
ck <- clusterData(exp = mat,
                  cluster.method = "kmeans",
                  cluster.num = 3)
x = visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = sample(rownames(mat),30,replace = F))
  
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
head(enrich,3)

enrich.KEGG <- enrichCluster(object = ck,
                             OrgDb = org.Hs.eg.db,
                             type = "KEGG",
                             organism = "hsa",
                             pvalueCutoff = 0.9,
                             topn = 5,
                             seed = 5201314)

palette = c("Grays","Blues2","Reds2","Purples2","Greens2")
lapply(seq_along(unique(enrich$group)), function(x){
  # go plot
  tmp <- enrich |> dplyr::filter(group == unique(enrich$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description,levels = tmp$Description)
  
  # plot
  p <-
    ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue),y = Description,fill = -log10(pvalue)),
             width = 0.75) +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~.,name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  # plot kegg
  tmp.kg <- enrich.KEGG |> dplyr::filter(group == unique(enrich.KEGG$group)[x]) |>
    dplyr::arrange(desc(pvalue))
  
  tmp.kg$Description <- factor(tmp.kg$Description,levels = tmp.kg$Description)
  
  # plot
  pk <-
    ggplot(tmp.kg) +
    geom_segment(aes(x = 0,xend = -log10(pvalue),y = Description,yend = Description),
                 lty = "dashed",linewidth = 0.75) +
    geom_point(aes(x = -log10(pvalue),y = Description,color = -log10(pvalue)),size = 5) +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    colorspace::scale_color_binned_sequential(palette = palette[x]) +
    ylab("") + xlab("-log10(pvalue)")
  
  # combine
  cb <- cowplot::plot_grid(plotlist = list(p,pk))
  
  return(cb)
}) -> gglist

x1 = gglist[[1]]
x2 = gglist[[2]]
x3 = gglist[[3]]
plt = ggarrange(x1,x2,x3, ncol = 1, labels = c("a", "b",'c'))
plt

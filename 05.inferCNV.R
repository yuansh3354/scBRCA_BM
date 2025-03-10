
####################################################
if(T){
  file.path='results'
  sce = readRDS('data/Combined.SeuratObj.RDS')
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
    saveRDS(my.infercnv_obj,file = 'results/infercnv_obj')
  }
}
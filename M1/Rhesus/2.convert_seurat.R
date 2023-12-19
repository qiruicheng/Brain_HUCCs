library(dplyr)
library(Seurat)
# library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(purrr)
options(Seurat.object.assay.version = 'v4')

sample_id <- list.files("results/")
dir <- "/picb/neurosys/qiruicheng/cellTypeEvolutionsWithBulkandSCdata/data/m1/rhesus/results/"
seu_m1 = lapply(sample_id, function(x){
  print(x)
  dir1 = paste0(dir,"/",x,"/filtered_feature_bc_matrix/")
  files = dir(dir1)
  
  #### counts
  # counts_pth = paste0(dir1,grep('mtx',files,value = T))
  # counts = Matrix::readMM(counts_pth)
  counts_pth = paste0(dir,"/",x,"/filtered_feature_bc_matrix.h5")
  counts = Read10X_h5(counts_pth)
  ### read genes info
  # gene_pth = paste0(dir1,grep('features',files,value = T))
  # gene = read.table(gene_pth,header=F)
  ### read cell info
  meta_pth = paste0(dir1,grep('barcodes',files,value = T))
  meta = read.delim(meta_pth,header = F)
  row.names(meta) = meta[,1]

  #### annot counts data
  # row.names(counts) = gene[,2]
  # colnames(counts) = meta[,1]
  
  rownames(meta) <- meta[,1]
  meta[,1] <- x
  colnames(meta) <- "orig.ident"
  #### Create SeuObj
  seu = CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 3, min.features = 200)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu = subset(seu, subset = nFeature_RNA > 200 & percent.mt < 5)
  Idents(seu) = seu$orig.ident
  seu@meta.data$cellID = rownames(seu@meta.data)
  return(seu)
})

names(seu_m1) <- sample_id

### functions
my_ref_int_for_largedata <- function(seu){
  

  seu_sep = lapply(seu, function(x){
    fin = NormalizeData(x) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 3000) %>% ScaleData()
    return(fin)
  })
  features <- SelectIntegrationFeatures(object.list = seu_sep,nfeatures = 3000)
  seu_sep <- lapply(X = seu_sep, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })

  seu_anchors = FindIntegrationAnchors(seu_sep, reduction = "rpca", reference=1, dims = 1:50, k.anchor = 30, anchor.features = 3000)
  seu_integrate = IntegrateData(seu_anchors, dims = 1:50)
  seu_integrate = seu_integrate %>% ScaleData() %>% RunPCA() %>% RunUMAP(reduction = 'pca', dims = 1:50) %>% FindNeighbors(reduction = "pca", dims = 1:50)

  return(seu_integrate)
}
# 20459 238699
# sc_m1 <- merge(seu_m1[[1]],y=c(seu_m1[[2]],seu_m1[[3]],seu_m1[[4]],seu_m1[[5]],seu_m1[[6]],seu_m1[[7]],seu_m1[[8]],
#                                seu_m1[[9]],seu_m1[[10]],seu_m1[[11]],seu_m1[[12]],seu_m1[[13]],seu_m1[[14]],seu_m1[[15]],seu_m1[[16]],
 #                               seu_m1[[17]],seu_m1[[18]],seu_m1[[19]],seu_m1[[20]],seu_m1[[21]],seu_m1[[22]],seu_m1[[23]],seu_m1[[24]],seu_m1[[25]]))
# rm(seu_m1)
# gc()
sc_m1 <- my_ref_int_for_largedata(seu_m1)
saveRDS(sc_m1,"rhesus.m1.rds")

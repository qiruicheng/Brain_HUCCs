library(dplyr)
library(Seurat)
# library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(purrr)
library(rlist)
options(Seurat.object.assay.version = 'v4')

sample_id <- list.files("results/")
dir <- "/picb/neurosys/qiruicheng/cellTypeEvolutionsWithBulkandSCdata/data/m1/chimpanzee/aggr_chimpanzee/outs/count"
counts_pth = paste0(dir,"/filtered_feature_bc_matrix.h5")
counts = Read10X_h5(counts_pth)

dir <- "/picb/neurosys/qiruicheng/cellTypeEvolutionsWithBulkandSCdata/data/m1/chimpanzee/results/"
meta = lapply(sample_id, function(x){
  print(x)
  dir1 = paste0(dir,"/",x,"/filtered_feature_bc_matrix/")
  files = dir(dir1)
  meta_pth = paste0(dir1,grep('barcodes',files,value = T))
  meta = read.delim(meta_pth,header = F)
  rownames(meta) <- meta[,1]
  meta[,1] <- x
  colnames(meta) <- "orig.ident"
  return(meta)
})
meta <- meta %>% list.rbind()%>%data.frame()
colnames(counts) <- rownames(meta)
#### Create SeuObj
seu_m1 = CreateSeuratObject(counts = counts, meta.data = meta, min.cells = 3, min.features = 200)
# seu_m1[["percent.mt"]] <- PercentageFeatureSet(seu_m1, pattern = "^MT-")
# seu_m1 = subset(seu_m1, subset = nFeature_RNA > 200 & percent.mt < 5)
Idents(seu_m1) = seu_m1$orig.ident
seu_m1@meta.data$cellID = rownames(seu_m1@meta.data)

my_ref_int_for_largedata <- function(seu){

  seu_sep = SplitObject(seu, split.by = "orig.ident")
  seu_sep = lapply(seu_sep, function(x){
    fin = x
    fin[["percent.mt"]] <- PercentageFeatureSet(fin, pattern = "^MT-")
    fin = subset(fin, subset = nFeature_RNA > 200 & percent.mt < 5)
    fin = NormalizeData(fin) %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 3000) %>% ScaleData()
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

seu_m1 <- my_ref_int_for_largedata(seu_m1)
saveRDS(seu_m1,"chimpanzee.m1.rds")
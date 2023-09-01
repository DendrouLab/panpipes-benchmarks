start = Sys.time()
print("started at")
print(start)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(tictoc)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")


print(Sys.time()-start)
obj <- readRDS("multimodal_pbmcv4.rds")
print("starting with RNA integrations")
#Warning message:
#  Key ‘originalexp_’ taken, using ‘rna_’ instead 
DefaultAssay(obj) <- "RNA"

#split by covariate to integrate on
tic()
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_id)
obj
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
print("after PCA time elapsed")
toc()
tic()
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
print("after UMAP time elapsed")
toc()
p1 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("sample_id", "celltype.l2"))

print("Harmony starts..")
tic()

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
print("Harmony ends")
toc()

print("scvi starts...")
tic()

obj <- IntegrateLayers(
  object = obj, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "~/miniconda3/envs/panpipes_env", verbose = FALSE
)
print("scvi ends")
toc()

tic()
obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")
print("scvi neighbours ends")
toc()


p2 <- DimPlot(
  obj,
  reduction = "umap.scvi",
  group.by = c("sample_id", "celltype.l2"),
  combine = FALSE
)

print("harmony neighbours starts")
tic()


obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p3 <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("sample_id", "celltype.l2"),
  combine = FALSE
)

print("harmony neighbours ends")
toc()
# ggsave(p1, file = "no_integration_umap.pdf")
# ggsave(p2, file = "scvi_integration_umap.pdf")
# ggsave(p3, file = "harmony_integration_umap.pdf")

start = Sys.time()
DefaultAssay(obj) <- "PROT"
print("protein integration starts")
print(start)
#correct proteins

tic()

obj[["PROT"]] <- split(obj[["PROT"]], f = obj$sample_id)
obj
obj <- NormalizeData(obj, normalization.method="CLR")
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
print("after PCA time elapsed")
toc()
tic()
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
print("after UMAP time elapsed")
toc()
p1 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = c("sample_id", "celltype.l2"))

print("Harmony starts..")
tic()

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  assay="PROT",
  orig.reduction = "pca", key="alternative_",new.reduction = "alternative",
  verbose = FALSE
)
print("Harmony ends")
toc()

print("harmony neighbours starts")
tic()


obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p3 <- DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("sample_id", "celltype.l2"),
  combine = FALSE
)

print("harmony neighbours ends")
toc()


print("running WNN")
tic()

bm <- FindMultiModalNeighbors(
  obj, reduction.list = list("harmony", "pca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight"
)

toc()

print("running multimodal umap")
tic()
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
print("finished")
toc()

p4 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
# ggsave(p4, file = "wnn_integration_umap.pdf")
saveRDS(bm, file = "integrated_object.rds")
print("Finished all")
print(Sys.time())
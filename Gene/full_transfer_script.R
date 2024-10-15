
library(data.table)
library(ggplot2)
library(monocle3)
library(dplyr)
library(Seurat)

options(future.globals.maxSize = 8000 * 1024^2)

full_mtx <- readRDS("~/aneuploidy/data/full_cds.RDS")
full_obs <- as.data.frame(colData(full_mtx))


seurat_cao_dis <- CreateSeuratObject(counts = assay(full_mtx), meta.data = full_obs)
seurat_cao_dis[["RNA"]] <- split(seurat_cao_dis[["RNA"]], f = seurat_cao_dis$origin)

seurat_cao_dis <- NormalizeData(seurat_cao_dis)
seurat_cao_dis <- FindVariableFeatures(seurat_cao_dis)
seurat_cao_dis <- ScaleData(seurat_cao_dis)
seurat_cao_dis <- RunPCA(seurat_cao_dis)

seurat_cao_dis <- FindNeighbors(seurat_cao_dis, dims = 1:50, reduction = "pca")
seurat_cao_dis <- FindClusters(seurat_cao_dis, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_cao_dis <- RunUMAP(seurat_cao_dis, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")

seurat_cao_dis <- IntegrateLayers(
  object = seurat_cao_dis, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

seurat_cao_dis <- FindNeighbors(seurat_cao_dis, reduction = "integrated.cca", dims = 1:50)
seurat_cao_dis <- FindClusters(seurat_cao_dis, resolution = 2, cluster.name = "cca_clusters")

integration_ref <- subset(seurat_cao_dis, origin %in% c("Cao"))
integration_query <- subset(seurat_cao_dis, origin %in% c("Disteche", "Ian"))

int_test <- FindTransferAnchors(reference = integration_ref, query = integration_query, dims = 1:50,
    reference.reduction = "integrated.cca")
predictions_int <- TransferData(anchorset = int_test, refdata = integration_ref$Main_cluster_name, dims = 1:50)
integration_query_int <- AddMetaData(integration_query, metadata = predictions_int)

use_df <- as.data.frame(seurat_cao_dis@meta.data)
# Select the columns "cell_unique" and "predicted.id" to create a new dataframe
new_df <- integration_query_int@meta.data %>% 
  select(cell_unique, predicted.id)

# Left join obs_df_main with the new dataframe on the column "cell_unique"
use_df <- left_join(use_df, new_df, by = "cell_unique")

seurat_cao_dis@meta.data$predicted_id <- use_df$predicted.id

final_obs <- as.data.frame(seurat_cao_dis@meta.data)
final_var <- as.data.frame(dis_cao_cds@rowRanges@elementMetadata)
rownames(final_var) <- final_var$ensembl
exp_norm <- cbind2(seurat_cao_dis@assays[["RNA"]]@layers[["data.Cao"]], seurat_cao_dis@assays[["RNA"]]@layers[["data.Disteche"]])
exp_norm <- cbind2(exp_norm, seurat_cao_dis@assays[["RNA"]]@layers[["data.Ian"]])

# creates CDS object
dis_cao_cds_final <- new_cell_data_set(exp_norm,
                         cell_metadata = final_obs,
                         gene_metadata = full_var)

saveRDS(dis_cao_cds_final, "~/aneuploidy/data/final_cds_full.RDS")

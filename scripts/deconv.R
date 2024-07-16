library(spacedeconv)

seurat_downsize <- subset(in_data, downsample = 500)  # picks only 2000 cells per cluster
sce_downsize <- as.SingleCellExperiment(seurat_downsize)

sce_downsize <- readRDS("1.rds")
spatial_obj <- readRDS("2.rds")

#rctd GOOD
deconv_rctd <- build_and_deconvolute(
  sce_downsize,
  spatial_obj,
  method = "rctd",
  cell_type_col = "celltype_subset",
  batch_id_col = NULL,
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE,
  n_cores = 8
)

#spotlight GOOD
deconv_spotlight <-build_and_deconvolute(
  sce_downsize,
  spatial_obj,
  method = "spotlight",
  cell_type_col = "celltype_subset",
  batch_id_col = NULL,
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE
)


#card GOOD
devtools::install_github('YingMa0107/CARD')
library(CARD)

deconv_card <-build_and_deconvolute(
  sce_downsize,
  spatial_obj,
  method = "card",
  cell_type_col = "celltype_subset",
  batch_id_col = "orig.ident",
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE
)




#spatialdwls
giotto_instructions$python_path <- "/usr/bin/python"


deconv_spatialdwls <-build_and_deconvolute(
  sce_downsize,
  spatial_obj,
  method = "spatialdwls",
  cell_type_col = "celltype_subset",
  batch_id_col = "orig.ident",
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE
)

#cell2location
deconv_cell2location <-build_and_deconvolute(
  sce_downsize,
  spatial_obj,
  method = "cell2location",
  cell_type_col = "celltype_subset",
  batch_id_col = "orig.ident",
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE
)




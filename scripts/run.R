# Rscript --vanilla ./scripts/run.R --scdata "data/scRNA_wu" --scmeta "data/scRNA_wu/metadata.csv" --outdir ./data/results/test/run_1 --seed 0 --downsize 500 --grain_lvl "celltype_major" --gene_column 1 --synth_dataset "artificial_regional_rare_celltype_diverse" --select_celltype "T-cells"
# source("/localdisk/home/s2600569/benchdeconv/benchdeconv/scripts/run.R")
start_time <- Sys.time()

library(argparser)

print("------STARTING benchdeconv------")
print(getwd())

#argparser args
print("---Parsing input args---")

# Create a parser
input_args <- arg_parser("benchdeconv: a benchmarking tool for spatial deconvolution methods.")

# Add command line arguments
input_args <- add_argument(input_args, "--scdata", help="Input single cell data directory. Requires count matrix barcodes, genes, counts in a sparse matrix, and metadata of cell annotations",
                           type="character",
                           default = "./data/scRNA_wu")

input_args <- add_argument(input_args, "--scmeta", help="Metadata for cell annotations", 
                           type="character",
                           default = "./data/scRNA_wu/metadata.csv")

input_args <- add_argument(input_args, "--outdir", help="Output directory for results",
                           type="character",
                           default = "./data/results")

input_args <- add_argument(input_args, "--seed", help="Setting seed for random sampling",
                           type = "numeric",
                           default = 1)

input_args <- add_argument(input_args, "--downsize", help = "No. of cells to downsize scRNA reference to. (0 for full dataset).",
                           type = "numeric",
                           default = 0)

input_args <- add_argument(input_args, "--grain_lvl", help="Column used for celltype annotations",
                           type = "character",
                           default = "celltype_major")

input_args <- add_argument(input_args, "--gene_column", help="Gene column in input scRNA data",
                           type = "numeric",
                           default = 1)

input_args <- add_argument(input_args, "--synth_dataset", help="Dataset setting for synthspots",
                           type = "character",
                           default = "artificial_regional_rare_celltype_diverse")

input_args <- add_argument(input_args, "--select_celltype", help="select_celltype option in synthspots",
                           type = "character",
                           default = "T-cells")

input_args <- add_argument(input_args, "--n_cells_max", help="max no. of cells per spot in synthspot",
                           type = "numeric",
                           default = 40)

input_args <- add_argument(input_args, "--min_cell_id_test", help="cell density frequency for control spots used to test minimum cell density for identification (0-1)",
                           type = "numeric",
                           default = 0)

input_args <- add_argument(input_args, "--select_celltype_min_id", help="Selected celltype for the minimum cell density for identification test, if set to 1",
                           type = "character",
                           default = "T-cells")

input_args <- add_argument(input_args, "--subtype", help="subtype to split on with input scRNA data",
                           type = "character",
                           default = "HER2+")

input_args <- add_argument(input_args, "--coords", help="coordinate files that have been drawn via Loupe Browser v8",
                           type = "character",
                           default = "./data/spot_coords/out1.csv")

input_args <- add_argument(input_args, "--coords_total", help="reference total coords from Loupe Browser v8",
                           type = "character",
                           default = "./data/spot_coords/Spatial-Projection.csv")


argv <- parse_args(input_args)
print("Settings:")
# test settings
argv <- list()
argv$scdata <- "data/scRNA_wu"
argv$scmeta <- "data/scRNA_wu/metadata.csv"
argv$outdir <- "./data/results/test/run_1"
argv$seed <- 0
argv$downsize <- 500
argv$grain_lvl <- "celltype_major"
argv$gene_column <- 1
argv$synth_dataset <- "artificial_diverse_overlap"
argv$select_celltype <- "T-cells"
argv$n_cells_max <- 40
argv$min_cell_id_test <- 0.1
argv$select_celltype_min_id <- "T-cells"
argv$subtype <- "HER2+"
argv$coords <- "./data/spot_coords/out1_mintest.csv"
argv$coords_total <- "./data/spot_coords/Spatial-Projection.csv"
print(argv)
print("Parsing done.")




# Set seed for reproducibility
print(paste0("---Setting seed = ", argv$seed, "---"))
set.seed(argv$seed) 
print("Seed setting done.")

print("---Sourcing benchdeconv.R---")
# install_local("../synthspot_devbuild_0.1", force = FALSE)
source("./scripts/benchdeconv.R")
print("Source done.")


#create dirs
print("---Creating dirs---")
#if (!dir.exists("./data")){
#  dir.create("./data", showWarnings = TRUE, recursive = TRUE)
#}
#if (!dir.exists("./data/results")){
#  dir.create("./data/results", showWarnings = TRUE, recursive = TRUE)
#}
if (!dir.exists(dirname(argv$outdir))){
  dir.create(dirname(argv$outdir), showWarnings = TRUE, recursive = TRUE)
}
if (!dir.exists("./data/spot_coords")){
  dir.create("./data/spot_coords", showWarnings = TRUE, recursive = TRUE)
}

print("Dir creation done.")



#import data
print("---Importing data---")
sc_seurat_meta_sce_split <- import_data_meta(data.dir = argv$scdata, 
                                       gene.column= argv$gene_column,
                                       project = "scRNA_humanbreastcancer",
                                       min.cells = 3,
                                       min.features = 200,
                                       meta.dir = argv$scmeta,
                                       grain_level = argv$grain_lvl,
                                       subtype = argv$subtype)
print("Import done.")


#generate the synth spots
print("---Generating synthetic spots---")
synthetic_visium_data <- generate_synthetic_visium_multi(seurat_obj = sc_seurat_meta_sce_split$seurat_obj_synth,
                                                         dataset_type = argv$synth_dataset, 
                                                         clust_var = "celltype_subset", 
                                                         n_regions = 5, 
                                                         max_n_region_spots = 175,
                                                         visium_mean = 30000, 
                                                         visium_sd = 8000,
                                                         select_celltype = argv$select_celltype,
                                                         n_cells_max = argv$n_cells_max,
                                                         min_cell_id_test = argv$min_cell_id_test,
                                                         select_celltype_min_id = argv$select_celltype_min_id)
print("Synthspot done.")

#get region coords
print("---Import region coords---")
selected_coords <- get_region_coords(region_file = argv$coords,
                                     total_file = argv$coords_total,
                                     export_file = "./data/spot_coords/regions_coords.csv")
print("Region coords import done.")

#build spatial obj with synth spots and region coords
print("---Building spatial objects---")
synth_spatialexp_counts <- assign_spotcoords_build_spatialexp(region_coords_file = "./data/spot_coords/regions_coords.csv",
                                                  synthetic_visium_data = synthetic_visium_data)
print("Spatial object creation done.")

#convert spatialexp to seurat spatial
print("---Converting st data format---")
spatial_obj_seurat <- spatialexp_to_seurat(spatial_obj = synth_spatialexp_counts$spatial_obj)
##SPATIAL PREPPED FOR DECONV
print("st data conversion done.")


if (argv$downsize >0){

  #downsize sc (overwrites current sc var)
  print("---Downsizing---")
  downsized_sc <- downsize_seurat_sce(sc_obj_seurat = sc_seurat_meta_sce_split$seurat_obj_train,
                                      meta = sc_seurat_meta_sce_split$meta_train,
                                      downsample = argv$downsize)
  sc_seurat_meta_sce_split$seurat_obj_train <- downsized_sc$seurat_obj
  sc_seurat_meta_sce_split$meta_train <- downsized_sc$meta
  sc_seurat_meta_sce_split$sce_obj_train <- downsized_sc$sce_obj
  print("DOWNSIZE FOR TESTING DONE. REMOVE IF NOT TESTING")
} 

# saveRDS(sce_downsize, file = "./data/rds/sce_downsize.rds")
# saveRDS(spatial_obj, file = "./data/rds/spatial_obj.rds")
# saveRDS(synthetic_visium_data, file = "./data/rds/synthetic_visium_data.RDS")
# saveRDS(region_coords, file = "./data/rds/spot_coords.RDS")
# saveRDS(filtered_counts, file = "./data/rds/filtered_counts.RDS")
# saveRDS(sc_obj_seurat, file = "./data/rds/sc_obj_seurat.rds")

##DECONVOLUTE
#also measure runtimes
print("---Deconvoluting---")
rctd_time_1 <- Sys.time()
deconv_rctd <- build_and_deconvolute(
  sc_seurat_meta_sce_split$sce_obj_train,
  synth_spatialexp_counts$spatial_obj,
  method = "rctd",
  cell_type_col = "celltype_subset",
  batch_id_col = NULL,
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE,
  n_cores = 8
)
rctd_time_2 <- Sys.time()
deconv_rctd <- true_celltype_colnames(deconv_rctd)

spotlight_time_1 <- Sys.time()
deconv_spotlight <-build_and_deconvolute(
  sc_seurat_meta_sce_split$sce_obj_train,
  synth_spatialexp_counts$spatial_obj,
  method = "spotlight",
  cell_type_col = "celltype_subset",
  batch_id_col = NULL,
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE
)
spotlight_time_2 <- Sys.time()
deconv_spotlight <- true_celltype_colnames(deconv_spotlight)

card_time_1 <- Sys.time()
deconv_card <-build_and_deconvolute(
  sc_seurat_meta_sce_split$sce_obj_train,
  synth_spatialexp_counts$spatial_obj,
  method = "card",
  cell_type_col = "celltype_subset",
  batch_id_col = "orig.ident",
  assay_sc = "counts",
  assay_sp = "counts",
  return_object = FALSE,
  verbose = TRUE
)
card_time_2 <- Sys.time()
deconv_card <- true_celltype_colnames(deconv_card)

print("Deconvolution done.")

##GENERATE STATS
print("---Generating stats---")
method_names <- c("rctd", "spotlight", "card")
methods <- list(deconv_rctd, deconv_spotlight, deconv_card)

#get rmsd
rmsd_all <- list(rmsd_table = data.frame(),
                 rmsd_table_mintest = data.frame())
i <- 1
for (method in methods){
  method = data.frame(method)
  rmsd_method <- getRMSD(prediction_fracs = method,
                  synthetic_visium_data = synthetic_visium_data,
                  method_annot = method_names[i],
                  min_test = argv$min_cell_id_test)
  rmsd_all$rmsd_table <- rbind(rmsd_all$rmsd_table, rmsd_method$rmsd_table)
  rmsd_all$rmsd_table_mintest <- rbind(rmsd_all$rmsd_table_mintest, rmsd_method$rmsd_table_mintest)
  
  i <- i+1
}
write.csv(x = rmsd_all$rmsd_table, file = paste0(argv$outdir, "rmsd.csv"))
if(dim(rmsd_all$rmsd_table_mintest)[1] != 0){
  write.csv(x = rmsd_all$rmsd_table_mintest, file = paste0(argv$outdir, "rmsd_mintest.csv"))
}
print("RMSD done.")

#get JSD
i <- 1
jsd_all <- list()
for (method in methods){
  method <- data.frame(method)
  jsd_method <- getJSD(prediction_fracs = method,
                       synthetic_visium_data = synthetic_visium_data,
                       method_annot = method_names[i],
                       min_test = argv$min_cell_id_test)
  jsd_all$mean <- c(jsd_all$mean, jsd_method$mean)
  jsd_all$jsd_table <- rbind(jsd_all$jsd_table, jsd_method$jsd_table)
  jsd_all$jsd_mean_mintest <- c(jsd_all$jsd_mean_mintest, jsd_method$mean)
  jsd_all$jsd_table_mintest <- rbind(jsd_all$jsd_table_mintest, jsd_method$jsd_table_mintest)
  i <- i+1
}
write.csv(x = jsd_all$jsd_table, file = paste0(argv$outdir, "jsd.csv"))
if(dim(rmsd_all$rmsd_table_mintest)[1] != 0){
  write.csv(x = jsd_all$jsd_table_mintest, file = paste0(argv$outdir, "jsd_mintest.csv"))
}
print("JSD done.")



##VISUALIZE
print("---Generating plots---")
#plot spatial scatter pie plot
i <- 1
for (method in methods){
  plot_spatial_scatter_pie(prediction_fracs = method,
                           selected_coords = selected_coords,
                           outfile= paste0( argv$outdir, "_spatial_scatterpie_", method_names[i], ".pdf" ),
                           scatterpie_alpha = 1,
                           pie_scale = 1)
  i <- i+1
}
#plot ground truth
plot_spatial_scatter_pie_truth(synthetic_visium_data = synthetic_visium_data,
                               selected_coords = selected_coords,
                               outfile= paste0(argv$outdir, "_truth_spatial_scatterpie.pdf"),
                               scatterpie_alpha = 1,
                               pie_scale = 1)
print("Spatial scatter pie done.")

#error heatmaps
#overall error heatmap

i <- 1
for (method in methods){
  pdf(paste0(argv$outdir, "rmsd_heatmap_ALL", method_names[i]  , ".pdf"))
  print(spatial_deconvolution_error_plot(synthetic_visium_data = synthetic_visium_data,
                                   predicted = method,
                                   coordinates = selected_coords))
  dev.off()
  i <- i+1
}

#per celltype error heatmaps
i <- 1
for (method in methods){
  spatial_deconvolution_error_heatmap(synthetic_visium_data = synthetic_visium_data,
                                      predicted = method,
                                      coordinates = selected_coords,
                                      outdir = paste0(argv$outdir, "rmsd_heatmap_", method_names[i], ".pdf"))
  i <- i+1
}

#final runtime
end_time <- Sys.time()
#get runtimes
runtime_vector <- c(end_time - start_time,
                    rctd_time_2 - rctd_time_1,
                    spotlight_time_2 - spotlight_time_1,
                    card_time_2 - card_time_1)
runtimes <- data.frame(method = c("benchdeconv", method_names),
                       runtimes = runtime_vector)
write.csv(runtimes, file = paste0(argv$outdir, "_runtimes.csv"))


print("------benchdeconv DONE------")



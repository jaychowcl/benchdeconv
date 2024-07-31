#!/exports/eddie/scratch/s2600569/envs/new_conda/bin/R

print("------STARTING benchdeconv------")
print(getwd())


#argparser args
print("---Parsing input args---")
library(argparser)
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

input_args <- add_argument(input_args, "--test", help = "Flag for running with reduced dataset. FOR TESTING ONLY. (0) off ,(1) on",
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

# Parse the command line arguments
argv <- parse_args(input_args)
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
if (!dir.exists(argv$outdir)){
  dir.create(argv$outdir, showWarnings = TRUE, recursive = TRUE)
}
print("Dir creation done.")



#import data
print("---Importing data---")
sc_seurat_meta <- import_data_meta(data.dir = argv$scdata, 
                                       gene.column= argv$gene_column,
                                       project = "scRNA_humanbreastcancer",
                                       min.cells = 3,
                                       min.features = 200,
                                       meta.dir = argv$scmeta,
                                       grain_level = argv$grain_lvl)
print("Import done.")

#split data for training and synthetic spot generation
print("---Splitting data---")
sc_seurat_meta_sce_split <- split_data(sc_obj_seurat = sc_seurat_meta$sc_obj_seurat,
                                   meta = sc_seurat_meta$meta,
                                   proportion = 0.5)

print("Data split done.")

#generate the synth spots
print("---Generating synthetic spots---")
synthetic_visium_data <- generate_synthetic_visium_multi(seurat_obj = sc_seurat_meta_sce_split$seurat_obj_synth,
                                                         dataset_type = argv$synth_dataset, 
                                                         clust_var = "celltype_subset", 
                                                         n_regions = 5, 
                                                         max_n_region_spots = 175,
                                                         visium_mean = 30000, 
                                                         visium_sd = 8000,
                                                         select_celltype = "T-cells")
print("Synthspot done.")

#get region coords
print("---Import region coords---")
selected_coords <- get_region_coords(region_file = "./data/spot_coords/out1.csv",
                                     total_file = "./data/spot_coords/Spatial-Projection.csv",
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


if (argv$test == 1){

##DOWNSIZING FOR TESTING
#downsize sc for testing (overwrites current sc var)
print("---DOWNSIZING. REMOVE BEFORE PROD.---")
downsized_sc <- downsize_seurat_sce(sc_obj_seurat = sc_seurat_meta_sce_split$seurat_obj_train,
                                    meta = sc_seurat_meta_sce_split$meta_train,
                                    downsample = 500)
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
print("---Deconvoluting---")
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
deconv_rctd <- true_celltype_colnames(deconv_rctd)

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
deconv_spotlight <- true_celltype_colnames(deconv_spotlight)

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
deconv_card <- true_celltype_colnames(deconv_card)

print("Deconvolution done.")

##GENERATE STATS
print("---Generating stats---")
method_names <- c("rctd", "spotlight", "card")
methods <- list(deconv_rctd, deconv_spotlight, deconv_card)

#get rmsd
rmsd_all <- data.frame()
i <- 1
for (method in methods){
  method = data.frame(method)
  rmsd_method <- getRMSD(prediction_fracs = method,
                  synthetic_visium_data = synthetic_visium_data,
                  method_annot = method_names[i])
  rmsd_all <- rbind(rmsd_all, rmsd_method)
  i <- i+1
}
write.csv(x = rmsd_all, file = paste0(argv$outdir, "rmsd.csv"))
print("RMSD done.")

#get JSD
i <- 1
jsd_all <- list()
for (method in methods){
  method <- data.frame(method)
  jsd_method <- getJSD(prediction_fracs = method,
                       synthetic_visium_data = synthetic_visium_data,
                       method_annot = method_names[i])
  jsd_all$mean <- c(jsd_all$mean, jsd_method$mean)
  jsd_all$jsd_table <- rbind(jsd_all$jsd_table, jsd_method$jsd_table)
  i <- i+1
}
write.csv(x = jsd_all, file = paste0(argv$outdir, "jsd.csv"))
print("JSD done.")



##VISUALIZE
print("---Generating plots---")
#plot spatial scatter pie plot
i <- 1
for (method in methods){
  plot_spatial_scatter_pie(prediction_fracs = method,
                           selected_coords = selected_coords,
                           outfile= paste0( argv$outdir, "spatial_scatterpie_", method_names[i], ".pdf" ),
                           scatterpie_alpha = 1,
                           pie_scale = 0.4)
  i <- i+1
}
#plot ground truth
plot_spatial_scatter_pie_truth(synthetic_visium_data = synthetic_visium_data,
                               selected_coords = selected_coords,
                               outfile= paste0(argv$outdir, "/truth_spatial_scatterpie.pdf"),
                               scatterpie_alpha = 1,
                               pie_scale = 0.4)
print("Spatial scatter pie done.")

#


print("------SHUTTING DOWN benchdeconv------")



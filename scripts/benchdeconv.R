library(DropletUtils)
library(devtools)
library(synthspot)
library(Seurat)
library(dplyr)
library(patchwork)
library(SpatialExperiment)
library(spacedeconv)
library(Matrix)
library(ggplot2)
library(SPOTlight)
library(philentropy)



# Read sc data
import_data_meta <- function(data.dir = "data/scRNA_wu", 
                             gene.column=1,
                             project = "scRNA_humanbreastcancer",
                             min.cells = 3,
                             min.features = 200,
                             meta.dir = "/localdisk/home/s2600569/project/data/scRNA_wu/metadata.csv",
                             grain_level = "celltype_major"
                             ){
  # Import data
  sc_obj_seurat.data <- Read10X(data.dir = data.dir, gene.column = gene.column)
  sc_obj_seurat <- CreateSeuratObject(counts = sc_obj_seurat.data, project = project, min.cells = min.cells, min.features = min.features)
  
  # Import metadata
  meta <- read.csv(meta.dir)
  
  # Format and add metadata
  rownames(meta) <- colnames(sc_obj_seurat)
  AddMetaData(sc_obj_seurat, meta)
  sc_obj_seurat@meta.data$celltype_subset <- meta[[grain_level]]
  
  # Set idents
  Idents(sc_obj_seurat) <- meta[[grain_level]]
  
  return(list(sc_obj_seurat = sc_obj_seurat, meta = meta))
}



# Split sc data
split_data <- function(sc_obj_seurat = sc_seurat_meta$sc_obj_seurat,
                       meta = sc_seurat_meta$meta,
                       proportion = 0.5){
  # Ensure the proportion is between 0 and 1
  if (proportion <= 0 | proportion >= 1) {
    stop("Proportion must be between 0 and 1 (exclusive).")
  }
  
  # Get the total number of cells
  total_cells <- ncol(sc_obj_seurat)
  
  # Calculate the number of cells for the first subset
  num_cells_1 <- round(total_cells * proportion)
  
  # Generate a random sample of cells for the first subset
  #set.seed(seed) # Setting a seed for reproducibility
  cell_indices <- sample(1:total_cells, size = num_cells_1)
  
  # Split the cells into two subsets
  cells_1 <- colnames(sc_obj_seurat)[cell_indices]
  cells_2 <- colnames(sc_obj_seurat)[-cell_indices]
  
  # Create two new Seurat objects
  seurat_obj_1 <- subset(sc_obj_seurat, cells = cells_1)
  seurat_obj_2 <- subset(sc_obj_seurat, cells = cells_2)
  
  # Split metadata
  meta_1 <- subset(meta, X %in% cells_1)
  meta_2 <- subset(meta, X %in% cells_2)
  
  # Convert to sce
  sce_obj_synth = as.SingleCellExperiment(seurat_obj_1, assay = "RNA")
  sce_obj_train = as.SingleCellExperiment(seurat_obj_2, assay = "RNA")
  
  return(list(seurat_obj_synth = seurat_obj_1,
              seurat_obj_train = seurat_obj_2,
              meta_synth = meta_1,
              meta_train = meta_2,
              sce_obj_synth = sce_obj_synth,
              sce_obj_train = sce_obj_synth))
}

# Run synthspots multiple times to fill st regions
generate_synthetic_visium_multi <- function(seurat_obj = sc_seurat_meta_split$seurat_obj_synth,
                                            dataset_type = "artificial_partially_dominant_celltype_diverse", 
                                            clust_var = "celltype_subset", 
                                            n_regions = 5, 
                                            max_n_region_spots = 175,
                                            visium_mean = 30000, 
                                            visium_sd = 8000,
                                            select_celltype = "random"){
  # Find no. of repeats to run
  freq_cells <- table(seurat_obj@meta.data[,clust_var]) # table of freqs
  min_celltype <- min(seurat_obj@meta.data[, clust_var]) # celltype with least cells
  cycles <- ceiling(max_n_region_spots / freq_cells[min_celltype]) # no. of repeats
  
  # If one cycle:
  if(cycles == 1){
    synthetic_visium_data <- generate_synthetic_visium(seurat_obj = seurat_obj, 
                                                       dataset_type = dataset_type, 
                                                       clust_var = clust_var, 
                                                       n_regions = n_regions, 
                                                       n_spots_min = freq_cells[min_celltype]-1,
                                                       n_spots_max = freq_cells[min_celltype],
                                                       visium_mean = visium_mean, 
                                                       visium_sd = visium_sd,
                                                       select_celltype = select_celltype)
  } 
  
  # If multiple cycles:
  if (cycles >1){
    synthetic_visium_data <- generate_synthetic_visium(seurat_obj = seurat_obj, 
                                                       dataset_type = dataset_type, 
                                                       clust_var = clust_var, 
                                                       n_regions = n_regions, 
                                                       n_spots_min = freq_cells[min_celltype]-1,
                                                       n_spots_max = freq_cells[min_celltype],
                                                       visium_mean = visium_mean, 
                                                       visium_sd = visium_sd,
                                                       select_celltype = select_celltype)
    for (i in 1:cycles-1){
      synthetic_visium_data_add <- generate_synthetic_visium(seurat_obj = seurat_obj, 
                                                         dataset_type = dataset_type, 
                                                         clust_var = clust_var, 
                                                         n_regions = n_regions, 
                                                         n_spots_min = freq_cells[min_celltype]-1,
                                                         n_spots_max = freq_cells[min_celltype],
                                                         visium_mean = visium_mean, 
                                                         visium_sd = visium_sd,
                                                         select_celltype = select_celltype)
      
    }
  }
  return(synthetic_visium_data)
}

#
get_region_coords <- function(region_file = "./data/spot_coords/out1.csv",
                              total_file = "./data/spot_coords/Spatial-Projection.csv",
                              export_file = "./data/spot_coords/regions_coords.csv"){
  # import data
  selected_spots <- read.csv(region_file)
  total_spots <- read.csv(total_file)
  
  # Select selected spot coords
  selected_coords <- subset(total_spots, Barcode %in% selected_spots$Barcode) # get coords of selected spots
  selected_coords <- cbind(selected_spots$out1, selected_coords) # add in region names to coords
  selected_coords <- subset(selected_coords, select = -Barcode) # remove temp barcodes
  selected_coords <- selected_coords[order(selected_coords$`selected_spots$out1`),] # group regions together
  
  # rename region names to spot names
  index =1
  name = 1
  current <- selected_coords$`selected_spots$out1`[1]
  for (spot in selected_coords$`selected_spots$out1`){
    if (spot != current ){
      current <- spot
      name = 1
      selected_coords$`selected_spots$out1`[index] <- paste0(spot, "_spot_", name)
    } else {
      selected_coords$`selected_spots$out1`[index] <- paste0(spot, "_spot_", name)
    }
    index = index + 1
    name = name + 1
    current <- spot
  }
  # export
  write.csv(selected_coords, file = export_file, row.names = FALSE)
  
  return(selected_coords)
}

# Assign coords to spots and build spatialexperiment object

assign_spotcoords_build_spatialexp <- function(region_coords_file = "./data/spot_coords/regions_coords.csv",
                                               synthetic_visium_data = synthetic_visium_data
                                               ){
  # import spot coordinates from LoupeR spatial barcodes
  region_coords <- read.csv(region_coords_file)
  rownames(region_coords) <- region_coords$selected_spots.out1
  
  # extract counts from synthspot
  counts <- synthetic_visium_data$counts
  
  # assign coords to synthspot spots
  col_counts <- colnames(counts)
  col_regions <- as.vector(region_coords$selected_spots.out1)
  choose_counts <- col_counts %in% col_regions
  filtered_counts <- counts[, choose_counts] # final counts matrix
  
  # build SPATIALEXPERIMENT spatial object
  spatial_obj <- SpatialExperiment(assay= list(counts = filtered_counts),
                                   colData=region_coords,
                                   spatialCoordsNames = c("X.Coordinate", "Y.Coordinate"))
  return(list(spatial_obj = spatial_obj,
              filtered_counts = filtered_counts))
  
}

#convert to spatial exp to seurat
spatialexp_to_seurat <- function(spatial_obj = synth_spatialexp_counts$spatial_obj){
  #extract from spatial_obj
  counts_seurat <- counts(spatial_obj)
  metadata_seurat <- colData(spatial_obj)
  spatial_coords_seurat <- spatialCoords(spatial_obj)
  colnames(spatial_coords_seurat) <- c("spatial_1","spatial_2")
  #create Seurat object
  spatial_obj_seurat <- CreateSeuratObject(counts = counts_seurat)
  #add metadata
  spatial_obj_seurat <- AddMetaData(object = spatial_obj_seurat, metadata = as.data.frame(metadata_seurat))
  #add spatial coordinates
  spatial_obj_seurat[['spatial']] <- CreateDimReducObject(embeddings = as.matrix(spatial_coords_seurat), key = "spatial_", assay = "RNA")
  
  return(spatial_obj_seurat)
}

# FOR TESTING: donwsize single cell so it runs faster (MAINLY FOR DEBUGGING AND DEVELOPMENT)
downsize_seurat_sce <- function(sc_obj_seurat = sc_seurat_meta_split$seurat_obj_train,
                                meta = sc_seurat_meta_split$meta_train,
                                downsample = 500){
  seurat_downsize <- subset(sc_obj_seurat, downsample = downsample) 
  sce_downsize <- as.SingleCellExperiment(seurat_downsize)
  
  #downsize metadata as well
  meta_downsize <- meta[rownames(meta) %in% rownames(sce_downsize@colData),]
  
  return(list(seurat_obj = seurat_downsize,
              sce_obj = sce_downsize,
              meta = meta_downsize))
}


##CALCULATING STATISTICS
#clean spacedeconv cols for rmsd
true_celltype_colnames <- function(prediction_table){
  cols <- colnames(prediction_table)
  cols <- sub(".*_", "", cols)
  colnames(prediction_table) <- cols
  return(prediction_table)
}

#calculate rmsd
getRMSD <- function(prediction_fracs = deconv_rctd,
                    synthetic_visium_data = synthetic_visium_data,
                    method_annot = "deconv_tool"){
  #get true fractions
  true_fracs <- subset(synthetic_visium_data$relative_spot_composition, select = -region)
  rownames(true_fracs) <- synthetic_visium_data$relative_spot_composition$name
  #Choose only included spots
  true_fracs <- subset(true_fracs, name %in% rownames(prediction_fracs), select = -name)
  
  #init rmsd vector
  rmsd <- 1:length(colnames(true_fracs))
  
  #calc rmsd
  for (i2 in 1:length(colnames(prediction_fracs))){#per celltype
    #print(i2)
    
    #init differences
    differences_squared <- 1:length(rownames(true_fracs))
    
    for (i in 1:length(rownames(prediction_fracs))){#per spot
      #print(i)
      difference <- true_fracs[i, i2] - prediction_fracs[i, i2]
      difference_squared <- difference^2
      differences_squared[i] <- difference_squared
    }
    
    mean_sum <- sum(differences_squared) / length(rownames(true_fracs))
    rmsd_celltype <- sqrt(mean_sum)
    rmsd[i2] <- rmsd_celltype
  }
  
  method_annot <- rep(method_annot, length(colnames(true_fracs)))
  rmsd_table <- data.frame(method = method_annot, celltype = colnames(prediction_fracs), rmsd = rmsd)
  return(rmsd_table)
}

#calculate JSD
getJSD <- function(prediction_fracs = deconv_rctd,
                   synthetic_visium_data = synthetic_visium_data,
                   method_annot = "deconv_tool"){
  #get true fractions
  true_fracs <- subset(synthetic_visium_data$relative_spot_composition, select = -region)
  rownames(true_fracs) <- synthetic_visium_data$relative_spot_composition$name
  #Choose only included spots
  true_fracs <- subset(true_fracs, name %in% rownames(prediction_fracs), select = -name)
  
  #jsd calc function
  calculate_jsd <- function(p, q) {
  M <- 0.5 * sum((p + q))
  jsd <- 0.5 * (sum(p * log(sum(p) / sum(M))) + sum(sum(q) * log(sum(q) / M)))
  return(jsd)
  }
  
  #get all jsd for each spot
  jsd_all <- seq(1:dim(true_fracs)[2])
  for (i in 1:dim(true_fracs)[2]){
    jsd_i <- calculate_jsd(true_fracs[, i], prediction_fracs[, i])
    jsd_all[i] <- jsd_i
  }
  
  #gather all jsd and calc jsd mean
  jsd_table <- data.frame(method = method_annot, celltype = colnames(true_fracs), jsd = jsd_all)
  jsd_mean <- mean(jsd_table$jsd)
  
  return(list(mean = jsd_mean,
              jsd_table = jsd_table))
}

##VISUALIZE
#plot spatial scatter pie
plot_spatial_scatter_pie <- function(prediction_fracs = deconv_rctd,
                                     selected_coords = selected_coords,
                                     outfile= "./data/results/rctd_spatial_scatterpie.pdf",
                                     scatterpie_alpha = 1,
                                     pie_scale = 0.4){
  rownames(selected_coords) <- selected_coords$`selected_spots$out1`
  spot_coords <- subset(selected_coords, select = -`selected_spots$out1`)
  # print(head(spot_coords))
  # print(head(prediction_fracs))
  
  spatialscatter <- plotSpatialScatterpie(
    x = spot_coords,
    y = prediction_fracs,
    cell_types = colnames(prediction_fracs),
    img = FALSE,
    scatterpie_alpha = scatterpie_alpha,
    pie_scale = pie_scale
  )
  
  # print(spatialscatter)
  pdf(outfile)
  print(spatialscatter)
  dev.off()
}

#visualize ground truth
plot_spatial_scatter_pie_truth <- function(synthetic_visium_data = synthetic_visium_data,
                                           selected_coords = selected_coords,
                                           outfile= "./data/results/truth_spatial_scatterpie.pdf",
                                           scatterpie_alpha = 1,
                                           pie_scale = 0.4){
  #prep rownames of selected coords
  rownames(selected_coords) <- selected_coords$`selected_spots$out1`
  spot_coords <- subset(selected_coords, select = -`selected_spots$out1`)
  #Extract truth and choose selected spots
  relative_spots_truth <- synthetic_visium_data$relative_spot_composition
  rownames(relative_spots_truth) <- relative_spots_truth$name
  relative_spots_truth <- relative_spots_truth[relative_spots_truth$name %in% rownames(spot_coords), ]
  
  spatialscatter <- plotSpatialScatterpie(
    x = spot_coords,
    y = relative_spots_truth,
    cell_types = colnames(subset(relative_spots_truth, select = c(-name, -region))),
    img = FALSE,
    scatterpie_alpha = 1,
    pie_scale = 0.4
  )
  
  pdf(outfile)
  print(spatialscatter)
  dev.off()
}

#gather rmsd and jsd stats 

outdir <- "data/results/all_runs/run_1"
gather_and_append_all_stats <- function(outdir = argv$outdir,
                                        finalfile = paste0(argv$outdir, "jsd_rmsd_all.csv")){
  infile_jsd <- read.csv(paste0(outdir, "jsd.csv"))
  infile_rmsd <- read.csv(paste0(outdir,"rmsd.csv"))
  
  
  
}

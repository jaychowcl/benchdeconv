library(ggplot2)

#gather runs
gather_runs <- function(indir = "./data/all_runs/",
                        prefix = "run_",
                        runs = c(1,50)){
  
  #gather rmsd and jsd stats from runs
  all_run_files <- c()
  for(i in runs[1]:runs[2]){
    runs_files <- list.files(path=indir,
                             pattern = paste0(prefix, i, "[rj]"),
                             full.names = TRUE)
    all_run_files <- c(all_run_files, runs_files)
                             
  }  
  
  if(length(all_run_files) == 0){
    stop("No files in indir!")
  }
  
    
  rmsd_files <- all_run_files[grep("rmsd.csv$", all_run_files)]
  jsd_files <- all_run_files[grep("jsd.csv$", all_run_files)]

  #read all files and put into single table
  final_tab <- data.frame()
  for(i in 1:length(rmsd_files)){
    #read files
    infile_rmsd <- read.csv(paste0(rmsd_files[i]))
    infile_jsd <- read.csv(paste0(jsd_files[i]))
    
    #put into table
    run_no <- regmatches(rmsd_files[i], regexpr("run_(\\d+)", rmsd_files[i]))
    new_tab <- data.frame(method = infile_rmsd$method,
                          celltype = infile_rmsd$celltype,
                          rmsd = infile_rmsd$rmsd,
                          jsd = infile_jsd$jsd,
                          run = run_no)
    #append
    final_tab <- rbind(final_tab, new_tab)
  }
  
  return(final_tab)
}

#gather mintest runs
gather_runs_mintest <- function(indir = "./data/all_runs/",
                        prefix = "run_",
                        runs = c(1,50)){
  
  #gather rmsd and jsd stats from runs
  all_run_files <- c()
  for(i in runs[1]:runs[2]){
    runs_files <- list.files(path=indir,
                             pattern = paste0(prefix, i, "[rj]"),
                             full.names = TRUE)
    all_run_files <- c(all_run_files, runs_files)
    
  }  
  
  if(length(all_run_files) == 0){
    stop("No files in indir!")
  }
  
  
  rmsd_files <- all_run_files[grep("rmsd_mintest.csv$", all_run_files)]
  jsd_files <- all_run_files[grep("jsd_mintest.csv$", all_run_files)]
  
  #read all files and put into single table
  final_tab <- data.frame()
  for(i in 1:length(rmsd_files)){
    #read files
    infile_rmsd <- read.csv(paste0(rmsd_files[i]))
    infile_jsd <- read.csv(paste0(jsd_files[i]))
    
    #put into table
    run_no <- regmatches(rmsd_files[i], regexpr("run_(\\d+)", rmsd_files[i]))
    new_tab <- data.frame(method = infile_rmsd$method,
                          celltype = infile_rmsd$celltype,
                          rmsd = infile_rmsd$rmsd,
                          jsd = infile_jsd$jsd,
                          run = run_no)
    #append
    final_tab <- rbind(final_tab, new_tab)
  }
  
  return(final_tab)
}


#plot boxplot
plot_boxplot <- function(indata = datasets_data,
                         plot_metric = "rmsd",
                         annot = "HER2+",
                         annot_col = "annot",
                         annot2 = "none",
                         annot2_col = "none",
                         cats = "method",
                         groups = "celltype"
                         ){
  if(annot != "all"){
    indata <- indata[indata[[annot_col]] == annot, ]
  }
  
  if(annot2 != "none"){
    indata <- indata[indata[[annot2_col]] == annot2, ]
  }
  
  # Create the horizontal boxplot for the selected plot_metric
  boxplot(indata[[plot_metric]] ~ indata[[cats]],
          horizontal = TRUE,
          xlab = plot_metric,
          ylab = "",
          main = paste(toupper(plot_metric), "Values by ", cats),
          las = 1,
          outline = FALSE)
  
  if(groups != "none"){
    #prep points
    points_tab <- data.frame()
    i=1
    for(cat in unique(indata[[cats]])){
      for(group in unique(indata[[groups]])){
        subset_values <- indata[[plot_metric]][indata[[cats]] == cat & indata[[groups]] == group]
        subset_value <- mean(subset_values)
        stderr <- sd(subset_values) / sqrt(length(subset_values))
        
        subset_table <- data.frame(method = cat,
                                   celltype = group,
                                   metric = subset_value,
                                   stderr = stderr)
        
        points_tab <- rbind(points_tab, subset_table)
      }
    }
    
    # Define unique cell types and assign symbols and colors
    cell_types <- unique(points_tab$celltype)
    symbols <- 1:length(cell_types) # Using different point symbols (1 to number of cell types)
    colors <- rainbow(length(cell_types)) # Assign different colors to each cell type
    
    # Create a mapping of cell types to symbols and colors
    celltype_symbol_map <- setNames(symbols, cell_types)
    celltype_color_map <- setNames(colors, cell_types)
    
    # Add  points with different symbols and colors for each cell type
    for (i in 1:nrow(points_tab)) {
      points(points_tab[["metric"]][i], as.numeric(factor(points_tab$method))[i], 
             pch = celltype_symbol_map[points_tab$celltype[i]], 
             col = celltype_color_map[points_tab$celltype[i]])
      
    }
  }
  # Add a smaller legend
  legend("right", legend = cell_types, pch = symbols, col = colors, title = groups, cex = 0.7)
  
}

#plot scatter

plot_scatter <- function(indata = mintest_data,
                         plot_metric = "rmsd",
                         method = "rctd",
                         annot = "B-cells",
                         annot_col = "select_celltype"){
  
  # indata <- indata[indata[[annot_col]] == annot, ]
  indata_all <- indata[indata$method == method, ]
  indata <- indata_all
  
  #get means
  mean_tab <- c()
  for(select_celltype in unique(indata$select_celltype)){
    for(density in unique(indata$density)){
      for(celltype in unique(indata$celltype)){
        
        # print(indata$rmsd[(indata$select_celltype == select_celltype) & (indata$density == density) & (indata$celltype == celltype)])
        rmsd_values <- indata$rmsd[(indata$select_celltype == select_celltype) & (indata$density == density) & (indata$celltype == celltype)]
        rmsd <- mean(rmsd_values,
                     na.rm = TRUE)
        rmsd_se <- sd(rmsd_values)/sqrt(length(rmsd_values))
        
        jsd_values <- indata$jsd[(indata$select_celltype == select_celltype) & (indata$density == density) & (indata$celltype == celltype)]
        jsd <-mean(jsd_values,
                   na.rm = TRUE)
        jsd_se <- sd(jsd_values)/sqrt(length(jsd_values))
       
        
        df_append <- data.frame(celltype = celltype,
                                rmsd = rmsd,
                                rmsd_se = rmsd_se,
                                jsd = jsd,
                                jsd_se = jsd_se,
                                density = density,
                                select_celltype = select_celltype) 
        mean_tab <- rbind(mean_tab, df_append)
      }
    }
  }
  mean_tab$select_celltype <- gsub("-", ".", mean_tab$select_celltype)
  #get only select_celltype
  mean_tab <- mean_tab[mean_tab$celltype == mean_tab$select_celltype,]
  
  # Create the plot
  ggplot(mean_tab, aes(x = density, y = rmsd, color = select_celltype)) +
    geom_line() +            # Line plot
    geom_point() +           # Add points
    geom_errorbar(aes(ymin = mean_tab$rmsd - mean_tab$rmsd_se, ymax = mean_tab$rmsd + mean_tab$rmsd_se), 
                  width = 0.3, 
                  color = "black", 
                  size = 0.8) +
    labs(title = "Multiple Groups Line Plot",
         x = "Density",
         y = plot_metric,
         color = "Group")   # Labels and title
  
  # 
  # plot(x = mean_tab$density[mean_tab$select_celltype == mean_tab$select_celltype[1]],
  #      y = mean_tab[[plot_metric]][mean_tab$select_celltype == mean_tab$select_celltype[1]],
  #      type = "b",
  #      xlab = "Density",
  #      ylab = plot_metric,
  #      pch=1,
  #      col = 1)
  # 
  # for(celltype in unique(mean_tab$celltype)[-1]){
  #   points()
  #   
  # }
  # 
  
}





###START
stats_args <- input_args <- arg_parser("run_stats.R: A tool for generating further stats after running benchdeconv through high throughput computing.")

stats_args <- add_argument(stats_args, "--indir", help="Input benchdeconv results directory.",
                           default = "./data/all_runs")
stats_args <- add_argument(stats_args, "--outdir", help="Output directory.",
                           default = "./data/results/")
argv <- parse_args(stats_args)


#experiments:
#t1-50 insize test
insize_annot_tags <- c(rep("none", 10),
                       rep(500, 10),
                       rep(1000, 10),
                       rep (2500, 10),
                       rep(5000, 10))
#gather indata 
insize_data <- gather_runs(indir = argv$indir,
                           prefix = "run_",
                           runs=c(1,50))






#t50-80 differ datasets
#gather data
datasets_data <- gather_runs(indir = argv$indir,
                           prefix = "run_",
                           runs=c(51,80))
dataset_annot_tags <- c(rep("HER2+", length(datasets_data$run[datasets_data$run == unique(datasets_data$run)[1]])*10),
                        rep("ER+", length(datasets_data$run[datasets_data$run == unique(datasets_data$run)[1]])*10),
                        rep("TNBC", length(datasets_data$run[datasets_data$run == unique(datasets_data$run)[1]])*10) )
                          
datasets_data$annot <- dataset_annot_tags


#plot methods per dataset
for(dataset_type in unique(dataset_annot_tags)){
  
  pdf(paste0(argv$outdir, dataset_type, "_datasets_boxplot_permethodpercelltype_rmsd.pdf"))
  plot_boxplot(indata = datasets_data,
               plot_metric = "rmsd",
               annot = dataset_type,
               annot_col = "annot",
               cats = "method",
               groups = "celltype")
  dev.off()
  
  pdf(paste0(argv$outdir, dataset_type, "_datasets_boxplot_permethodpercelltype_jsd.pdf"))
  plot_boxplot(indata = datasets_data,
               plot_metric = "jsd",
               annot = dataset_type,
               annot_col = "annot",
               cats = "method",
               groups = "celltype")
  dev.off()
  
}
#plot datasets
pdf(paste0(argv$outdir,"datasets_boxplot__perdatasetpermethod_rmsd.pdf"))
plot_boxplot(indata = datasets_data,
             plot_metric = "rmsd",
             annot = "all",
             annot_col = "annot",
             cats = "annot",
             groups = "method")
dev.off()

pdf(paste0(argv$outdir,"datasets_boxplot__perdatasetpermethod_jsd.pdf"))
plot_boxplot(indata = datasets_data,
             plot_metric = "jsd",
             annot = "all",
             annot_col = "annot",
             cats = "annot",
             groups = "method")
dev.off()
#plot celltypes 
for(method in unique(datasets_data$method)){
  pdf(paste0(argv$outdir, "datasets_boxplot_percelltpyepermethod_rmsd_",method, ".pdf"))
  plot_boxplot(indata = datasets_data,
               plot_metric = "rmsd",
               annot = "TNBC",
               annot_col = "annot",
               annot2 = method,
               annot2_col = "method",
               cats = "celltype",
               groups = "none")
  dev.off()
}




#t80-180 min density for id
#gather the data 
mintest_data <- gather_runs_mintest(indir = argv$indir,
                                    prefix = "run_",
                                    runs=c(81,180))

#gather annot
mintest_density_tags <- c(rep(0.05, length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*5),
                          rep(0.1, length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*5),
                          rep(0.2, length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*5),
                          rep(0.5, length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*5),
                          rep(0.8, length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*5))
mintest_data$density <- mintest_density_tags

mintest_celltype_tags <- c(rep("B-cells", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25),
                           rep("Plasmablasts", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25),
                           rep("T-cells", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25),
                           rep("Myeloid", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25))
mintest_data$select_celltype <- mintest_celltype_tags

#create scatter
pdf(paste0(argv$outdir, "mintest_scatter", ".pdf"))
plot_scatter(indata = mintest_data,
             plot_metric = "rmsd",
             method = "rctd",
             annot = "B-cells",
             annot_col = "select_celltype")
dev.off()



# insize_stats <- gather_runs(indir = "./data/results/all_runs",
#                                       prefix = "run",
#                                       runs = c(1,50),
#                                       annot_tags = insize_annot_tags)
# 
# 

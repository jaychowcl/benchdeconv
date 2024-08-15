library(ggplot2)
library(argparser)
library(dunn.test)
library(FSA)         # For Dunn test
library(ggsignif)    # For adding statistical annotations
library(dplyr)       # For data manipulation

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
                         groups = "celltype",
                         ggplot_type = "0"
                         ){
  if(annot != "all"){
    indata <- indata[indata[[annot_col]] == annot, ]
  }
  
  if(annot2 != "none"){
    indata <- indata[indata[[annot2_col]] == annot2, ]
  }
  
  if(ggplot_type == "0"){
    # pdf("testhere.pdf")
    # Create the horizontal boxplot for the selected plot_metric
    boxplot(indata[[plot_metric]] ~ indata[[cats]],
            horizontal = TRUE,
            xlab = plot_metric,
            ylab = "",
            main = paste(toupper(plot_metric), "Values by ", cats),
            las = 1,
            outline = FALSE,
            range = 0,
            ylim = c(0, max(indata[[plot_metric]])))
   # dev.off()
    
    if(groups != "none"){
      #prep points
      points_tab <- data.frame()
      i=1
      for(cat_i in unique(indata[[cats]])){
        for(group in unique(indata[[groups]])){
          subset_values <- indata[[plot_metric]][indata[[cats]] == cat_i & indata[[groups]] == group]
          subset_value <- mean(subset_values)
          stderr <- sd(subset_values) / sqrt(length(subset_values))
          
          subset_table <- data.frame(method = cat_i,
                                     celltype = group,
                                     metric = subset_value,
                                     stderr = stderr)
          
          points_tab <- rbind(points_tab, subset_table)
        }
      }
      
      # Define unique cell types and assign symbols and colors
      cell_types <- unique(points_tab$celltype)
      symbols_i <- 1:length(cell_types) # Using different point symbols (1 to number of cell types)
      colors_i <- rainbow(length(cell_types)) # Assign different colors to each cell type
      
      # Create a mapping of cell types to symbols and colors
      celltype_symbol_map <- setNames(symbols_i, cell_types)
      celltype_color_map <- setNames(colors_i, cell_types)
      
      # Add  points with different symbols and colors for each cell type
      for (i in 1:nrow(points_tab)) {
        points(points_tab[["metric"]][i], as.numeric(factor(points_tab$method))[i], 
               pch = celltype_symbol_map[points_tab$celltype[i]], 
               col = celltype_color_map[points_tab$celltype[i]])
        
      } 
      
    # Add a smaller legend
    legend("right", legend = cell_types, pch = symbols_i, col = colors_i, title = groups, cex = 0.7)
  
      
    } else {
      cell_types <- unique(indata$celltype)
    }
  } else if(ggplot_type != 0){
    print("doing ggplot")
    # Perform Dunn test
    dunn_results <- dunnTest(indata[[plot_metric]] ~ indata[[cats]], method = "bonferroni")
    
    # Extract the significant comparisons
    significant_comparisons <- dunn_results$res[dunn_results$res$P.adj < 0.05, c("Comparison", "P.adj")]
    
    # Create a boxplot with ggplot2 to allow for easy annotation
    p <- ggplot(indata, aes_string(x = cats, y = plot_metric)) +
      geom_boxplot(outlier.shape = NA) +
      coord_flip() +
      labs(x = "", y = plot_metric, title = paste(toupper(plot_metric), "Values by", cats)) +
      theme_classic()
    
    # Add Dunn test results to the plot
    if (nrow(significant_comparisons) > 0) {
      # Prepare the data for annotation
      comparisons <- strsplit(significant_comparisons$Comparison, " - ")
      
      for (i in seq_along(comparisons)) {
        p <- p + geom_signif(comparisons = list(comparisons[[i]]),
                             annotations = paste0("p = ", round(significant_comparisons$P.adj[i], 3)),
                             map_signif_level = FALSE,
                             y_position = max(indata[[plot_metric]]) * (1 + 0.05 * i))
      }
    }
    
    # Display the plot
    print(p)
    
    # Save the plot
    ggsave(filename = paste0(argv$outdir, "datasets_boxplot_percelltpyepermethod", method, metric, ".png"), plot = p, width = 12, height = 6)
    
    
    }
  
}

#plot scatter

plot_scatter_mintest <- function(indata = mintest_data,
                         method = "rctd"){
  
  # indata <- indata[indata[[annot_col]] == annot, ]
  #plot_metric = "rmsd", annot = "B-cells",  annot_col = "select_celltype"
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
  
  return(mean_tab)

}

#plot scatter insize
plot_scatter_insize <- function(indata = insize_data,
                                method = "rctd"){
  
  # indata <- indata[indata[[annot_col]] == annot, ]
  # #plot_metric = "rmsd", annot = "B-cells",  annot_col = "select_celltype"
  # indata_all <- indata[indata$method == method, ]
  # indata <- indata_all
  # 
  #get means
  mean_tab <- data.frame()
  for(method in unique(indata$method)){
    for(n_cell in unique(indata$n_cells)){
        
        # print(indata$rmsd[(indata$select_celltype == select_celltype) & (indata$density == density) & (indata$celltype == celltype)])
        rmsd_values <- na.omit(indata$rmsd[(indata$n_cells == n_cell) & (indata$method == method)])
        rmsd <- mean(rmsd_values,
                     na.rm = TRUE)
        rmsd_se <- sd(rmsd_values)/sqrt(length(rmsd_values))
        
        jsd_values <- na.omit(indata$jsd[(indata$n_cells == n_cell)& (indata$method == method)])
        jsd <-mean(jsd_values,
                   na.rm = TRUE)
        jsd_se <- sd(jsd_values)/sqrt(length(jsd_values))
        
        method_i <- method
        
        df_append <- data.frame(rmsd = rmsd,
                                rmsd_se = rmsd_se,
                                jsd = jsd,
                                jsd_se = jsd_se,
                                n_cells = n_cell,
                                method = method_i) 
        mean_tab <- rbind(mean_tab, df_append)
      
    }
  }
  # mean_tab$select_celltype <- gsub("-", ".", mean_tab$select_celltype)
  # #get only select_celltype
  # mean_tab <- mean_tab[mean_tab$celltype == mean_tab$select_celltype,]
  # 
  return(mean_tab)
  
}

#table of jsd and rmsd
create_pivot_dataframe <- function(data, value_column) {
  
  # Create the new dataframe
  new_df <- data %>%
    select(celltype, run, !!sym(value_column)) %>%  # Select the relevant columns dynamically
    spread(key = run, value = !!sym(value_column))  # Reshape data: spread 'run' as columns, dynamic value column
  
  # Set the row names as cell types
  rownames(new_df) <- new_df$celltype
  
  # Drop the celltype column since it's now row names
  new_df <- new_df %>% select(-celltype)
  
  # Return the new dataframe
  return(new_df)
}






###START
stats_args <- input_args <- arg_parser("run_stats.R: A tool for generating further stats after running benchdeconv through high throughput computing.")

stats_args <- add_argument(stats_args, "--indir", help="Input benchdeconv results directory.",
                           default = "./data/all_runs")
stats_args <- add_argument(stats_args, "--outdir", help="Output directory.",
                           default = "./data/results/")
argv <- parse_args(stats_args)

if (!dir.exists(argv$outdir)){
  dir.create(argv$outdir, showWarnings = TRUE, recursive = TRUE)
}

# argv$indir <- "./data/all_runs"
# argv$outdir <- "./data/results/"

#save all data
alldata <- gather_runs(indir = argv$indir,
                       prefix = "run_",
                       runs=c(1,180))
write.csv(alldata, file = paste0(argv$outdir, "allruns.csv"))

alldata$run <- gsub("run_", "", alldata$run)
for(method_i in unique(alldata$method)){
  for (metric in c("rmsd", "jsd")){
    subdata <- alldata[alldata$method == method_i,]
    newdata <- create_pivot_dataframe(data = subdata,
                                      value_column = metric) 
    newdata <- newdata [, order(as.numeric(names(newdata)))]
    write.csv(newdata, file = paste0(argv$outdir, "allruns_", method_i, "_", metric,".csv"))
    
  }
}


#experiments:
#t1-50 insize test

#gather indata 
insize_data <- gather_runs(indir = argv$indir,
                           prefix = "run_",
                           runs=c(1,50))
insize_annot_tags <- c(rep(7200, length(insize_data$run[insize_data$run == unique(insize_data$run)[1]])*10),
                       rep(500, length(insize_data$run[insize_data$run == unique(insize_data$run)[1]])*10),
                       rep(1000, length(insize_data$run[insize_data$run == unique(insize_data$run)[1]])*10),
                       rep(2500, length(insize_data$run[insize_data$run == unique(insize_data$run)[1]])*10),
                       rep(5000, length(insize_data$run[insize_data$run == unique(insize_data$run)[1]])*10))

insize_data$n_cells <- insize_annot_tags
insize_data <- na.omit(insize_data)

#plot scatter
#prep data
mean_tab <- plot_scatter_insize(indata = insize_data,
                                       method = "none")

# Calculate the average RMSD for each method per n_cells
rmsd_plot <- ggplot(mean_tab, aes(x = n_cells, y = rmsd, color = method)) +
  geom_line() +            # Line plot
  geom_point() +           # Add points
  xlim(0, max(mean_tab$n_cells)) +
  geom_errorbar(aes(ymin = mean_tab$rmsd - mean_tab$rmsd_se, ymax = mean_tab$rmsd + mean_tab$rmsd_se), 
                width = 0.2, 
                size = 0.8) +
  labs(title = "Multiple Groups Line Plot",
       x = "n_cells",
       y = "RMSD",
       color = "Group") + # Labels and title
  theme_minimal()
ggsave(paste0(argv$outdir, "_insize_scatter_rmsd", ".png"),
       plot = rmsd_plot, width = 6, height = 4, dpi = 300)


#jsd

jsd_plot <- ggplot(mean_tab, aes(x = n_cells, y = jsd, color = method)) +
  geom_line() +            # Line plot
  geom_point() +           # Add points
  xlim(0, max(mean_tab$n_cells)) +
  geom_errorbar(aes(ymin = mean_tab$jsd - mean_tab$jsd_se, ymax = mean_tab$jsd + mean_tab$jsd_se), 
                width = 0.2, 
                size = 0.8) +
  labs(title = "Multiple Groups Line Plot",
       x = "n_cells",
       y = "JSD",
       color = "Group") +  # Labels and title
  theme_minimal()

ggsave(paste0(argv$outdir, "_insize_scatter_jsd", ".png"),
       plot = jsd_plot, width = 6, height = 4, dpi = 300)










#t50-80 differ datasets
#gather data
datasets_data <- gather_runs(indir = argv$indir,
                           prefix = "run_",
                           runs=c(51,80))
dataset_annot_tags <- c(rep("HER2+", length(datasets_data$run[datasets_data$run == unique(datasets_data$run)[1]])*10),
                        rep("ER+", length(datasets_data$run[datasets_data$run == unique(datasets_data$run)[1]])*10),
                        rep("TNBC", length(datasets_data$run[datasets_data$run == unique(datasets_data$run)[1]])*10) )
                          
datasets_data$annot <- dataset_annot_tags
datasets_data <- na.omit(datasets_data)

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
#methods per dataset KW and dunn
for(dataset_type in unique(dataset_annot_tags)){
  datasets_data_dataset <- datasets_data[datasets_data$annot == dataset_type, ]
  dunn_methodperdataset <- capture.output(dunn.test(datasets_data_dataset$rmsd, datasets_data_dataset$method, method="bh", kw = TRUE,
                                            label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
  my_list <- dunn.test(datasets_data_dataset$rmsd, datasets_data_dataset$method, method="bh", kw = TRUE,
                       label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
  writeLines(capture.output(print(my_list)), paste0(argv$outdir, dataset_type,"_dunn_methodsperdataset_rmsd_all.txt"))
  
  write(x = dunn_methodperdataset,
            file = paste0(argv$outdir, dataset_type,"_dunn_methodsperdataset_rmsd.txt"))
}
for(dataset_type in unique(dataset_annot_tags)){
  datasets_data_dataset <- datasets_data[datasets_data$annot == dataset_type, ]
  datasets_data_methodperdataset <- capture.output(dunn.test(datasets_data_dataset$jsd, datasets_data_dataset$method, method="bh", kw = TRUE,
                                                    label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
  my_list <- dunn.test(datasets_data_dataset$jsd, datasets_data_dataset$method, method="bh", kw = TRUE,
                       label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
  writeLines(capture.output(print(my_list)), paste0(argv$outdir, dataset_type,"_dunn_methodsperdataset_jsd_all.txt"))
  write(x = dunn_methodperdataset,
        file = paste0(argv$outdir, dataset_type,"_dunn_methodsperdataset_jsd.txt"))
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

#kw and dunn tests per dataset
datasets_data_dataset <- datasets_data
dunn_perdataset <- capture.output(dunn.test(datasets_data_dataset$rmsd, datasets_data_dataset$annot, method="bh", kw = TRUE,
                                          label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
my_list <- dunn.test(datasets_data_dataset$rmsd, datasets_data_dataset$annot, method="bh", kw = TRUE,
                     label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)

writeLines(capture.output(print(my_list)), paste0(argv$outdir, "_dunn_perdataset_rmsd_all.txt"))
write(x = dunn_perdataset,
          file = paste0(argv$outdir, "_dunn_perdataset_rmsd.txt"))

dunn_perdataset <- capture.output(dunn.test(datasets_data_dataset$jsd, datasets_data_dataset$annot, method="bh", kw = TRUE,
                                            label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
my_list <- dunn.test(datasets_data_dataset$jsd, datasets_data_dataset$annot, method="bh", kw = TRUE,
                     label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
writeLines(capture.output(print(my_list)), paste0(argv$outdir, "_dunn_perdataset_jsd_all.txt"))
write(x = dunn_perdataset,
      file = paste0(argv$outdir, "_dunn_perdataset_jsd.txt"))


#plot celltypes
try(
for(metric in c("rmsd", "jsd")){
  for(method in unique(datasets_data$method)){#NOT WORKING
    print(method)
    # pdf(paste0(argv$outdir, "datasets_boxplot_percelltpyepermethod", method, metric, ".pdf"))
    plot_boxplot(indata = datasets_data,
                 plot_metric = metric,
                 annot = "TNBC",
                 annot_col = "annot",
                 annot2 = method,
                 annot2_col = "method",
                 cats = "celltype",
                 groups = "none",
                 ggplot_type = "1")
    # dev.off()
  }
}
)
#kw and dunn per celltypes
for(dataset_type in unique(datasets_data$method)){
  datasets_data_dataset <- datasets_data[datasets_data$method == dataset_type, ]
  dunn_methodperdataset <- capture.output(dunn.test(datasets_data_dataset$rmsd, datasets_data_dataset$celltype, method="bh", kw = TRUE,
                                                    label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
  my_list <- dunn.test(datasets_data_dataset$rmsd, datasets_data_dataset$celltype, method="bh", kw = TRUE,
                       label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
  writeLines(capture.output(print(my_list)), paste0(argv$outdir, dataset_type,"_dunn_percelltype_rmsd_all.txt"))
  write(x = dunn_methodperdataset,
        file = paste0(argv$outdir, dataset_type,"_dunn_percelltype_rmsd.txt"))
}
for(dataset_type in unique(datasets_data$method)){
  datasets_data_dataset <- datasets_data[datasets_data$method == dataset_type, ]
  dunn_methodperdataset <- capture.output(dunn.test(datasets_data_dataset$jsd, datasets_data_dataset$celltype, method="bh", kw = TRUE,
                                                    label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
  my_list <- dunn.test(datasets_data_dataset$jsd, datasets_data_dataset$celltype, method="bh", kw = TRUE,
                       label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
  writeLines(capture.output(print(my_list)), paste0(argv$outdir, dataset_type,"_dunn_percelltype_jsd_all.txt"))
  write(x = dunn_methodperdataset,
        file = paste0(argv$outdir, dataset_type,"_dunn_percelltype_jsd.txt"))
}




#t80-180 min density for id
#gather the data 
mintest_data <- gather_runs_mintest(indir = argv$indir,
                                    prefix = "run_",
                                    runs=c(81,180))

#gather annot
mintest_density_tags <- c(rep(c(rep(0.05, 135),
                          rep(0.1, 135),
                          rep(0.2, 135),
                          rep(0.5, 135),
                          rep(0.8, 135)), 4))
# mintest_data$density <- mintest_density_tags
##debugstart
allrunlist <- c(paste0("run_", seq(81,180)))
notindata <- setdiff(allrunlist, unique(mintest_data$run))
print(notindata)
mintest_density_tags <- mintest_density_tags[-c(2106:2132)]
mintest_data$density <- mintest_density_tags

##debug end



mintest_celltype_tags <- c(rep("B-cells", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25),
                           rep("Plasmablasts", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25),
                           rep("T-cells", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*25),
                           rep("Myeloid", length(mintest_data$run[mintest_data$run == unique(mintest_data$run)[1]])*24))
mintest_data$select_celltype <- mintest_celltype_tags
mintest_data <- na.omit(mintest_data)


#create scatter
for(method_select in unique(mintest_data$method)){
  method_select <- as.character(method_select)
  #prep data
  mean_tab <- plot_scatter_mintest(indata = mintest_data,
                           method = method_select)
  
  # pdf(paste0(argv$outdir, method_select, "_mintest_scatter_rmsd", ".pdf"))
  rmsd_plot <- ggplot(mean_tab, aes(x = density, y = rmsd, color = select_celltype)) +
    geom_line() +            # Line plot
    geom_point() +       
    xlim(0, max(mean_tab$density)) +
    geom_errorbar(aes(ymin = mean_tab$rmsd - mean_tab$rmsd_se, ymax = mean_tab$rmsd + mean_tab$rmsd_se), 
                  width = 0.05, 
                  size = 0.8) +
    labs(title = "Multiple Groups Line Plot",
         x = "Density",
         y = "RMSD",
         color = "Group") +
    theme_minimal()
    
  
  ggsave(paste0(argv$outdir, method_select, "_mintest_scatter_rmsd", ".png"),
         plot = rmsd_plot, width = 6, height = 4, dpi = 300)
  

  # dev.off()
  
  #jsd
  # pdf(paste0(argv$outdir, method_select, "_mintest_scatter_jsd", ".pdf"))
  jsd_plot <- ggplot(mean_tab, aes(x = density, y = jsd, color = select_celltype)) +
    geom_line() +            # Line plot
    geom_point() +           # Add points
    xlim(0, max(mean_tab$density)) +
    geom_errorbar(aes(ymin = mean_tab$jsd - mean_tab$jsd_se, ymax = mean_tab$jsd + mean_tab$jsd_se), 
                  width = 0.05, 
                  size = 0.8) +
    labs(title = "Multiple Groups Line Plot",
         x = "Density",
         y = "JSD",
         color = "Group") +
    theme_minimal()
  
  ggsave(paste0(argv$outdir, method_select, "_mintest_scatter_jsd", ".png"),
         plot = jsd_plot, width = 6, height = 4, dpi = 300)
  

  # dev.off()

}

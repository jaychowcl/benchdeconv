library(dunn.test)
library(argparser)

#gather all stats 
gather_all_stats <- function(indir = "./data/results/all_runs",
                             outdir = "./data/results/all_runs",
                             mintest = 0,
                             runs = c(1,10)){

  #get all file names for rmsd and jsd
  if(mintest != 0){
  rmsd_files <- list.files(path=indir,
                             pattern = "rmsd_mintest\\.csv$",
                             full.names = TRUE)
    jsd_files <- list.files(path = indir,
                            pattern = "jsd_mintest\\.csv$",
                            full.names = TRUE)
  } else {
    
    rmsd_files <- list.files(path=indir,
                             pattern = "rmsd\\.csv$",
                             full.names = TRUE)
    jsd_files <- list.files(path = indir,
                            pattern = "jsd\\.csv$",
                            full.names = TRUE)
  }
  
  #iterate through all files and create final all table
  final_tab <- data.frame()
  for(i in 1:length(rmsd_files)){
    infile_rmsd <- read.csv(paste0(rmsd_files[i]))
    infile_jsd <- read.csv(paste0(jsd_files[i]))
    
    run_no <- regmatches(rmsd_files[i], regexpr("run_(\\d+)", rmsd_files[i]))
    
    if(mintest == 0){
    new_tab <- data.frame(method = infile_rmsd$method,
                          celltype = infile_rmsd$celltype,
                          rmsd = infile_rmsd$rmsd,
                          jsd = infile_jsd$jsd_table.jsd,
                          run = run_no)
    } else {
      new_tab <- data.frame(method = infile_rmsd$method,
                            celltype = infile_rmsd$celltype,
                            rmsd = infile_rmsd$rmsd,
                            jsd = infile_jsd$jsd_table.jsd,
                            run = run_no,
                            density = mintest)
    }
    
    final_tab <- rbind(final_tab, new_tab)
    
  }
  
  #now avg of runs
  rmsd_list <- c()
  jsd_list <- c()
  method_list <- c()
  if(mintest == 0){
  for (method in unique(final_tab$method)){
    for (celltype in unique(final_tab$celltype)){
      total_rmsd <- final_tab[final_tab$method == method & final_tab$celltype == celltype, "rmsd"]
      total_jsd <- final_tab[final_tab$method == method & final_tab$celltype == celltype, "jsd"]
      
      rmsd_list <- c(rmsd_list, mean(total_rmsd))
      jsd_list <- c(jsd_list, mean(total_jsd))
      method_list <- c(method_list, method)
    }
  }
  
  summary_tab <- data.frame(method = method_list,
                            celltype = unique(final_tab$celltype),
                            rmsd_avg = rmsd_list,
                            jsd_avg = jsd_list)
  }
  
  #now method stats
  method_tab_all <- c()
  if(mintest == 0){
    for (method in unique(final_tab$method)){
      total_rmsd <- summary_tab[summary_tab$method == method, "rmsd_avg"]
      total_jsd <- summary_tab[summary_tab$method == method, "jsd_avg"]
    
      method_tab <- data.frame(method = method,
                               rmsd_avg = mean(total_rmsd),
                               jsd_avg = mean(total_jsd))
      method_tab_all <- rbind(method_tab_all, method_tab)
    
    }
  }
  
  if(mintest == 0){
  write.csv(final_tab, file = paste0(outdir, "/all_stats.csv"))
  write.csv(summary_tab, file = paste0(outdir, "/SUMMARY.csv"))
  write.csv(method_tab_all, file = paste0(outdir, "/SUMMARY_methods.csv"))
  } else {
    write.csv(final_tab, file = paste0(outdir, "/all_stats_mintest.csv"))
    write.csv(summary_tab, file = paste0(outdir, "/SUMMARY_mintest.csv"))
    write.csv(method_tab_all, file = paste0(outdir, "/SUMMARY_methods_mintest.csv"))
  }
  
  return(list(final_tab = final_tab,
              summary_tab = summary_tab,
              method_tab = method_tab_all))
}

#plot boxplot

plot_boxplots <- function(final_tab = stats_tabs$final_tab,
                               plot_type = "rmsd"){
  # Ensure the data is a data frame
  final_tab <- as.data.frame(final_tab)
  
  # Check if the plot_type is a valid column name
  if (!plot_type %in% colnames(final_tab)) {
    stop("Invalid plot_type. Please provide a valid column name.")
  }
  
  # Define unique cell types and assign symbols and colors
  cell_types <- unique(final_tab$celltype)
  symbols <- 1:length(cell_types) # Using different point symbols (1 to number of cell types)
  colors <- rainbow(length(cell_types)) # Assign different colors to each cell type
  
  # Create a mapping of cell types to symbols and colors
  celltype_symbol_map <- setNames(symbols, cell_types)
  celltype_color_map <- setNames(colors, cell_types)
  
  # Create the horizontal boxplot for the selected plot_type
  boxplot(final_tab[[plot_type]] ~ final_tab$method,
          horizontal = TRUE,
          xlab = plot_type,
          ylab = "Method",
          main = paste(toupper(plot_type), "Values by Method"),
          las = 1)
  
  # Add jittered points with different symbols and colors for each cell type
  for (i in 1:nrow(final_tab)) {
    points(final_tab[[plot_type]][i], jitter(as.numeric(factor(final_tab$method))[i]), 
           pch = celltype_symbol_map[final_tab$celltype[i]], 
           col = celltype_color_map[final_tab$celltype[i]])
  }
  
  # Add a smaller legend
  legend("topright", legend = cell_types, pch = symbols, col = colors, title = "Cell Type", cex = 0.7)
}

#Kruskal Wallis Tests
do_kruskal_wallis_dunn <- function(stats_tabs = stats_tabs,
                                   outdir = "./data/results/all_runs/SUMMARY"){
  rmsd_kruskal <- kruskal.test(stats_tabs$final_tab$rmsd ~ stats_tabs$final_tab$method)
  jsd_kruskal <- kruskal.test(stats_tabs$final_tab$jsd ~ stats_tabs$final_tab$method)
  
  # Initialize lists to store results
  method_names <- c()
  rmsd_stats <- c()
  rmsd_pvalues <- c()
  jsd_stats <- c()
  jsd_pvalues <- c()
  
  # Extract overall test results and store in lists
  method_names <- c(method_names, "Overall")
  rmsd_stats <- c(rmsd_stats, rmsd_kruskal$statistic)
  rmsd_pvalues <- c(rmsd_pvalues, rmsd_kruskal$p.value)
  jsd_stats <- c(jsd_stats, jsd_kruskal$statistic)
  jsd_pvalues <- c(jsd_pvalues, jsd_kruskal$p.value)
  
  # Initialize an empty list to store intra-method test results
  intra_methods <- list()
  
  # Loop through each unique method and perform Kruskal-Wallis test for rmsd and jsd within each method
  for (method in unique(stats_tabs$final_tab$method)) {
    method_tabs <- stats_tabs$final_tab[stats_tabs$final_tab$method == method, ]
    intra_methods[[method]] <- list(
      rmsd = kruskal.test(method_tabs$rmsd ~ method_tabs$celltype),
      jsd = kruskal.test(method_tabs$jsd ~ method_tabs$celltype)
    )
    # Append intra-method test results to the lists
    method_names <- c(method_names, method)
    rmsd_stats <- c(rmsd_stats, intra_methods[[method]]$rmsd$statistic)
    rmsd_pvalues <- c(rmsd_pvalues, intra_methods[[method]]$rmsd$p.value)
    jsd_stats <- c(jsd_stats, intra_methods[[method]]$jsd$statistic)
    jsd_pvalues <- c(jsd_pvalues, intra_methods[[method]]$jsd$p.value)
  }
  
  # Create a dataframe from the lists
  results_df <- data.frame(
    method = method_names,
    rmsd_statistic = rmsd_stats,
    rmsd_pvalue = rmsd_pvalues,
    jsd_statistic = jsd_stats,
    jsd_pvalue = jsd_pvalues
  )
  
  # Do Dunn tests
  dunn_rmsd_tab <- capture.output(dunn.test(stats_tabs$final_tab$rmsd, stats_tabs$final_tab$method, method="bh", kw = TRUE,
                         label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
  dunn_rmsd <- dunn.test(stats_tabs$final_tab$rmsd, stats_tabs$final_tab$method, method="bh", kw = TRUE,
                         label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
  dunn_jsd_tab <- capture.output(dunn.test(stats_tabs$final_tab$jsd, stats_tabs$final_tab$method, method="bh", kw = TRUE,
                                           label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
  dunn_jsd <- dunn.test(stats_tabs$final_tab$jsd, stats_tabs$final_tab$method, method="bh", kw = TRUE,
                        label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
  
  # Per intramethod dunn
  intra_methods_dunn <- list()
  for (method in unique(stats_tabs$final_tab$method)) {
    method_tabs <- stats_tabs$final_tab[stats_tabs$final_tab$method == method, ]
    
    dunn_rmsd_tab_intra <- capture.output(dunn.test(method_tabs$rmsd, method_tabs$celltype, method="bh", kw = TRUE,
                                              label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
    dunn_rmsd_intra <- dunn.test(method_tabs$rmsd, method_tabs$celltype, method="bh", kw = TRUE,
                           label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
    dunn_jsd_tab_intra <- capture.output(dunn.test(method_tabs$jsd, method_tabs$celltype, method="bh", kw = TRUE,
                                             label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05))
    dunn_jsd_intra <- dunn.test(method_tabs$jsd, method_tabs$celltype, method="bh", kw = TRUE,
                          label = TRUE, wrap = TRUE, table = TRUE, alpha = 0.05)
    
    
    intra_methods_dunn[[method]] <- list(
      dunn_rmsd_tab = dunn_rmsd_tab_intra,
      dunn_rmsd = dunn_rmsd_intra,
      dunn_jsd_tab = dunn_jsd_tab_intra,
      dunn_jsd = dunn_jsd_intra
    )
  }
  
  #write to file
  write.csv(results_df, file = paste0(outdir, "/kruskal_tests.csv"))
  
  for (method in unique(stats_tabs$final_tab$method)){
    writeLines(intra_methods_dunn[[method]]$dunn_rmsd_tab, paste0(outdir, "/rmsd_dunn_", method, ".txt"))
    writeLines(intra_methods_dunn[[method]]$dunn_jsd_tab, paste0(outdir, "/jsd_dunn_", method, ".txt"))
  }

  return(list(kruskal_wallis = results_df,
         dunn = intra_methods_dunn))
}


# Rscript ./scripts/run_stats.R --indir ./data/results/all_runs --outdir ./data/results/all_runs/SUMMARY
##get args
stats_args <- input_args <- arg_parser("run_stats.R: A tool for generating further stats after running benchdeconv through high throughput computing.")

stats_args <- add_argument(stats_args, "--indir", help="Input benchdeconv results directory.",
                           default = "./data/results/all_runs")

stats_args <- add_argument(stats_args, "--outdir", help="Output results directory",
                           default = "./data/results/all_runs/SUMMARY")

stats_args <- add_argument(stats_args, "--mintest", help="Flag for doing stats for min tests as well",
                           default = 0)

argvstats <- parse_args(stats_args)


##run
#make dirs
if (!dir.exists(argvstats$outdir)){
  dir.create(argvstats$outdir, showWarnings = TRUE, recursive = TRUE)
}

#get stats tables
stats_tabs <- gather_all_stats(indir = argvstats$indir,
                              outdir = argvstats$outdir,
                              mintest = 0)
# and _mintest tables
if(argv$mintest !=0){
stats_tabs_mintest <- gather_all_stats(indir = argvstats$indir,
                                       outdir = argvstats$outdir,
                                       mintest = 1)
}

#make rmsd boxplot
pdf(paste0(argvstats$outdir, "/boxplot_rmsd.pdf"))
plot_boxplots(final_tab = stats_tabs$final_tab,
              plot_type = "rmsd")
dev.off()
#make jsd boxplot
pdf(paste0(argvstats$outdir, "/boxplot_jsd.pdf"))
plot_boxplots(final_tab = stats_tabs$final_tab,
              plot_type = "jsd")
dev.off()


#Kruskal Wallis and dunn Tests
kwd_tests <- do_kruskal_wallis_dunn(stats_tabs = stats_tabs,
                                    outdir = argvstats$outdir)




  
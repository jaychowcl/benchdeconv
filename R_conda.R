#!/exports/eddie/scratch/s2600569/envs/benchdeconv_env/bin/R

options(repos = c(CRAN = "https://cloud.r-project.org/"))

if (!require("BiocManager", quietly = TRUE))#
	    install.packages("BiocManager")

install.packages('devtools')
library(devtools)

install_version("Matrix", version = "1.6.5", repos = "http://cran.us.r-project.org")


BiocManager::install("DropletUtils")

install.packages('SeuratObject')
install.packages('Seurat')
install.packages("igraph")
install.packages(c("BPCells", "presto", "glmGamPoi"))

BiocManager::install("DropletUtils")
BiocManager::install("SpatialExperiment")
BiocManager::install("SingleCellExperiment")


install.packages("pak")
pak::pkg_install("omnideconv/spacedeconv", dependencies=TRUE)

#install.packages("spacexr") # RCTD
#devtools::install_github('YMa-lab/CARD') # CARD
#install.packages("SPOTlight") # SPOTlight

BiocManager::install("scran") # spatialdwls
BiocManager::install("BiocParallel")# spatialdwls
BiocManager::install("genefilter")# spatialdwls
BiocManager::install("DESeq2")# spatialdwls
devtools::install_github('edsgard/trendsceek')# spatialdwls
install.packages('multinet')# spatialdwls
install.packages('FactoMineR')# spatialdwls

BiocManager::install("SpatialDecon")# spatialdecon

install.packages("philentropy")#jsd stats 


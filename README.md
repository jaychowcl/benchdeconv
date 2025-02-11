[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![Picture2](https://github.com/user-attachments/assets/2f5e7a57-e7d8-4828-aaa9-5e9ec0ac177c)

A benchmarking pipeline for spatial deconvolution methods. 

## Description

Spatial deconvolution is the prediction of cell type identities and fractions within each spot of spatially resolved transcriptomic datasets, either involving direct statistical inference or machine learning 
to produce a signature matrix or to directy predict cell fractions. Here we provide a benchmarking pipeline to test the performance of different spatial deconvolution methods with a focus on tumour specific metrics.

Please see Understanding_the_Immune_Tumour_Microenvironment_by_Integrating_Single-Cell_and_Spatial_Transcriptomics.pdf for full dissertation. 

### General Overview
1. Pipeline takes annotated single cell RNAseq datasets and each splits sample into two: a)synthetic scRNA dataset (testing dataset) and b)training dataset with each containing the same number of cell types.
2. The synthic scRNA dataset is used to produce synthetic spatial transcriptomic datasets with known ground truth fractions using Synthspot [(Sangaram et al., 2024)](https://pubmed.ncbi.nlm.nih.gov/38787371/) to produce the testing dataset.
3. The training dataset is used to train relevent deconvolution methods.
4. The synthetic spatial transcriptomic datasets are then deconvoluting via the different methods. (RCTD, SPOTlight, CARD are used by default. See below for adding new methods.)
5. Predicted cell fractions are compared with ground truth fractions via various metrics.


## Getting Started

1. Install using steps below
   
(Optional): Modify 10X genomics spot coordinates in ./data/spot_coords/out1.csv to set spatial transciptomic plate size

(Optional): Modify _mintest_ spots to set spot locations of known cell fractions for testing minimum cell density for identification in ./data/spot_coords/out1_mintest.csv

3. Gather annotated single cell RNAseq data. Requires count matrix barcodes, genes, counts in a sparse matrix, and metadata of cell annotations.
4. Run run.R (see --help for flags) for importing data, generating synthetic datasets, deconvolution, and generating preliminary statistics.
5. Run run_stats.R (see --help for flags) to generate statistics across experiments.

### For parallelization
1. Use run_hpc.sh to run in a high performance computing enviornment for parallelization. Current settings follow disseratation experiments.
2. run run_stats.R (see --help for flags) to collect all HPC runs and generate statistics.

Results and statsitcs will be reported in data/results unless specified by --outdir.

### Adding new deconvolution methods
1. Place deconvolution steps in the run.R ##DECONOVLUTION section
2. Add deconovlution name into method_names variable in same section
3. Add deconvolution result output dataframe variable into method_list variable in same section

All proceding steps will now take into account the new deconvolution method. 


### Dependencies

See .yaml and R_installation.R for packages.

### Installing

1. Clone repository into workspace
2. Use benchdeconv.yaml to install dependencies into a conda environment
3. Run R_installation.R to install additional R packages into the conda environment
4. Installation done!


## Authors

[Jay Chow (Chi Lung)](https://github.com/jaychowcl)

[Dr. Florent Petitprez (Supervisor)](https://edwebprofiles.ed.ac.uk/profile/florent-petitprez)

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License


## References

* [Synthspot](https://github.com/saeyslab/synthspot)
* [Spacedeconv](https://github.com/omnideconv/spacedeconv)
* [RCTD](https://github.com/dmcable/spacexr)
* [CARD](https://github.com/YMa-lab/CARD)
* [SPOTlight](https://github.com/MarcElosua/SPOTlight)

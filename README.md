# Benchdeconv

A benchmarking pipeline for spatial deconvolution methods. 

## Description

Spatial deconvolution is the prediction of cell type identities and fractions within each spot of spatially resolved transcriptomic datasets, either involving direct statistical inference or machine learning 
to produce a signature matrix or to directy predict cell fractions. Here we provide a benchmarking pipeline to test the performance of different spatial deconvolution methods with a focus on tumour specific metrics.

Please see Understanding_the_Immune_Tumour_Microenvironment_by_Integrating_Single-Cell_and_Spatial_Transcriptomics.pdf for full dissertation. 

General Overview:
1. Pipeline takes annotated single cell RNAseq datasets and each sample into two: a)synthetic scRNA dataset and b)training dataset with each containing the same number of cell types.
2. The synthic scRNA dataset is used to produce synthetic spatial transcriptomic datasets with known ground truth fractions using Synthspot [(Sangaram et al., 2024)](https://pubmed.ncbi.nlm.nih.gov/38787371/).
3. The training dataset is used to train relevent deconvolution methods.
4. Predicted cell fractions are compared with ground truth fractions via various metrics.


## Getting Started

1. Install using steps below
(Optional): Modify 10X genomics spot coordinates in ./data/spot_coords/out1.csv to set spatial transciptomic plate size
(Optional): Modify _mintest_ spots to set spot locations of known cell fractions for testing minimum cell density for identification in ./data/spot_coords/out1_mintest.csv
2. Gather annotated single cell RNAseq data. Requires count matrix barcodes, genes, counts in a sparse matrix, and metadata of cell annotations.
3. Run run.R (see --help for flags) for importing data, generating synthetic datasets, deconvolution, and generating preliminary statistics.
(Optional): Run run_hpc.sh to run in a high performance computing enviornment for parallelization. Current settigns 
 


### Dependencies

* Describe any prerequisites, libraries, OS version, etc., needed before installing program.
* ex. Windows 10

### Installing

1. Clone repository into workspace
2. Use benchdeconv.yaml to install dependencies into a conda environment
3. Run R_installation.R to install additional R packages into the conda environment
4. Installation done!


### Executing program

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

ex. Dominique Pizzie  
ex. [@DomPizzie](https://twitter.com/dompizzie)

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)

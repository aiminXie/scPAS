# scPAS : Single-Cell Phenotype-Associated Subpopulation identifier
 A tool for identifying phenotype-associated cell subpopulations from single-cell sequencing data by integrating bulk data
### Introduction ###
`scPAS` is a new tool which enables the quantitative estimation of the strength of association between each cell in scRNA-seq data and a phenotype by constructing a network-regularized sparse regression model. This model integrates bulk RNA-seq data with phenotype information and the gene-gene similarity network from single-cell data. 

The workflow of scPAS is shown in the following Figure:

<p align="center">
<img src=Flow_diagram.png height="492" width="700">
</p>

### Installation ###
* System requirements: scPAS is developed under R (version >= 4.0.5).
* scPAS R package can be installed from Github using devtools: `devtools::install_github("aiminXie/scPAS")`.

#### Installation of other dependencies if not automatically installed
- Install [Seurat](https://github.com/satijalab/seurat) using `install.packages('Seurat')`.
- Install [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) using `install.packages('Rcpp')`.
- Install [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html) using `install.packages('Matrix')`.
- Install [preprocessCore](https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html) using `if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install( "preprocessCore")`.

### Tutorial ###
[scPAS Tutorial](https://github.com/aiminXie/scPAS/scPAS_Tutorial.html)

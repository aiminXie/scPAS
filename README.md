# scPAS ：Single-Cell Phenotype-Associated Subpopulation identifier
 A tool for identifying Phenotype-Associated cell Subpopulations from single-cell sequencing data by integrating bulk data
### Introduction ###
`scPAS` is a new tool which enables the quantitative estimation of the strength of association between each cell in scRNA-seq data and a phenotype by constructing a network-regularized sparse regression model. This model integrates bulk RNA-seq data with phenotype information and the gene-gene similarity network from single-cell data. 

The workflow of scPAS is shown in the following Figure:

<p align="center">
<img src=Flow_diagram.png height="500" width="500">
</p>

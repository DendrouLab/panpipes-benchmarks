# panpipes-benchmarks

This repo contains the code used to run panpipes benchmarks. 
To allow a fair comparison, equivalent methods where run, namely scvi, harmony and PCA for no batch correction.

We also ran a panpipes only benchmark for time and resources with scalability of number of cells in the following datasets:
1) Seurat Multimodal PBMC (RNA.ADT) post-QC at 10K, 50K , 100K, full dataset ([Seurat v4 paper](https://doi.org/10.1016/j.cell.2021.04.048))
2) TARUS dataset (RNA only) at 500K and 1M ([ biorxiv preprint](https://doi.org/10.1101/2023.05.05.539635))

The R script was submitted using the provided bash submission script and elapsed times are printed to a log file. The time is then calculated as the difference in seconds between the start and the end of the computation. 


Panpipes integration was run calling
`panpipes integration make full` and counting the run time until the process `collate_integration_outputs` which is comparable with the end of the provided R script as it's the step of the pipeline that collates the output of all the integration methods run in parallel.

#!/bin/bash

#SBATCH --mem=90GB
#SBATCH -c 6
#SBATCH --job-name=seurat_bm
#SBATCH --gres=gpu:1

echo "********************************************************"
echo "Job ID.             : "$SLURM_JOB_ID          
echo "TASK ID.            : "$SLURM_TASK_ID          
echo "Run on host         : "`hostname`
echo "Operating system    : "`uname -s`
echo "Username            : "`whoami`
echo "Started at          : "`date`
echo "********************************************************"

echo "activate environment"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate panpipes_env
echo 

echo "Starting reticulate benchmark " `date`
Rscript check_reticulate_times.R > reticulate_bm.log 2>&1

echo "Starting seurat benchmark " `date`
Rscript seurat_multimodal_v5.R > seurat_bm.log 2>&1
echo "********************************************************"
echo "Finished at          : "`date`
echo "********************************************************"

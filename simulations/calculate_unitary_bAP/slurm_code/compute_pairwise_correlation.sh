#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --account=def-akhadra
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2023a
matlab -nodisplay -r "correlation_uAP_pairwise_correlations"
# srun matlab -nodisplay -r "compute_unitaryAP_beluga"
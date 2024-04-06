#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --account=def-akhadra
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

module load matlab/2023a
matlab -nodisplay -r "compute_unitaryAP_beluga"
# srun matlab -nodisplay -r "compute_unitaryAP_beluga"
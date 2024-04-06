#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-akhadra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=synthetic_spikes.log

module load matlab/2020a
srun matlab -nodisplay -r "compute_cortical_area"
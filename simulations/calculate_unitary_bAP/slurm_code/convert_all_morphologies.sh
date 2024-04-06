#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

module load python/3.8.10
module load mpi4py
module load scipy-stack
module load neuron

folder=/lustre04/scratch/nbrake/data/simulations/unitary_AP

cd $folder
for d in $folder/*; do
    cd $d
    python /lustre04/scratch/nbrake/code/simulate_blue_brain/generate_morphology.py $d
    cd ..
done
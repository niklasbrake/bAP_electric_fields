#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
## SBATCH --output=simulation.log

module load python/3.8.10
module load mpi4py
module load scipy-stack
module load neuron
pip install /home/nbrake/pkgs/LFPykit-0.5.1-py3-none-any.whl

# Only run on models that have changed EI ratios
folder=/lustre04/scratch/nbrake/data/simulations/unitary_AP
cd $folder
file=/lustre04/scratch/nbrake/code/simulate_blue_brain/changedEI.txt

lineNumber=${SLURM_ARRAY_TASK_ID}
((lineNumber--))
while IFS="" read -r d || [ -n "$d" ]
do
  if [ $((lineNumber % 10)) -eq 0 ]; then
      cd "$d"
      python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py $d
      cd ..
  fi
  ((lineNumber++))
done < $file
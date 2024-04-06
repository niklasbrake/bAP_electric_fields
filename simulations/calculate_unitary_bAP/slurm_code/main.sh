#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=simulation-%x.%j.log

module load python/3.8.10
module load mpi4py
module load scipy-stack
module load neuron
pip install /home/nbrake/pkgs/LFPykit-0.5.1-py3-none-any.whl

folder=/lustre04/scratch/nbrake/data/simulations/unitary_AP

# cd $folder
# for d in $folder/*; do
#     cd "$d"
#     nrnivmodl ./mechanisms
#     for i in $(seq 1 10);
#     do
#         python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py $d $i False
#         # python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py $d $i True
#     done
#     cd ..
# done


# cd $folder
# for d in $folder/*; do
#     cd "$d"
#     # nrnivmodl ./mechanisms
#     python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py $d
#     cd ..
# done

cd $folder
if [ ${SLURM_ARRAY_TASK_ID} == 1 ]
then
  x=$folder/L1*
elif [ ${SLURM_ARRAY_TASK_ID} == 2 ]
then
  x=$folder/L23*
elif [ ${SLURM_ARRAY_TASK_ID} == 3 ]
then
  x=$folder/L4*
elif [ ${SLURM_ARRAY_TASK_ID} == 4 ]
then
  x=$folder/L5*
elif [ ${SLURM_ARRAY_TASK_ID} == 5 ]
then
  x=$folder/L6*
fi

for d in $x; do
  simulations=$d/matlab_recordings/*mat
  toSimulate=true
  for file in $simulations; do
    time=`date -r $file +%s`
    if [ "$time" -gt 1711497600 ]
    then
      toSimulate=false
    fi
  done
  if [ "$toSimulate" == true ]
  then
    cd "$d"
    python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py $d
    cd ..
  fi
done

# Only run on models that have changed EI ratios
# cd $folder
# file=/lustre04/scratch/nbrake/code/simulate_blue_brain/changedEI.txt
# while IFS="" read -r d || [ -n "$d" ]
# do
#   cd "$d"
#   python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py $d
#   cd ..
# done < $file


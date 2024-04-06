cd /lustre04/scratch/nbrake/data/simulations/blue_brain/sampled_mtypes
for /D %%i in (*) do (
    cd %%i
    nrnivmodl ./mechanisms
    for /l %%x in (1, 1, 10) do (
        python /lustre04/scratch/nbrake/code/simulate_blue_brain/run_LFPy.py %%i %%x
    )
    cd ..
)
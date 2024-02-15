L1_DAC_bNAC219_3cd E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\other_code\representative_models
E:
for /D %%i in (*) do (
    cd %%i
    for /l %%x in (1, 1, 10) do (
        python C:\Users\brake\Documents\GitHub\spiking_amplitude_estimation\simulations\blue_neurons\run_LFPy.py %%i %%x
    )
    cd ..
)
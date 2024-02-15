cd E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\other_code\neuron_models\
E:
for /D %i in (*) do (
    cd %i
    nrnivmodl mechanisms
    python C:\Users\brake\Documents\GitHub\spiking_amplitude_estimation\simulations\blue_neurons\generate_morphology.py %i
    cd ..
)
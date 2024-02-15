SET file0=nbrake@beluga.computecanada.ca:/lustre04/scratch/nbrake/data/simulations/unitary_AP
SET file1=E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models

rem download all files from folder
rem cd %file1%
rem E:
rem for /D %i in (*) do (
rem     rem scp -r %file0%/%i/matlab_recordings/ %file1%\%i
rem     scp %file0%/%i/morphData.mat %file1%/%i/
rem )

rem only execute on modified simulations
SET file=E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\changedEI.txt
for /F "tokens=*" %i in (%file%) do (
    scp -r %file0%/%i/matlab_recordings/ %file1%\%i
)
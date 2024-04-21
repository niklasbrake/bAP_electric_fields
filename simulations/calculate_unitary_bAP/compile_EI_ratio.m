folder = '/lustre04/scratch/nbrake/data/simulations/unitary_AP_EI_ratio';
F = dir(folder); F = F(3:end);

EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};

data = struct();
for i = 1:length(F)
    folder0 = fullfile(folder,F(i).name,'matlab_recordings');
    data.(F(i).name) = struct();
    data.(F(i).name).passive = struct('dipoles',zeros(160001,3,10),'voltage',zeros(160001,10),'time',zeros(160001,10));
    data.(F(i).name).active = struct('dipoles',zeros(160001,3,10),'voltage',zeros(160001,10),'time',zeros(160001,10));

    for j = 1:length(EI_vec)
        passiveFile = ['synaptic_input_EI' EI_vec{j} '_passive.mat'];
        load(fullfile(folder0,passiveFile));
        data.(F(i).name).passive.dipoles(:,:,j) = dipoles;
        data.(F(i).name).passive.voltage(:,j) = voltage;
        data.(F(i).name).passive.time = time;
        
        activeFile = ['synaptic_input_EI' EI_vec{j} '.mat'];
        load(fullfile(folder0,activeFile));
        data.(F(i).name).active.dipoles(:,:,j) = dipoles;
        data.(F(i).name).active.voltage(:,j) = voltage;
        data.(F(i).name).active.time = time;
    end
end

save('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/EI_ratio.mat','data');
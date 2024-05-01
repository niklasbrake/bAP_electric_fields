SSExy_D = zeros(8e3,1e3);
Pxy_D = zeros(8e3,1e3);
count_D = zeros(1,1e3);
for i = 1:20
       data = load(sprintf('/lustre04/scratch/nbrake/data/simulation_analyzed/cross_spectra/chain%d.mat',i))
       SSExy_D = SSExy_D+data.SSExy_D;
       Pxy_D = Pxy_D+data.Pxy_D;
       count_D = count_D+data.count_D;
end

clearvars data i
save('/lustre04/scratch/nbrake/data/simulation_analyzed/cross_spectra/all_chains.mat')
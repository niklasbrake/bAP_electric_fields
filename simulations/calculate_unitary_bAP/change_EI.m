folder = '/lustre04/scratch/nbrake/data/simulations/unitary_AP';

load('/lustre04/scratch/nbrake/code/simulate_blue_brain/newEI.mat')

for i = 1:length(idcs)
    j = idcs(i);
    fprintf('%s, %.2f, EI: %f -> %f\n',files{j},N(j),EI(j),newEI)
    csvwrite(fullfile(folder,files{j},'EI_ratio.csv'),newEI);
end

fid = fopen('/lustre04/scratch/nbrake/code/simulate_blue_brain/changedEI.txt','w');
fprintf(fid,'%s\n',files{idcs});
fclose(fid)
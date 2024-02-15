folder1 = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\other_code\neuron_models';
folder2 = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\other_code\representative_models';

x = dir(folder1)
x = x(3:end)
names = {x(:).name};
mTypes = unique(cellfun(@(x) x(1:end-2),names,'UniformOutput',false))

for i = 1:length(mTypes)
    idcs = find(cellfun(@(x) ~isempty(strfind(x,mTypes{i})),names));
    idx = idcs(randsample(length(idcs),1));

    copyfile(fullfile(folder1,names{idx}),fullfile(folder2,names{idx}))
end
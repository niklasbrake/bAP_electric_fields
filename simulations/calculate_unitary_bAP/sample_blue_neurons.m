function bs_sample = sample_blue_neurons(N)
    load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\unitary_AP_PSD.mat','mtype');
    load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\mtype_abundance.mat','mtype_abundance');
    [J,ID] = findgroups(mtype);
    for i = 1:length(ID)
        abundance(i) = mtype_abundance(ID{i},:).Abundance;
    end
    mtype_abundance.count = splitapply(@(x) sum(x),ones(size(J)),J)';
    scaledAbundance = mtype_abundance.Abundance./mtype_abundance.count;
    prop = scaledAbundance(J);
    F = cumsum(prop)/sum(prop);
    [~,idcs] = unique(F);
    R = rand(N,1);
    bs_sample = interp1(F(idcs),idcs,R,'next','extrap');
end
folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models';
F = dir(folder); F = F(3:end);

N = zeros(length(F),3);
asym_idx = zeros(length(F),3);
for i = 1:length(F)
    file = fullfile(folder,F(i).name,'morphData.mat');
    [asym_idx(i,:),N(i,:)] = main(file);
end

function [asym_idx,N] = main(file)
    load(file)
    connections = connections+1;


    xSoma = data(data(:,2)==1,3:6);
    temp = data(:,5);
    data(:,3:5) = data(:,3:5) - mean(xSoma(:,1:3));

    count = 0;
    for i = 1:length(segs)
        idcs = segs{i};
        for j = 2:length(idcs)
            count = count+1;
        end
    end

    SA = zeros(count,1);
    x0 = zeros(count,3);
    count = 1;
    for i = 1:length(segs)
        idcs = segs{i};
        xSegment = data(idcs,3:6);
        for j = 2:length(idcs)
            h = norm(xSegment(j,1:3)-xSegment(j-1,1:3));
            r1 = xSegment(j,4);
            r2 = xSegment(j-1,4);
            SA(count) = 1/3*pi*h*(r1^2+r1*r2+r2^2);
            x0(count,:) = 0.5*(xSegment(j,1:3)+xSegment(j-1,1:3));
            count=count+1;
        end
        % asym_idx(i,:) = sum(x0.*SA);
    end
    asym_idx = sum(x0.*SA)./std(x0);
    N = std(x0);
end
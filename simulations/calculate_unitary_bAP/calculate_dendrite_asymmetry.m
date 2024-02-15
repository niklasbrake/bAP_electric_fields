folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models';
F = dir(folder); F = F(3:end);

N = zeros(length(F),3);
asym_idx = zeros(length(F),3);
for i = 1:length(F)
    file = fullfile(folder,F(i).name,'morphData.mat');
    [asym_idx(i,:),N(i,:)] = main(file);
end


% V = max(abs(squeeze(savedUnitaryAP(:,2,:))));
% X = [N(:),abs(asym_idx(:,2)),V(:)];
% I = (10*(X(:,2)-80)+(X(:,1)-1e3) < X(:,3));

function [asym_idx,N] = main(file)
    load(file)
    connections = connections+1;


    xSoma = data(data(:,2)==1,3:6);
    temp = data(:,5);
    data(:,3:5) = data(:,3:5) - mean(xSoma(:,1:3));
    N = data(:,6);
    asym_idx = sum(data(:,3:5).*N);%./sum(N);
    N = sum(N);
    % return;
    span = std(data(:,3:5));

    asym_idx = zeros(length(segs),3);
    for i = 1:length(segs)
        idcs = segs{i};
        xSegment = data(idcs,3:6);
        SA = zeros(length(idcs)-1,1);
        x0 = zeros(length(idcs)-1,3);
        for j = 2:length(idcs)
            h = norm(xSegment(j,1:3)-xSegment(j-1,1:3));
            r1 = xSegment(j,4);
            r2 = xSegment(j-1,4);
            % SA(j-1) = pi*(r1+r2)*sqrt(h^2+(r1-r2)^2);
            SA(j-1) = 1/3*pi*h*(r1^2+r1*r2+r2^2);
            x0(j-1,:) = 0.5*(xSegment(j,1:3)+xSegment(j-1,1:3));
        end
        asym_idx(i,:) = sum(x0.*SA./span);
    end
    asym_idx = sum(asym_idx);
    N = span;
end
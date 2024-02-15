folder = 'E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\resources\morphology_proportions';
F = dir(folder);
F = F(3:end);
nc = [6,1,9,3,9,4,9,5,9];

for i = 5

    I0 = imread(fullfile(folder,F(i).name));
    I = double(reshape(I0,[],3))/255;

    idcs = find(and(vecnorm(I-mean(I,2),2,2)>0.1,vecnorm(I,2,2)>0.1));
    [i2,j,jc] = unique(round(I(idcs,:),2),'rows');
    [K,xi] = mvksdensity(I(idcs,:),i2,'BandWidth',0.05);
    idcs2 = find(K>5);
    [g,C] = kmeans(i2(idcs2,:),nc(i));

    subCount = zeros(length(idcs2),1);
    for l = 1:length(idcs2)
        subCount(l) = sum(ismember(jc,idcs2(l)));
    end
    fig = figureNB(14,7);
    subplot(1,2,1);
        imshow(I0)
    % subplot(2,2,3);
        % S = scatter3(I(:,1),I(:,2),I(:,3),5,[0.5,0.5,0.5],'filled');
        % S.MarkerFaceAlpha = 0.1;
        % hold on;
        % scatter3(i2(idcs2,1),i2(idcs2,2),i2(idcs2,3),10,I(idcs(j(idcs2)),:),'filled')

    subplot(1,2,2)
        count = splitapply(@sum,subCount(:),g(:));
        P = pie(count);
        for l = 1:2:2*nc(i)
            P(l).FaceColor = C((l-1)/2+1,:);
        end

    drawnow;
    print(fig,fullfile('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\manuscript\figures\workspace\figure1\morph_prop_estimations',F(i).name),'-dpng','-r600')
    % close(fig);
end
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cortical_area\pairwise_distance.mat')
figureNB;
plotwitherror(rValues,total_area,'SE','color','k','LineWidth',1)
xlim([0,200])
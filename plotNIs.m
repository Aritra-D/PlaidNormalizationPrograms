% Work on NIs and modeling

load('C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\all_N20_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_gridType_Microelectrode_UnitID0.mat');

xAll{1} = squeeze(firingRateData.analysisDataST(:,1,:,:));
for i=1:4
    xAll{i+1} = squeeze(energyData.analysisDataST{i}); 
    xAll{i+1} = 10.^xAll{i+1};
end

edges = [0:0.1:1.5 inf];
edgeCenters = [0.1:0.1:1.5 1.6];

for i=1:5
    x = xAll{i};
    NI = squeeze(x(:,1,5) ./ (x(:,1,1) + x(:,5,5)));
    n = histcounts(NI,edges);
    
    subplot(2,3,i);
    bar(edgeCenters,n);
    title(num2str([i median(NI)]));
end
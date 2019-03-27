% spikeData
close all; clear;

load('E:\spikeDataST.mat')
cValsUnique = [0 6.25 12.5 25 50];
colors = jet(5);

for iElec = 1: size(spikeDataST,1)
   clear spikeRateData_elec
   spikeRateData_elec = squeeze(spikeDataST(iElec,:,:)); 
   plot(cValsUnique,spikeRateData_elec(end,:), 'Marker','o','LineWidth',2,'color',colors(end,:,:))
   hold on;
   plot(cValsUnique,diag(flipud(spikeRateData_elec)),'Marker','o','LineWidth',2,'color','k')
   title(['ElecNum: ' num2str(iElec)])
   pause
   hold off;
end
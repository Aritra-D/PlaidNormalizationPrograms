function combineData_OriTuningProtocol_PlaidProtocol(monkeyName,folderSourceString,timeRangeForComputation)
if ~exist('folderSourceString','var') 
   folderSourceString = 'M:\Data\PlaidNorm\';
end

% Variable Parameters
spikeCutoff = 20;
snrCutoff = 2;
% timeRangeForComputation = [0.25 0.5]; % expressed in second
% timeRangeForComputationBL = -0.05+[-diff(timeRangeForComputation) 0];
dRange = [0 0.75];
tapers_MT = [1 1]; % parameters for MT analysis
removeERPFlag = 1;
% normalizateSpikeDataFlag = 0;

% Fixed parameters
folderSourceString_Project = strtok(folderSourceString,'\');
folderSave = fullfile(folderSourceString_Project,'Projects\PlaidNormalizationProject\savedData_Figures');
if ~exist(folderSave,'dir')
    mkdir(folderSave)
end
gridType = 'Microelectrode';
getSpikeElectrodesFlag = 1;
combineUniqueElectrodeData = 0;
unitID = 0;

% if removeERPFlag ==0
%     LFPdataProcessingMethod = 'Evoked Response';
% elseif removeERPFlag ==1
%     LFPdataProcessingMethod = 'Induced Response';
% end

% oriSelectiveFlag = 0;
fileSave1 = fullfile(folderSave,[monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) '.mat']); 

if exist(fileSave1,'file')
    disp(['Loading file ' fileSave1]);
    
    load(fileSave1);
else
   error('No data file found') 
end

for iElec = 1:size(firingRateData.analysisDataST,1)
firingRateST_Ori1(iElec,:) = squeeze(firingRateData.analysisDataST(iElec,:,end,:))';
firingRateST_Ori2(iElec,:) = flip(squeeze(firingRateData.analysisDataST(iElec,:,:,1)))'; %#ok<*AGROW>
end

firingRateST_Ori1 = cat(2,firingRateST_Ori1,oriTuningData.OriPairFR(:,1));
firingRateST_Ori2 = cat(2,firingRateST_Ori2,oriTuningData.OriPairFR(:,2));

cValsUnique = [0 6.25 12.5 25 50 100];
colors = jet(5);

for iElec = 1: size(firingRateST_Ori1,1)
   plot(cValsUnique,firingRateST_Ori1(iElec,:),'Marker','o','LineWidth',2,'color',colors(end,:,:))
   hold on;
   plot(cValsUnique,firingRateST_Ori2(iElec,:),'Marker','o','LineWidth',2,'color',colors(1,:,:))
   title(['ElecNum: ' num2str(iElec)])
   pause
   hold off;
end


end
function plotFigures(monkeyName,folderSourceString,timeRangeForComputation)
if ~exist('folderSourceString','var') 
   folderSourceString = 'M:\Data\PlaidNorm\';
end
close all; % closes any open figure to avoid any overlaying issues

% Variable Parameters
spikeCutoff = 20;
snrCutoff = 2;
% timeRangeForComputation = [0.25 0.5]; % expressed in second
timeRangeForComputationBL = -0.05+[-diff(timeRangeForComputation) 0];
dRange = [0 0.75];
tapers_MT = [1 1]; % parameters for MT analysis
removeERPFlag = 1;
normalizateSpikeDataFlag = 0;

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
% contrastIndexList{1} = [1 1]; % Plaid (Ori 1 contrast: 0% and Ori 2 contrast: 50%)
% contrastIndexList{2} = [5 5]; % Plaid (Ori 1 contrast: 50% and Ori 2 contrast: 0%)
freqRanges{1} = [8 12]; % alpha
freqRanges{2} = [30 80]; % gamma
freqRanges{3} = [104 250]; % hi-gamma
freqRanges{4} = [16 16];  % SSVEP

timeRangeParameters.blRange = timeRangeForComputationBL;
timeRangeParameters.stRange = timeRangeForComputation;
timeRangeParameters.erpRange = [0.05 0.2];

if removeERPFlag ==0
    LFPdataProcessingMethod = 'Evoked Response';
elseif removeERPFlag ==1
    LFPdataProcessingMethod = 'Induced Response';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% display properties %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 1 (Spike data from Orientation Selective Elecs-
% Pref (x-axis) Null (y-axis) axes); spike data along
% increasing contrasts of Pref Ori with contrasts of  null Ori presented 
% by different colors; absolute color map 5x5; relative color map 5x5; 
% CRF (no Norm, Pref,Pref+Null, Avg, Null) 
hFigure1 = figure(1);
set(hFigure1,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig1.hPlot1 = getPlotHandles(1,5,[0.15 0.65 0.7 0.2],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot1);
hPlotsFig1.hPlot2 = getPlotHandles(1,3,[0.15 0.23 0.7 0.2578],0.1188,0.05,0); 

% Figure 2 (Spike data from all Elecs)
hFigure2 = figure(2);
set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig2.hPlot1 = getPlotHandles(1,5,[0.15 0.65 0.7 0.2],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot1);
hPlotsFig2.hPlot2 = getPlotHandles(1,3,[0.15 0.23 0.7 0.2578],0.1188,0.05,0);
% 
% % Figure 3 (PSDs (A); deltaPSDs (B) along increasing contrast of Ori 1 
% % with contrasts of Ori 2 presented in different colors;
% % absolute color map 5x5; relative color map 5x5; 
% % CRF (no Norm, Pref,Pref+Null, Avg, Null) 
% % for alpha (C), gamma (D) and high-gamma (E) for ERP-subtracted data
hFigure3 = figure(3);
set(hFigure3,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig3.hPlot1 = getPlotHandles(2,5,[0.25 0.69 0.5 0.28],0.01,0.01,1); 
linkaxes(hPlotsFig3.hPlot1(1,:)); linkaxes(hPlotsFig3.hPlot1(2,:));
hPlotsFig3.hPlot2 = getPlotHandles(3,3,[0.25 0.08 0.5 0.51],0.15,0.04,1);
% 
% % Figure 4 (PSDs (A); deltaPSDs (B) along increasing contrast of Ori 1 
% % with contrasts of Ori 2 presented in different colors;
% % absolute color map 5x5; relative color map 5x5; 
% % CRF (no Norm, Pref,Pref+Null, Avg, Null) 
% % for SSVEP (C) for non-ERP subtracted data
hFigure4 = figure(4);
set(hFigure4,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig4.hPlot1 = getPlotHandles(2,5,[0.15 0.5 0.7 0.4],0.01,0.01,1); 
linkaxes(hPlotsFig4.hPlot1(1,:));linkaxes(hPlotsFig4.hPlot1(2,:));
hPlotsFig4.hPlot2 = getPlotHandles(1,3,[0.15 0.1 0.7 0.2578],0.1188,0.05,0);
% 
% 
% % Figure 5: Comparison of fitted parameters for FR, alpha, gamma, high
% % gamma, SSVEP--- Three plots for each NI, sigma and alpha fitted
% % parameters
hFigure5 = figure(5);
set(hFigure5,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig5.hPlot1 = getPlotHandles(6,4,[0.2 0.05 0.6 0.9],0.05,0.04,1); 

%%%%%%%%%%%%%%%%%%%%% Get Session Details for Monkey(s) %%%%%%%%%%%%%%%%%%%
fileNameStringListAll = getFileNameStringList(monkeyName,gridType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES 1 & 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Get data for all Good Electrodes %%%%%%%%%%%%%%%%%%%%
disp('all Electrodes:')
oriSelectiveFlag = 0;
electrodeList_All = getElectrodesList(fileNameStringListAll,oriSelectiveFlag,folderSourceString);

fileSave1 = fullfile(folderSave,[monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) '.mat']); 

if exist(fileSave1,'file')
    disp(['Loading file ' fileSave1]);
    load(fileSave1);
    NI_Data_allElecsInduced = NI_Data; %#ok<NODEF>
else
    % get Data all Session for monkey(s) for all Electrodes
    [erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,~] = ...
    getData(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); %#ok<ASGLU>
    save(fileSave1,'erpData','firingRateData','fftData','energyData','oriTuningData','NI_Data')
end

% Put plot Functions for figures 1,3
plotData_spikes(hPlotsFig1,firingRateData,normalizateSpikeDataFlag) % spikes for static gratings, Fig 1
rescaleData(hPlotsFig1.hPlot1,-0.1,0.5,getYLims(hPlotsFig1.hPlot1),14);
rescaleData(hPlotsFig1.hPlot2(3),0,50,getYLims(hPlotsFig1.hPlot2(3)),14);
folderSave_Figs = fullfile(folderSourceString_Project,'Projects\PlaidNormalizationProject\Figures');
if ~exist(folderSave_Figs,'dir')
    mkdir(folderSave_Figs)
end
FigName1 = fullfile(folderSave_Figs,['Figure 1_' monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID)]);
saveas(hFigure1,[FigName1 '.fig'])
saveas(hFigure1,[FigName1,'.tif'])
    

plotData_energy(hPlotsFig3,energyData) % alpha, gamma, hi-gamma for static gratings, Fig 3; 
rescaleData(hPlotsFig3.hPlot1,0,250,getYLims(hPlotsFig3.hPlot1),12);
rescaleData(hPlotsFig3.hPlot1(1,:),0,250,[-1.5 3.5],12);
rescaleData(hPlotsFig3.hPlot1(2,:),0,250,[-4 10],12);
set(findall(hFigure3,'-property','FontSize'),'FontSize',13);
FigName3 = fullfile(folderSave_Figs,['Figure 3_' monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) ]);
saveas(hFigure3,[FigName3 '.fig'])
saveas(hFigure3,[FigName3,'.tif'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Get data for Orientation selective Good Electrodes %%%%%%%%%%%
disp('Orientation-Tuned Electrodes:')
oriSelectiveFlag = 1; % Fig 2 spike data for orientation selective elecs
electrodeList_OriTuned = getElectrodesList(fileNameStringListAll,oriSelectiveFlag,folderSourceString);

clear erpData firingRateData fftData energyData NI_Data % clear data for all Elecs

% Orientation-tuned electrode data (only required for Figure 2)
fileSave2 = fullfile(folderSave,[monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_oriTunedElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) '.mat']); 

if exist(fileSave2,'file')
    disp(['Loading file ' fileSave2]);
    load(fileSave2);
    NI_Data_OriTunedElecs = NI_Data;
else
    % get Data all Session for particular monkey or both combined for 
    % ori-tuned Electrodes
    [erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,~] = ...
    getData(folderSourceString,fileNameStringListAll,electrodeList_OriTuned,timeRangeParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); %#ok<ASGLU>
    save(fileSave2,'erpData','firingRateData','fftData','energyData','oriTuningData','NI_Data')
end

% Plotting Functions
plotData_spikes(hPlotsFig2,firingRateData,normalizateSpikeDataFlag) % spikes for static gratings from ori-selective electrodes, Fig 1
rescaleData(hPlotsFig2.hPlot1,-0.1,0.5,getYLims(hPlotsFig2.hPlot1),14);
rescaleData(hPlotsFig2.hPlot2(3),0,50,getYLims(hPlotsFig2.hPlot2(3)),14);
FigName2 = fullfile(folderSave_Figs,['Figure 2_' monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_oriTunedElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID)] );
saveas(hFigure2,[FigName2 '.fig'])
saveas(hFigure2,[FigName2,'.tif'])
% rescaleData(hPlotsFig3.hPlot1,0,250,getYLims(hPlotsFig3.hPlot1),12);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

removeERPFlag = 0; % Evoked Response
LFPdataProcessingMethod = 'Evoked Response';

fileSave3 = fullfile(folderSave,[monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) '.mat']); 
if exist(fileSave3,'file')
    disp(['Loading file ' fileSave3]);
    load(fileSave3);
    NI_Data_allElecsEvoked = NI_Data;
else
    oriSelectiveFlag = 0;
    % get Data all Session for particular monkey or both combined for all
    % Electrodes
    [erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,~] = ...
    getData(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); %#ok<ASGLU>
    save(fileSave3,'erpData','firingRateData','fftData','energyData','oriTuningData','NI_Data')
end

plotData_SSVEP(hPlotsFig4,energyData) % SSVEP Evoked, Fig 4;
rescaleData(hPlotsFig4.hPlot1,0,24,getYLims(hPlotsFig4.hPlot1),14);
rescaleData(hPlotsFig4.hPlot1(1,:),0,24,getYLims(hPlotsFig4.hPlot1(1,:)),14);
rescaleData(hPlotsFig4.hPlot1(2,:),0,24,getYLims(hPlotsFig4.hPlot1(2,:)),14);
rescaleData(hPlotsFig4.hPlot2(3),0,50,getYLims(hPlotsFig4.hPlot2(3)),14);
FigName4 = fullfile(folderSave_Figs,['Figure 4_' monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) ]);

saveas(hFigure4,[FigName4 '.fig'])
saveas(hFigure4,[FigName4,'.tif'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotData_histogram(hPlotsFig5, NI_Data_allElecsInduced,NI_Data_OriTunedElecs,NI_Data_allElecsEvoked)
FigName5 = fullfile(folderSave_Figs,['Figure 5_' monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff)...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) ]);

saveas(hFigure5,[FigName5 '.fig'])
saveas(hFigure5,[FigName5,'.tif'])
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plotData_spikes(hPlot,data,NormalizeDataFlag)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;

if NormalizeDataFlag
    % Normalize ERP and spike Data 
    %(fft or energy data need not be normalized because of log conversion
    data = normalizeData(data);    
end

% mean PSTH data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x
% timeVals
psthData = squeeze(mean(data.data,1));

% spike rate data: con_Ori2 (rows) x con_Ori2 (columns)
spikeRateDataST = squeeze(mean(data.analysisDataST,1));
spikeRateDataBL = squeeze(mean(data.analysisData_cBL,1));

% mean spike rate for two orthogonal pairs gratings when presented alone (or when the contrast of the other grating is 0%) 
for iElec = 1:size(data.analysisDataST,1)
    avgSpikeRateDataST_elecwise(iElec,:) = (squeeze(data.analysisDataST(iElec,:,end,:))' + flip(squeeze(data.analysisDataST(iElec,:,:,1)))')/2;
    sumSpikeRateDataST_elecwise(iElec,:) = squeeze(data.analysisDataST(iElec,:,end,:))' + flip(squeeze(data.analysisDataST(iElec,:,:,1)))';
end

avgSpikeRateDataST = mean(avgSpikeRateDataST_elecwise,1);
sem_avgSpikeRateDataST = std(avgSpikeRateDataST_elecwise,[],1)./sqrt(size(avgSpikeRateDataST_elecwise,1));

sumSpikeRateDataST = mean(sumSpikeRateDataST_elecwise,1);
sem_sumSpikeRateDataST = std(sumSpikeRateDataST_elecwise,[],1)./sqrt(size(sumSpikeRateDataST_elecwise,1));

% avgSpikeRateDataST = (spikeRateDataST(end,:)+ flip(spikeRateDataST(:,1))')/2;

diff_spikeRateData = spikeRateDataST - spikeRateDataBL;

sem_spikeRate = squeeze(std(squeeze(data.analysisDataST),[],1)./sqrt(size(data.analysisDataST,1)));

% computing N.I. population
for iElec= 1:size(data.analysisDataST,1)
    clear spikeRateElecVals_absolute spikeRateElecVals_relative
    spikeRateElecVals_absolute =  squeeze(data.analysisDataST(iElec,1,:,:));
    spikeRateElecVals_relative =  squeeze(data.analysisDataST(iElec,1,:,:))-squeeze(data.analysisData_cBL(iElec,1,:,:));
    NI_population_spikeRateAbsolute(iElec) = spikeRateElecVals_absolute(1,5)/(((spikeRateElecVals_absolute(1,1)+spikeRateElecVals_absolute(5,5)))/2)-1;
    NI_population_spikeRateRelative(iElec) = spikeRateElecVals_relative(1,5)/(((spikeRateElecVals_relative(1,1)+spikeRateElecVals_relative(5,5)))/2)-1;
end

%
OutlierVals = [-15 15];
NI_population_outlier = find(NI_population_spikeRateAbsolute<OutlierVals(1) | NI_population_spikeRateAbsolute>OutlierVals(2));
NI_population_outlierVals = NI_population_spikeRateAbsolute(NI_population_outlier);
NI_population_spikeRateAbsolute = NI_population_spikeRateAbsolute(setdiff(1:length(NI_population_spikeRateAbsolute),NI_population_outlier));
fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
    num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )

% remove Outlier elecs (add as a function)
OutlierVals = [-15 15];
NI_population_outlier = find(NI_population_spikeRateRelative<OutlierVals(1) | NI_population_spikeRateRelative>OutlierVals(2));
NI_population_outlierVals = NI_population_spikeRateRelative(NI_population_outlier);
NI_population_spikeRateRelative = NI_population_spikeRateRelative(setdiff(1:length(NI_population_spikeRateRelative),NI_population_outlier));
fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
    num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )


% PSTH plots
colors = jet(length(cValsUnique));
cFlipped_Indices = flip(1:length(cValsUnique2)); 

for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlot.hPlot1(1,c_Ori1),data.timeVals,squeeze(psthData(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(1,c_Ori1),'on');
    end
end

set(hPlot.hPlot1(1),'XLim',[-0.1 0.5]);
set(hPlot.hPlot1(1),'YLim',[0 90]);
tickLengthPlot = 2*get(hPlot.hPlot1(1),'TickLength');
xlabel(hPlot.hPlot1(1),'Time (s)')
ylabel(hPlot.hPlot1(1),'Spike rate(spike/s)')
displayRange(hPlot.hPlot1,[0.2 0.4],getYLims(hPlot.hPlot1),'k');


for i = 1:length(cValsUnique)
    set(hPlot.hPlot1(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end

% Color coded Plots of Spike Rates
imagesc(spikeRateDataST,'parent',hPlot.hPlot2(1));
colorBar_absSpikeRate = colorbar(hPlot.hPlot2(1)); 
colorYlabelHandle = get(colorBar_absSpikeRate,'Ylabel');
set(colorYlabelHandle,'String','Absolute Spike Rate (spikes/s)','fontSize',14);
plotPos = get(hPlot.hPlot2(1),'Position');
set(hPlot.hPlot2(1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
title(hPlot.hPlot2(1),['Mean NI: ',num2str(round(mean(NI_population_spikeRateAbsolute),2))],'fontWeight','bold');
% caxis(hPlot.hPlot2(1),[0 20]);
set(hPlot.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(1),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(1),'Contrast of Ori 2(%)');



imagesc(diff_spikeRateData,'parent',hPlot.hPlot2(2));
colorBar_rlvSpikeRate = colorbar(hPlot.hPlot2(2)); 
colorYlabelHandle = get(colorBar_rlvSpikeRate,'Ylabel');
set(colorYlabelHandle,'String','Change in Spike Rate (spikes/s)','fontSize',14);
plotPos = get(hPlot.hPlot2(2),'Position');
set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
title(hPlot.hPlot2(2),['Mean NI: ',num2str(round(mean(NI_population_spikeRateRelative),2))],'fontWeight','bold');
% caxis(hPlot.hPlot2(2),[0 20]);
set(hPlot.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(2),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(2),'Contrast of Ori 2(%)');


% CRF
errorbar(cValsUnique,spikeRateDataST(end,:),sem_spikeRate(end,:),...
    'Marker','o','LineWidth',2,'color',colors(end,:,:),'parent',hPlot.hPlot2(3))
hold(hPlot.hPlot2(3),'on');
errorbar(cValsUnique,flip(spikeRateDataST(:,1)),sem_spikeRate(:,1),...
    'Marker','o','LineWidth',2,'color',colors(1,:,:),'parent',hPlot.hPlot2(3))
errorbar(cValsUnique,diag(flipud(spikeRateDataST)),diag(flipud(sem_spikeRate)),'Marker','o','LineWidth',2,'color','k','parent',hPlot.hPlot2(3));
errorbar(cValsUnique,avgSpikeRateDataST,sem_avgSpikeRateDataST,'Marker','o','LineStyle','--','LineWidth',2,'color',[0.5 0.5 0.5],'parent',hPlot.hPlot2(3));
errorbar(cValsUnique,sumSpikeRateDataST,sem_sumSpikeRateDataST,'Marker','o','LineStyle','--','LineWidth',2,'color',[0.8 0.8 0.8],'parent',hPlot.hPlot2(3));


hold(hPlot.hPlot2(3),'off');
text(0.5,0.2,'cOri 2: 0%','color',colors(end,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.5,0.1,'cOri 1 = cOri 2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
set(hPlot.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlot.hPlot2(3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(3),'Absolute Spike Rate (spike/s)');

% axis(hPlot.hPlot2,'square');
% fix symmetry of axes boundary
plotPos = get(hPlot.hPlot1(3),'Position');
plotPos2 = get(hPlot.hPlot2(2),'Position');
set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos2(2) plotPos2(3) plotPos2(4)]);

plotPos = get(hPlot.hPlot1(5),'Position');
plotPos2 = get(hPlot.hPlot2(3),'Position');
set(hPlot.hPlot2(3),'Position',[plotPos(1) plotPos2(2) plotPos2(3)-(plotPos(1)-plotPos2(1)) plotPos2(4)]);


end
function plotData_energy(hPlot,data)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;
num_freqRanges = length(data.analysisDataST)-1; % alpha, gamma, hi-gamma Induced power is being plotted

% mean energy data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x
% freqVals
energyVsFrequencyDataST = squeeze(mean(data.dataST(:,1,:,:,:),1));
energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,1,:,:,:),1));
dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);

% PSD plots
colors = jet(length(cValsUnique));
cFlipped_Indices = flip(1:length(cValsUnique2)); 

for c_Ori1 = 1:length(cValsUnique)
    plot(hPlot.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
    hold(hPlot.hPlot1(1,c_Ori1),'on');
    plot(hPlot.hPlot1(2,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)-energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
    hold(hPlot.hPlot1(2,c_Ori1),'on');
end

for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlot.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataST(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(1,c_Ori1),'on');
        plot(hPlot.hPlot1(2,c_Ori1),data.freqVals,squeeze(dEnergyVsFrequencyData(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(2,c_Ori1),'on');
    end
end




% set(hPlot.hPlot1(1,1),'XLim',[0 250]);
set(hPlot.hPlot1(1,1),'YLim',[-1.5 3.5]);
set(hPlot.hPlot1(2,1),'YLim',[-4 10]);
set(hPlot.hPlot1(1,1),'XLim',[0 250]);
set(hPlot.hPlot1(2,1),'XLim',[0 250]);
tickLengthPlot = 2*get(hPlot.hPlot1(1),'TickLength');

% displayRange(hPlot.hPlot1,[0.2 0.4],getYLims(hPlot.hPlot1),'k');


for i = 1:length(cValsUnique)
    set(hPlot.hPlot1(1,i),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlot.hPlot1(2,i),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
xlabel(hPlot.hPlot1(2,1),'Frequency (Hz)','fontSize',12)

ylabel(hPlot.hPlot1(1,1),'log_1_0 (Power)','fontSize',12)
ylabel(hPlot.hPlot1(2,1),{'Change in','Power (dB)'},'fontSize',12)

% energy data: con_Ori2 (rows) x con_Ori2 (columns)
for i = 1: num_freqRanges

    clear energyDataST energyDataBL NI_population_energyAbsolute NI_population_energyRelative
    energyDataST = squeeze(mean(data.analysisDataST{i},1));
    energyDataBL = squeeze(mean(data.analysisData_cBL{i},1));

%     sem_EnergyDataST = squeeze(std(squeeze(data.analysisDataST{i}),[],1)./sqrt(size(data.analysisDataST{i},1)));
    dEnergyData = 10*(energyDataST - energyDataBL); %across elecs
    sem_dEnergyData = squeeze(std(10*(data.analysisDataST{i}-data.analysisDataBL{i}),[],1)./sqrt(size(data.analysisDataST{i},1)));
    
    % computing N.I. population
    for iElec= 1:size(data.analysisDataST{i},1)
        clear spikeRateElecVals_absolute spikeRateElecVals_relative
        energyData_Elec_absolute =  squeeze(data.analysisDataST{i}(iElec,:,:));
        energyData_Elec_relative =  10*(squeeze(data.analysisDataST{i}(iElec,:,:))-squeeze(data.analysisData_cBL{i}(iElec,:,:)));
        NI_population_energyAbsolute(iElec) = energyData_Elec_absolute(1,5)/(((energyData_Elec_absolute(1,1)+energyData_Elec_absolute(5,5)))/2)-1;
        NI_population_energyRelative(iElec) = energyData_Elec_relative(1,5)/(((energyData_Elec_relative(1,1)+energyData_Elec_relative(5,5)))/2)-1;
    end

    OutlierVals = [-15 15];
    NI_population_outlier = find(NI_population_energyAbsolute<OutlierVals(1) | NI_population_energyAbsolute>OutlierVals(2));
    NI_population_outlierVals = NI_population_energyAbsolute(NI_population_outlier);
    NI_population_energyAbsolute = NI_population_energyAbsolute(setdiff(1:length(NI_population_energyAbsolute),NI_population_outlier));
    fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
        num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )

    % remove Outlier elecs (add as a function)
    OutlierVals = [-15 15];
    NI_population_outlier = find(NI_population_energyRelative<OutlierVals(1) | NI_population_energyRelative>OutlierVals(2));
    NI_population_outlierVals = NI_population_energyRelative(NI_population_outlier);
    NI_population_energyRelative = NI_population_energyRelative(setdiff(1:length(NI_population_energyRelative),NI_population_outlier));
    fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
        num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )
    % remove Outlier elecs (add as a function)
%     OutlierVals = [-15 15];
%     NI_population_outlier = find(NI_population_spikeRateRelative<OutlierVals(1) | NI_population_spikeRateRelative>OutlierVals(2));
%     NI_population_outlierVals = NI_population_spikeRateRelative(NI_population_outlier);
%     NI_population_spikeRateRelative = NI_population_spikeRateRelative(setdiff(1:length(NI_population_spikeRateRelative),NI_population_outlier));
%     fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
%         num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )



    % Color coded Plots of energyData

    imagesc(energyDataST,'parent',hPlot.hPlot2(i,1));
    imagesc(dEnergyData,'parent',hPlot.hPlot2(i,2));

    
    % CRF
    errorbar(cValsUnique,dEnergyData(end,:),sem_dEnergyData(end,:),...
        'Marker','o','LineWidth',2,'color',colors(end,:,:),'parent',hPlot.hPlot2(i,3))
    hold(hPlot.hPlot2(i,3),'on');
    errorbar(cValsUnique,diag(flipud(dEnergyData)),diag(flipud(sem_dEnergyData)),'Marker','o','LineWidth',2,'color','k','parent',hPlot.hPlot2(i,3));
    hold(hPlot.hPlot2(i,3),'off');
    
    % Setting axes properties
    colorBar_absPSD = colorbar(hPlot.hPlot2(i,1)); 
    colorYlabelHandle = get(colorBar_absPSD,'Ylabel');
    set(colorYlabelHandle,'String','log_1_0(Power)','fontSize',12);
    plotPos = get(hPlot.hPlot2(i,1),'Position');
    set(hPlot.hPlot2(i,1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.0450 plotPos(4)]);
    title(hPlot.hPlot2(i,1),['Mean NI: ',num2str(round(mean(NI_population_energyAbsolute),2))],'fontSize',12,'fontWeight','bold');
    % caxis(hPlot.hPlot2(1),[0 4]);
    set(hPlot.hPlot2(i,1),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    if i == 3
        set(hPlot.hPlot2(i,1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
    else
        set(hPlot.hPlot2(i,1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',[],'YTickLabel',[]);
    end

    if i == 3
        xlabel(hPlot.hPlot2(i,1),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(i,1),'Contrast of Ori 2(%)');
    end


    colorBar_rlvPSD = colorbar(hPlot.hPlot2(i,2)); 
    colorYlabelHandle = get(colorBar_rlvPSD,'Ylabel');
    set(colorYlabelHandle,'String',{'\Delta Power (dB)'},'fontSize',12);
    plotPos = get(hPlot.hPlot2(i,2),'Position');
    set(hPlot.hPlot2(i,2),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.0450 plotPos(4)]);
    title(hPlot.hPlot2(i,2),['Mean NI: ',num2str(round(mean(NI_population_energyRelative),2))],'fontSize',12,'fontWeight','bold');
    % caxis(hPlot.hPlot2(2),[0 4]);
    set(hPlot.hPlot2(i,2),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    if i==3
        set(hPlot.hPlot2(i,2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
    else
        set(hPlot.hPlot2(i,2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',[],'YTickLabel',[]);
    end
    
    if i == 3
        xlabel(hPlot.hPlot2(i,2),'Contrast of Ori 1(%)','fontSize',12);ylabel(hPlot.hPlot2(i,2),'Contrast of Ori 2(%)','fontSize',12);
    end

    % % CRF
    if i == 3
        text(0.2,0.2,'cOri 2: 0%','color',colors(end,:,:),'fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(i,3))
        text(0.2,0.1,'cOri 1 = cOri 2','color','k','fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(i,3))
    end
    set(hPlot.hPlot2(i,3),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    if i==3
        set(hPlot.hPlot2(i,3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
    else
        set(hPlot.hPlot2(i,3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',[]);
    end
    
    if i == 3
        xlabel(hPlot.hPlot2(i,3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(i,3),'Change in Power(dB)');
    end

    % axis(hPlot.hPlot2,'square');
    % fix symmetry of axes boundary
%     plotPos = get(hPlot.hPlot1(2,3),'Position');
%     plotPos2 = get(hPlot.hPlot2(2),'Position');
%     set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos2(2) plotPos2(3) plotPos2(4)]);
% 
    plotPos = get(hPlot.hPlot2(i,2),'Position');
    plotPos2 = get(hPlot.hPlot2(i,3),'Position');
    set(hPlot.hPlot2(i,3),'Position',[plotPos2(1)-(plotPos(3)-plotPos2(3)) plotPos2(2) plotPos(3) plotPos2(4)]);
end
end
function plotData_SSVEP(hPlot,data)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;
% num_freqRanges = length(data.analysisDataST)-1; % alpha, gamma, hi-gamma Induced power is being plotted

% mean energy data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x
% freqVals
energyVsFrequencyDataST = squeeze(mean(data.dataST(:,2,:,:,:),1));
energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,2,:,:,:),1));
dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);

% PSD plots
colors = jet(length(cValsUnique));
cFlipped_Indices = flip(1:length(cValsUnique2)); 

for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlot.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataST(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(1,c_Ori1),'on');
        plot(hPlot.hPlot1(2,c_Ori1),data.freqVals,squeeze(dEnergyVsFrequencyData(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(2,c_Ori1),'on');
    end
end

for c_Ori1 = 1:length(cValsUnique)
plot(hPlot.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
plot(hPlot.hPlot1(2,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)-energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
end


% energy data: con_Ori2 (rows) x con_Ori2 (columns)

clear energyDataST energyDataBL NI_population_energyAbsolute NI_population_energyRelative
energyDataST = squeeze(mean(data.analysisDataST{4},1));
energyDataBL = squeeze(mean(data.analysisData_cBL{4},1));

%     sem_EnergyDataST = squeeze(std(squeeze(data.analysisDataST{i}),[],1)./sqrt(size(data.analysisDataST{i},1)));
dEnergyData = 10*(energyDataST - energyDataBL); %across elecs
sem_dEnergyData = squeeze(std(10*(data.analysisDataST{4}-data.analysisDataBL{4}),[],1)./sqrt(size(data.analysisDataST{4},1)));

% computing N.I. population
for iElec= 1:size(data.analysisDataST{4},1)
    clear spikeRateElecVals_absolute spikeRateElecVals_relative
    energyData_Elec_absolute =  squeeze(data.analysisDataST{4}(iElec,:,:));
    energyData_Elec_relative =  10*(squeeze(data.analysisDataST{4}(iElec,:,:))-squeeze(data.analysisData_cBL{4}(iElec,:,:)));
    NI_population_energyAbsolute(iElec) = energyData_Elec_absolute(1,5)/(((energyData_Elec_absolute(1,1)+energyData_Elec_absolute(5,5)))/2)-1;
    NI_population_energyRelative(iElec) = energyData_Elec_relative(1,5)/(((energyData_Elec_relative(1,1)+energyData_Elec_relative(5,5)))/2)-1;
end

% remove Outlier elecs (add as a function)
%     OutlierVals = [-15 15];
%     NI_population_outlier = find(NI_population_spikeRateRelative<OutlierVals(1) | NI_population_spikeRateRelative>OutlierVals(2));
%     NI_population_outlierVals = NI_population_spikeRateRelative(NI_population_outlier);
%     NI_population_spikeRateRelative = NI_population_spikeRateRelative(setdiff(1:length(NI_population_spikeRateRelative),NI_population_outlier));
%     fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
%         num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )
    OutlierVals = [-15 15];
    NI_population_outlier = find(NI_population_energyAbsolute<OutlierVals(1) | NI_population_energyAbsolute>OutlierVals(2));
    NI_population_outlierVals = NI_population_energyAbsolute(NI_population_outlier);
    NI_population_energyAbsolute = NI_population_energyAbsolute(setdiff(1:length(NI_population_energyAbsolute),NI_population_outlier));
    fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
        num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )

    % remove Outlier elecs (add as a function)
    OutlierVals = [-15 15];
    NI_population_outlier = find(NI_population_energyRelative<OutlierVals(1) | NI_population_energyRelative>OutlierVals(2));
    NI_population_outlierVals = NI_population_energyRelative(NI_population_outlier);
    NI_population_energyRelative = NI_population_energyRelative(setdiff(1:length(NI_population_energyRelative),NI_population_outlier));
    fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
        num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )



% Color coded Plots of energyData

imagesc(energyDataST,'parent',hPlot.hPlot2(1));
imagesc(dEnergyData,'parent',hPlot.hPlot2(2));


% CRF
errorbar(cValsUnique,dEnergyData(end,:),sem_dEnergyData(end,:),...
    'Marker','o','LineWidth',2,'color',colors(end,:,:),'parent',hPlot.hPlot2(3))
hold(hPlot.hPlot2(3),'on');
errorbar(cValsUnique,diag(flipud(dEnergyData)),diag(flipud(sem_dEnergyData)),'Marker','o','LineWidth',2,'color','k','parent',hPlot.hPlot2(3));
hold(hPlot.hPlot2(3),'off');

% Seeting axes properties
% set(hPlot.hPlot1(1),'XLim',[0 24]);
% set(hPlot.hPlot1(2),'XLim',[0 24]);
% set(hPlot.hPlot1(1,1),'YLim',[-1.5 3.5]);
% set(hPlot.hPlot1(2,1),'YLim',[-4 12]);
% set(hPlot.hPlot1(1),'YLim',[-1 4]);
tickLengthPlot = 2*get(hPlot.hPlot1(1),'TickLength');
xlabel(hPlot.hPlot1(2,1),'Frequency (Hz)')
ylabel(hPlot.hPlot1(1,1),'log_1_0(Power)')
ylabel(hPlot.hPlot1(2,1),{'Change in','Power(dB)'})
displayRange(hPlot.hPlot1(1,:),[16 16],getYLims(hPlot.hPlot1(1,:)),'k');
displayRange(hPlot.hPlot1(2,:),[16 16],getYLims(hPlot.hPlot1(2,:)),'k');
for i = 1:length(cValsUnique)
    set(hPlot.hPlot1(1,i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlot.hPlot1(2,i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end

% Color coded Plots of Spike Rates

colorBar_absPSD = colorbar(hPlot.hPlot2(1)); 
colorYlabelHandle = get(colorBar_absPSD,'Ylabel');
set(colorYlabelHandle,'String','log_1_0(Power)','fontSize',14);
plotPos = get(hPlot.hPlot2(1),'Position');
set(hPlot.hPlot2(1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
title(hPlot.hPlot2(1),['Mean NI: ',num2str(round(mean(NI_population_energyAbsolute),2))],'fontWeight','bold');
% caxis(hPlot.hPlot2(1),[0 4]);
set(hPlot.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(1),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(1),'Contrast of Ori 2(%)');


colorBar_rlvPSD = colorbar(hPlot.hPlot2(2)); 
colorYlabelHandle = get(colorBar_rlvPSD,'Ylabel');
set(colorYlabelHandle,'String','Change in Power(dB)','fontSize',14);
plotPos = get(hPlot.hPlot2(2),'Position');
set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
title(hPlot.hPlot2(2),['Mean NI: ',num2str(round(mean(NI_population_energyRelative),2))],'fontWeight','bold');
% caxis(hPlot.hPlot2(2),[0 4]);
set(hPlot.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(2),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(2),'Contrast of Ori 2(%)');


% % CRF
text(0.5,0.2,'cOri 2: 0%','color',colors(end,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.5,0.1,'cOri 1 = cOri 2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
set(hPlot.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlot.hPlot2(3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(3),'Change in Power(dB)');

% axis(hPlot.hPlot2,'square');
% fix symmetry of axes boundary
plotPos = get(hPlot.hPlot1(2,3),'Position');
plotPos2 = get(hPlot.hPlot2(2),'Position');
set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos2(2) plotPos2(3) plotPos2(4)]);

plotPos = get(hPlot.hPlot1(2,5),'Position');
plotPos2 = get(hPlot.hPlot2(3),'Position');
set(hPlot.hPlot2(3),'Position',[plotPos(1) plotPos2(2) plotPos2(3)-(plotPos(1)-plotPos2(1)) plotPos2(4)]);
end

function plotData_histogram(hPlotsFig5, NI_Data_allElecsInduced,NI_Data_OriTunedElecs,NI_Data_allElecsEvoked)
% 1st Row (Spike NI- allElecs) Column 1 & 2 gives NI of absolute and
% relative measures respectively
data = NI_Data_allElecsInduced.firingRate_ST;
histogram(hPlotsFig5.hPlot1(1,1),data,min(data):0.2:max(data));

data = NI_Data_allElecsInduced.dfiringRate;
histogram(hPlotsFig5.hPlot1(1,2),data,min(data):0.2:max(data));

% 2nd Row (Spike NI- oriTunedElecs) Column 1 & 2 gives NI of absolute and
% relative measures respectively
data = NI_Data_OriTunedElecs.firingRate_ST;
histogram(hPlotsFig5.hPlot1(2,1),data,min(data):0.2:max(data));
data = NI_Data_OriTunedElecs.dfiringRate;
histogram(hPlotsFig5.hPlot1(2,2),data,min(data):0.2:max(data));

% 3rd Row (alpha power induced- allElecs) Column 1 & 2 gives NI of absolute and
% relative measures respectively
data = NI_Data_allElecsInduced.energy_ST{1};
histogram(hPlotsFig5.hPlot1(3,1),data);
data = NI_Data_OriTunedElecs.denergy{1};
histogram(hPlotsFig5.hPlot1(3,2),data,min(data):0.2:max(data));

% 4th Row (gamma power induced- allElecs) Column 1 & 2 gives NI of absolute and
% relative measures respectively
data = NI_Data_allElecsInduced.energy_ST{2};
histogram(hPlotsFig5.hPlot1(4,1),data,min(data):0.2:max(data));
data = NI_Data_allElecsInduced.denergy{2};
histogram(hPlotsFig5.hPlot1(4,2),data,min(data):0.2:max(data));


% 5th Row (hi-gamma power induced- allElecs) Column 1 & 2 gives NI of absolute and
% relative measures respectively
data = NI_Data_allElecsInduced.energy_ST{3};
histogram(hPlotsFig5.hPlot1(5,1),data);
data = NI_Data_allElecsInduced.denergy{3};
histogram(hPlotsFig5.hPlot1(5,2),data,min(data):0.2:max(data));

% 6th Row (SSVEP power evoked- allElecs) Column 1 & 2 gives NI of absolute and
% relative measures respectively
data = NI_Data_allElecsEvoked.energy_ST{4};
histogram(hPlotsFig5.hPlot1(6,1),data);
data = NI_Data_allElecsEvoked.denergy{4};
histogram(hPlotsFig5.hPlot1(6,2),data,min(data):0.2:max(data));

title(hPlotsFig5.hPlot1(1,1),'NI-Absolute')
title(hPlotsFig5.hPlot1(1,2),'NI-Relative')
ylabel(hPlotsFig5.hPlot1(1,1),'Spikes')
ylabel(hPlotsFig5.hPlot1(2,1),'Spikes-OriElec')
ylabel(hPlotsFig5.hPlot1(3,1),'alpha')
ylabel(hPlotsFig5.hPlot1(4,1),'gamma')
ylabel(hPlotsFig5.hPlot1(5,1),'hi-gamma')
ylabel(hPlotsFig5.hPlot1(6,1),'SSVEP')
end


% Get Color String
% function [colorString, colorNames] = getColorString
% 
% colorNames = 'brkgcmy';
% colorString = 'blue|red|black|green|cyan|magenta|yellow';
% 
% end

% Get FileNamesList
function fileNameStringListAll = getFileNameStringList(monkeyName,gridType)

if strcmp(monkeyName,'all')
    clear monkeyName
    monkeyName{1} = 'alpaH'; monkeyName{2} = 'kesariH';
else
    monkeyName = mat2cell(monkeyName,1);
end

fileNameStringList = cell(1,2);
for i=1:size(monkeyName,2)
    clear expDates protocolNames
    [expDates,protocolNames,~]= dataInformationPlaidNorm(monkeyName{i},gridType,0);
    numSessions = size(protocolNames,2);
    tmpFileNameList = cell(1,numSessions);
    for j = 1:numSessions
    tmpFileNameList{j} = [monkeyName{i} expDates{j} protocolNames{j}];
    end
    fileNameStringList{i} = tmpFileNameList;
end

pos=1;
for i=1:size(monkeyName,2)
    for j=1:length(fileNameStringList{i})
        fileNameStringListAll{pos} = fileNameStringList{i}{j}; %#ok<*AGROW>
        pos=pos+1;
    end
end
end

% Get ElectrodesList
function ElectrodeArrayListAll = getElectrodesList(fileNameStringTMP,oriSelectiveFlag,folderSourceString)

% [~,tmpElectrodeArrayList,~] = getGoodElectrodesDetails(fileNameStringTMP,oriSelectiveFlag,folderSourceString);

gridType = 'microelectrode';

numSessions = length(fileNameStringTMP);
tmpElectrodeStringList = cell(1,numSessions);
tmpElectrodeArrayList = cell(1,numSessions);
numElecs = 0;

Monkey1_ExpDates = dataInformationPlaidNorm('alpaH',gridType,0);
Monkey1_SessionNum = length(Monkey1_ExpDates);
% Monkey2_ExpDates = dataInformationPlaidNorm('kesariH',gridType,0);
% Monkey2_SessionNum = length(Monkey2_ExpDates);

for i = 1:numSessions
    clear monkeyName
    if strcmp(fileNameStringTMP{i}(1:5),'alpaH') 
        monkeyName = 'alpaH';
        expDate = fileNameStringTMP{i}(6:11);
        protocolName = fileNameStringTMP{i}(12:end);
    elseif strcmp(fileNameStringTMP{i}(1:7),'kesariH')
        monkeyName = 'kesariH';
        expDate = fileNameStringTMP{i}(8:13);
        protocolName = fileNameStringTMP{i}(14:end);
    end
    if i == 1
        disp(['MonkeyName: ' ,monkeyName])
    elseif i == Monkey1_SessionNum+1 % 13 Sessions are from alpaH; 9 Sessions from kesariH;
        disp(['MonkeyName: ' ,monkeyName])
    end
    versionNum = 2;
    [tmpElectrodeStringList{i},tmpElectrodeArrayList{i},goodElectrodes] = getGoodElectrodesSingleSession(monkeyName,expDate,protocolName,gridType,folderSourceString,oriSelectiveFlag,versionNum);
    numElecs = numElecs+length(goodElectrodes);
%     allElecs = numElecs;
end
% if length(tmpElectrodeStringList)> 1
%    clear tmpElectrodeStringList
%    tmpElectrodeStringList = {['all (N=' num2str(allElecs) ')']};
% end

% ElectrodeStringListAll = tmpElectrodeStringList;
ElectrodeArrayListAll = tmpElectrodeArrayList;
end

% Normalize data for ERP and Spike data
function normData = normalizeData(x)
for iElec = 1:size(x.data,1)
    for t = 1:size(x.data,2)
        normData.data(iElec,t,:,:,:) = x.data(iElec,t,:,:,:)./max(max(max(abs(x.data(iElec,t,:,:,:)))));
        normData.analysisDataBL(iElec,t,:,:) = x.analysisDataBL(iElec,t,:,:)./max(max(abs(x.analysisDataBL(iElec,t,:,:))));
        normData.analysisData_cBL(iElec,t,:,:) = x.analysisData_cBL(iElec,t,:,:)./max(max(abs(x.analysisData_cBL(iElec,t,:,:))));
        normData.analysisDataST(iElec,t,:,:) = x.analysisDataST(iElec,t,:,:)./max(max(abs(x.analysisDataST(iElec,t,:,:))));
        normData.timeVals = x.timeVals;
        normData.N = x.N;
    end
end
end

% Get data for a single electrode in a selected session
% function data = getDataSingleElec(data,electrodeNum,analysisMeasure)
%     if analysisMeasure == 1 || analysisMeasure ==2
%     data.data = data.data(electrodeNum,:,:,:,:);
%     data.analysisDataBL = data.analysisDataBL(electrodeNum,:,:,:);
%     data.analysisDataST = data.analysisDataST(electrodeNum,:,:,:);
%     
%     elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6||analysisMeasure == 7
%         data.dataBL = data.dataBL(electrodeNum,:,:,:,:);
%         data.dataST = data.dataST(electrodeNum,:,:,:,:);
%         for i = 1:length(data.analysisDataST)
%             data.analysisDataBL{i} = data.analysisDataBL{i}(electrodeNum,:,:,:);
%             data.analysisDataST{i} = data.analysisDataST{i}(electrodeNum,:,:,:);
%         end
%     end
% end


% Draw lines for timeTange or FreqRange
function displayRange(plotHandles,range,yLims,colorName)
[nX,nY] = size(plotHandles);
%yLims = getYLims(plotHandles);

yVals = yLims(1):(yLims(2)-yLims(1))/100:yLims(2);
xVals1 = range(1) + zeros(1,length(yVals));
xVals2 = range(2) + zeros(1,length(yVals));

for i=1:nX
    for j=1:nY
        hold(plotHandles(i,j),'on');
        plot(plotHandles(i,j),xVals1,yVals,'color',colorName);
        plot(plotHandles(i,j),xVals2,yVals,'color',colorName);
    end
end

end

% get Y limits for an axis
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end

% Rescale data
function rescaleData(plotHandles,xMin,xMax,yLims,labelSize)

[numRows,numCols] = size(plotHandles);
% labelSize=14;
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j==1
                set(plotHandles(i,j),'fontSize',labelSize);
            elseif j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        elseif (rem(i,2)==0 && j~=1)
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        elseif (rem(i,2)==1 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        elseif (rem(i,2)==1 && j~=1)
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        else 
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

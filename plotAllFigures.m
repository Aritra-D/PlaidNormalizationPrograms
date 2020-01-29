% This program displays all the figures for this Plaid Normalization paper
% Response matrix is drawn in such a way that rows correspond to increasing
% contrasts of component Grating 2 in positive y direction and columns
% correspond to contrasts of component Grating 1 in positive x direction

% monkeyName can be either 'alpaH','kesariH','all'
% FolderSourceString is the path of all Extracted dataset
% timeRangeForComputation is the Stimulus period for which analysis is
% being performed
% colorScheme = 'color' or 'grayScale'.

function plotAllFigures(monkeyName,folderSourceString,timeRangeForComputation,combineUniqueElectrodeData,colorScheme)
if ~exist('folderSourceString','var')
    if strcmp(getenv('username'),'Aritra') || strcmp(getenv('username'),'Lab Computer-Aritra')
        folderSourceString = 'E:\data\PlaidNorm\';
    elseif strcmp(getenv('username'),'Supratim Ray')
        folderSourceString = 'M:\CommonData\Non_Standard\PlaidNorm\';
    end
end

close all; % closes any open figure to avoid any overlaying issues

% Variable Parameters
elecParams.spikeCutoff = 15;
elecParams.snrCutoff = 2;
elecParams.dRange = [0 0.75];
elecParams.unitID = 0;
elecParams.oriSelectiveFlag = 0;
elecParams.getSpikeElectrodesFlag = 1;
tapers_MT = [1 1]; % parameters for MT analysis
removeERPFlag = 0;

% Fixed parameters
gridType = 'Microelectrode';
normalizateSpikeDataFlag = 0;
timeRangeForComputationBL = -0.05+[-diff(timeRangeForComputation) 0];
folderSourceString_Project = strtok(folderSourceString,'\');

% path for Figures to be saved
folderSave = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures');
if ~exist(folderSave,'dir')
    mkdir(folderSave)
end

freqRanges{1} = [8 12]; % alpha
freqRanges{2} = [30 80]; % gamma
freqRanges{3} = [104 250]; % hi-gamma
freqRanges{4} = [16 16];  % SSVEP

timeRangeParameters.blRange = timeRangeForComputationBL;
timeRangeParameters.stRange = timeRangeForComputation;
timeRangeParameters.erpRange = [0.05 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% display properties %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 1: Example Single electrode PSTH data for a single session:
hFigure1 = figure(1);
set(hFigure1,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig1.hPlot1 = getPlotHandles(5,5,[0.2 0.4 0.55 0.55],0.01,0.01,0); linkaxes(hPlotsFig1.hPlot1);
hPlotsFig1.hPlot2 = getPlotHandles(2,5,[0.2 0.08 0.55 0.22],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

textH{1} = getPlotHandles(1,1,[0.12 0.97 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.12 0.23 0.01 0.01]);

textString = {'A','B'};
for i = 1:2
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end

% % Figure 2: Example Single electrode TF data for a single session:
hFigure2 = figure(2);
set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig2.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.65],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot1);
hPlotsFig2.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.13],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot2);

textH{1} = getPlotHandles(1,1,[0.12 0.97 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.12 0.23 0.01 0.01]);

textString = {'A','B'};
for i = 1:2
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end

% Figure 3 (Spike data from Orientation Selective Elecs-
% Pref (x-axis) Null (y-axis) axes); spike data along
% increasing contrasts of Pref Ori with contrasts of  null Ori presented
% by different colors; absolute color map 5x5; relative color map 5x5;
% CRF (no Norm, Pref,Pref+Null, Avg, Null)
hFigure3 = figure(3);
set(hFigure3,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig3.hPlot1 = getPlotHandles(1,5,[0.15 0.65 0.7 0.2],0.01,0.01,1); linkaxes(hPlotsFig3.hPlot1);
hPlotsFig3.hPlot2 = getPlotHandles(1,3,[0.15 0.23 0.7 0.2578],0.1188,0.05,0);

textH{1} = getPlotHandles(1,1,[0.1 0.9 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.1 0.5378 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.38 0.5378 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.65 0.5378 0.01 0.01]);

textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end

% % Figure 4 (PSDs (A); deltaPSDs (B) along increasing contrast of Ori 1
% % with contrasts of Ori 2 presented in different colors;
% % absolute color map 5x5; relative color map 5x5;
% % CRF (no Norm, Pref,Pref+Null, Avg, Null)
% % for alpha (C), gamma (D) and high-gamma (E) for ERP-subtracted data
hFigure4 = figure(4);
set(hFigure4,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig4.hPlot1 = getPlotHandles(2,5,[0.25 0.69 0.5 0.28],0.01,0.01,1);
linkaxes(hPlotsFig4.hPlot1(1,:)); linkaxes(hPlotsFig4.hPlot1(2,:));
hPlotsFig4.hPlot2 = getPlotHandles(2,3,[0.25 0.08 0.5 0.51],0.07,0.07,1);

textH{1} = getPlotHandles(1,1,[0.2 0.97 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.2 0.6 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.4 0.6 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.59 0.6 0.01 0.01]);

textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end


%
% % Figure 5 (PSDs (A); deltaPSDs (B) along increasing contrast of Ori 1
% % with contrasts of Ori 2 presented in different colors;
% % absolute color map 5x5; relative color map 5x5;
% % CRF (no Norm, Pref,Pref+Null, Avg, Null)
% % for SSVEP (C) for non-ERP subtracted (Evoked) data
hFigure5 = figure(5);
set(hFigure5,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig5.hPlot1 = getPlotHandles(2,5,[0.15 0.5 0.7 0.4],0.01,0.01,1);
linkaxes(hPlotsFig5.hPlot1(1,:));linkaxes(hPlotsFig5.hPlot1(2,:));
hPlotsFig5.hPlot2 = getPlotHandles(1,3,[0.15 0.1 0.7 0.2578],0.1188,0.05,0);
%
%
textH{1} = getPlotHandles(1,1,[0.1 0.9 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.1 0.36 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.38 0.36 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.65 0.36 0.01 0.01]);

textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
% % Figure 6: Comparison of NI,fitted parameters (alpha, sigma), percentage explained variance for spikeRate, gamma, high
% % gamma, SSVEP
hFigure6 = figure(6);
set(hFigure6,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig6.hPlot1 = getPlotHandles(2,2,[0.2 0.05 0.6 0.9],0.05,0.1,1);

textH{1} = getPlotHandles(1,1,[0.15 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.15 0.45 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.48 0.95 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.48 0.45 0.01 0.01]);

textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end


%%%%%%%%%%%%%%%%%%%%% Get Session Details for Monkey(s) %%%%%%%%%%%%%%%%%%%
fileNameStringListAll = getFileNameStringList(monkeyName,gridType);

%%%%%%%%%%%%%%%%%%%%% Get data for all Good Electrodes %%%%%%%%%%%%%%%%%%%%
disp('all Electrodes:')
[electrodeList_All,elecIDs,electrodeListString_All,elecInfo,numElecs] = getElectrodesList(fileNameStringListAll,spikeElecParams,timeRangeForComputation_spikes,folderSourceString);
disp([num2str(numElecs) ' Good Electrodes'])

cne = 0;
fileSave = fullfile(folderSave,[monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) ...
    '_allElecs_T' num2str(round(1000*timeRangeForComputation_spikes(1))) '_' num2str(round(1000*timeRangeForComputation_spikes(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(cne) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) '.mat']);

if exist(fileSave,'file')
    disp(['Loading file ' fileSave]);
    load(fileSave);
else
    % get Data all Session for monkey(s) for all Electrodes
    [erpData,firingRateData,fftData,energyData,energyDataTF,oriTuningData,NI_Data,~,N]  = ...
        getData(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,elecParams,removeERPFlag);
    save(fileSave,'elecInfo','erpData','firingRateData','fftData','energyData','energyDataTF','oriTuningData','NI_Data','N')
end

%%%%%%%%%%%%%%%%%%%%%% Process Dataset for plots %%%%%%%%%%%%%%%%%%%%%%%%%%
fileSave2 = fullfile(folderSave,[monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) ...
    '_allElecs_T' num2str(round(1000*timeRangeForComputation_spikes(1))) '_' num2str(round(1000*timeRangeForComputation_spikes(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) '.mat']);

if exist(fileSave,'file')
    disp(['Loading file ' fileSave]);
    load(fileSave);
else
[erpData_combined,firingRateData_combined,fftData_combined,energyData_combined,energyDataTF_combined,oriTuningData,N,elecIDs]  = ...
    combineData(erpData,firingRateData,fftData,energyData,energyDataTF,oriTuningData,N,elecParams);
    save(fileSave2,'elecIDs','erpData','firingRateData','fftData','energyData','energyDataTF','oriTuningData','NI_Data','N')
end

% Figures
if strcmp(colorScheme,'greyscale')||strcmp(colorScheme,'grayscale')
    folderSave_Figs = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\Figures\Grayscale\FitV2');
elseif strcmp(colorScheme,'color')
    folderSave_Figs = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\Figures\Color\FitV2');
end

if ~exist(folderSave_Figs,'dir')
    mkdir(folderSave_Figs)
end

FigName1 = fullfile(folderSaveFigs,['Figure 1_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_elec' num2str(elecNum)...
    '_T' num2str(round(1000*timeRangeParameters(1).stRange(1))) '_' num2str(round(1000*timeRangeParameters(1).stRange(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) ]);


FigName2 = fullfile(folderSaveFigs,['Figure 2_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_elec' num2str(elecNum)...
    '_T' num2str(round(1000*timeRangeParameters(2).stRange(1))) '_' num2str(round(1000*timeRangeParameters(2).stRange(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) ]);

elecNum = 29;
plotExampleElectrodeData(hPlotsFig1,hPlotsFig2,elecNum,elecInfo,firingRateData,energyData,energyDataTF,colorScheme)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combineUniqueElectrodeData = 1;
plotData_spikes(hPlotsFig3,firingRateData,timeRangeForComputation_spikes,elecInfo,normalizateSpikeDataFlag,combineUniqueElectrodeData,colorScheme) % spikes for static gratings, Fig 1

FigName3 = fullfile(folderSave_Figs,['Figure 3_' monkeyName '_N' num2str(spikeElecParams.spikeCutoff) '_S' num2str(spikeElecParams.snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation_spikes(1))) '_' num2str(round(1000*timeRangeForComputation_spikes(2))) ...
    '_d' num2str(spikeElecParams.dRange(1)) '_' num2str(spikeElecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(spikeElecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(spikeElecParams.unitID)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES 4 & 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotData_energy(hPlotsFig4,energyData,elecInfo,combineUniqueElectrodeData,colorScheme) % alpha, gamma, hi-gamma for static gratings, Fig 3;

set(findall(hFigure4,'-property','FontSize'),'FontSize',13);
FigName4 = fullfile(folderSave_Figs,['Figure 4_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) ]);


plotData_SSVEP(hPlotsFig5,energyData,elecInfo,combineUniqueElectrodeData,colorScheme) % SSVEP Evoked, Fig 4;
rescaleData(hPlotsFig5.hPlot1,0,24,getYLims(hPlotsFig5.hPlot1),14);
rescaleData(hPlotsFig5.hPlot1(1,:),0,24,getYLims(hPlotsFig5.hPlot1(1,:)),14);
rescaleData(hPlotsFig5.hPlot1(2,:),0,24,getYLims(hPlotsFig5.hPlot1(2,:)),14);
rescaleData(hPlotsFig5.hPlot2(3),0,50,getYLims(hPlotsFig5.hPlot2(3)),14);
FigName5 = fullfile(folderSave_Figs,['Figure 5_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) ]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotData_histogram(hPlotsFig5,NI_Data_allElecsEvoked,fittingParams)
Criterion = 'none';
plotBars(hPlotsFig6,fileSave1,fileSave3,colorScheme,Criterion)
FigName6 = fullfile(folderSave_Figs,['Figure 6_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff)...
    '_spikesT' num2str(round(1000*timeRangeForComputation_spikes(1))) '_' num2str(round(1000*timeRangeForComputation_spikes(2))) ...
    '_lfpT' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_Criterion_' Criterion '_' gridType '_UnitID' num2str(elecParams.unitID) ]);

saveas(hFigure1,[FigName1 '.fig'])
saveas(hFigure1,[FigName1,'.tif'])
saveas(hFigure2,[FigName2 '.fig'])
saveas(hFigure2,[FigName2,'.tif'])
saveas(hFigure3,[FigName3 '.fig'])
saveas(hFigure3,[FigName3,'.tif'])
saveas(hFigure4,[FigName4 '.fig'])
saveas(hFigure4,[FigName4,'.tif'])
saveas(hFigure5,[FigName5 '.fig'])
saveas(hFigure5,[FigName5,'.tif'])
saveas(hFigure6,[FigName6 '.fig'])
saveas(hFigure6,[FigName6,'.tif'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Accessory Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotData_spikes(hPlot,data,timeRangeForComputation,elecInfo,NormalizeDataFlag,combineUniqueElectrodeData,colorScheme)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;

if NormalizeDataFlag
    data = normalizeData(data);
end

% mean PSTH data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x timeVals
% psthData1 = squeeze(mean(data.data,1));
psthData1 = data.data;

% spike rate data: con_Ori2 (rows) x con_Ori2 (columns)
spikeRateDataST_elecwise = squeeze(data.analysisDataST);
spikeRateDataBL_elecwise = squeeze(data.analysisData_cBL);
dSpikeRateData_elecwise1 = spikeRateDataST_elecwise - spikeRateDataBL_elecwise;


if combineUniqueElectrodeData
    dSpikeRateData_elecwise = combineData(dSpikeRateData_elecwise1,elecInfo);
else
    dSpikeRateData_elecwise = dSpikeRateData_elecwise1;
end

if combineUniqueElectrodeData
    psthData = squeeze(mean(combineTimeSeriesData(psthData1,elecInfo),1));
else
    dSpikeRateData_elecwise = squeeze(mean(psthData1,1));
end


% dSpikeRateData_elecwise = allDataTMP;

m_dSpikeRateData =  squeeze(mean(dSpikeRateData_elecwise,1));
sem_dspikeRate = squeeze(std(dSpikeRateData_elecwise,[],1)./sqrt(size(dSpikeRateData_elecwise,1)));
median_dSpikeRateData = squeeze(median(dSpikeRateData_elecwise,1));

% computing N.I. population
NI_dSpikeData_Cohen =  (dSpikeRateData_elecwise(:,1,1) + dSpikeRateData_elecwise(:,5,5))./squeeze(dSpikeRateData_elecwise(:,1,5));
NI_dSpikeData_Ray =  (2*dSpikeRateData_elecwise(:,1,5)./squeeze((dSpikeRateData_elecwise(:,1,1))+ (dSpikeRateData_elecwise(:,5,5))))-1;
mNI_dSpikeData_Cohen = mean(NI_dSpikeData_Cohen,1);
median_dSpikeData_Cohen = median(NI_dSpikeData_Cohen,1);
mNI_dSpikeData_Ray = mean(NI_dSpikeData_Ray,1);
median_dSpikeData_Ray = median(NI_dSpikeData_Ray,1);
% PSTH plots
if strcmp(colorScheme,'color')
    colors = jet(length(cValsUnique));
elseif strcmp(colorScheme,'grayscale')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end



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
displayRange(hPlot.hPlot1,[timeRangeForComputation(1) timeRangeForComputation(2)],getYLims(hPlot.hPlot1),'k');


for i = 1:length(cValsUnique)
    set(hPlot.hPlot1(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end

% color-coded plots of change in Spike Rates
imagesc(median_dSpikeRateData,'parent',hPlot.hPlot2(1));
colorBar_rlvSpikeRate = colorbar(hPlot.hPlot2(1));
colorYlabelHandle = get(colorBar_rlvSpikeRate,'Ylabel');
set(colorYlabelHandle,'String','Change in Spike Rate (spikes/s)','fontSize',14);
plotPos = get(hPlot.hPlot2(1),'Position');
set(hPlot.hPlot2(1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
% title(hPlot.hPlot2(1),['Median NI_C_o_h_e_n: ',num2str(median_dSpikeData_Cohen) ', Median NI_R_a_y: ',num2str(median_dSpikeData_Ray)],'fontWeight','bold');
set(hPlot.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(1),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(2),'Contrast of Ori 2(%)');

% NI population histogram
edges = [-inf -0.1:0.1:1.6 inf];
edgeCenters = [-0.1 0:0.1:1.6 1.7];
xlim(hPlot.hPlot2(2),[-0.1 3]);
n = histcounts(NI_dSpikeData_Cohen,edges);
barHandle = bar(hPlot.hPlot2(2),edgeCenters,n);
set(hPlot.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hPlot.hPlot2(2),['Median NI: ',num2str(round(median_dSpikeData_Cohen,2))],'fontWeight','bold');
xlabel(hPlot.hPlot2(2),'Normalization Index');ylabel(hPlot.hPlot2(2),'Number of electrodes');
displayRange(hPlot.hPlot2(2),[1 2],getYLims(hPlot.hPlot2(2)),'k')
rescaleData(hPlot.hPlot2(2),-0.5,2.3,getYLims(hPlot.hPlot2(2)),14);



% color scheme of medianChangeinFR and population histogram
if strcmp(colorScheme,'grayscale')
    colormap(hPlot.hPlot2(1),'gray');
    set(barHandle,'facecolor',[0.5 0.5 0.5]);
end

% CRF plots
versionNum = 2;
fittingParams = normalizationModelFit(dSpikeRateData_elecwise,versionNum);
cList = [0 6.25 12.5 25 50]; %[0 1 2 4 8]/16;
mp = squeeze(mean(fittingParams.dataP,1));

parMP = getParametersPlaid(mp,versionNum);
[dp,pmp] = getResponseMatrixPlaid(parMP,mp,versionNum);
expVar = 1 - (dp/sum((mp(:)-mean(mp(:))).^2));

plot(hPlot.hPlot2(3),cList,mp(5,:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
hold (hPlot.hPlot2(3),'on');

plot(hPlot.hPlot2(3),cList,flip(mp(:,1)),'color',[0.8 0.8 0.8],'marker','s','linestyle','none');

hold (hPlot.hPlot2(3),'on');
plot(hPlot.hPlot2(3),cList,[mp(5,1) mp(4,2) mp(3,3) mp(2,4) mp(1,5)],'color','k','marker','v','linestyle','none');
plot(hPlot.hPlot2(3),cList,pmp(5,:),'color',[0.5 0.5 0.5]);
plot(hPlot.hPlot2(3),cList,flip(pmp(:,1),1),'color',[0.8 0.8 0.8]);

plot(hPlot.hPlot2(3),cList,[pmp(5,1) pmp(4,2) pmp(3,3) pmp(2,4) pmp(1,5)],'color','k');
title(hPlot.hPlot2(3),[{['\alpha = ' num2str(parMP(3),3) ', \sigma = ' num2str(parMP(4),3)]} {['ExpVar = ' num2str(round(100*expVar)) '%']}]);
hold(hPlot.hPlot2(3),'off');
text(0.7,0.25,'Grating 1','color',[0.5 0.5 0.5],'fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.7,0.2,'Grating 2','color',[0.8 0.8 0.8],'fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.7,0.15,'Plaid','color','k','fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(3))
set(hPlot.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlot.hPlot2(3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(3),'Change in  Spike Rate (spike/s)');

% fix symmetry of axes boundary
plotPos = get(hPlot.hPlot1(3),'Position');
plotPos2 = get(hPlot.hPlot2(2),'Position');
set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos2(2) plotPos(3) plotPos2(4)]);

plotPos = get(hPlot.hPlot1(5),'Position');
plotPos2 = get(hPlot.hPlot2(3),'Position');
set(hPlot.hPlot2(3),'Position',[plotPos(1) plotPos2(2) plotPos2(3)-(plotPos(1)-plotPos2(1)) plotPos2(4)]);

rescaleData(hPlot.hPlot1,-0.1,0.5,getYLims(hPlot.hPlot1),14);
rescaleData(hPlot.hPlot2(3),0,50,getYLims(hPlot.hPlot2(3)),14);

end

function plotData_energy(hPlot,data,elecInfo,combineUniqueElectrodeData,colorScheme)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;
num_freqRanges = length(data.analysisDataST); % gamma, hi-gamma Induced power is being plotted; alpha and SSVEP power ignored

% mean energy data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x freqVals
% energyVsFrequencyDataST = squeeze(mean(data.dataST(:,1,:,:,:),1));
% energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,1,:,:,:),1));
% dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);

if combineUniqueElectrodeData
    energyVsFrequencyDataST = squeeze(mean(combineTimeSeriesData(data.dataST(:,1,:,:,:),elecInfo),1));
    energyVsFrequencyDataBL = squeeze(mean(combineTimeSeriesData(data.data_cBL(:,1,:,:,:),elecInfo),1));
    dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
else
    energyVsFrequencyDataST = squeeze(mean(data.dataST(:,1,:,:,:),1));
    energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,1,:,:,:),1));
    dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
end







% PSD and deltaPSD plots
if strcmp(colorScheme,'color')
    colors = jet(length(cValsUnique));
elseif strcmp(colorScheme,'grayscale')|| strcmp(colorScheme,'greyscale')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end

cFlipped_Indices = flip(1:length(cValsUnique2));

% plotting baseline PSD and baseline deltaPSD plots
for c_Ori1 = 1:length(cValsUnique)
    plot(hPlot.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
    hold(hPlot.hPlot1(1,c_Ori1),'on');
    plot(hPlot.hPlot1(2,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)-energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
    hold(hPlot.hPlot1(2,c_Ori1),'on');
end

% plotting stimulus PSD and deltaPSD plots
for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlot.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataST(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(1,c_Ori1),'on');
        plot(hPlot.hPlot1(2,c_Ori1),data.freqVals,squeeze(dEnergyVsFrequencyData(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot.hPlot1(2,c_Ori1),'on');
    end
end

% set axes limits and properties for PSD and delta PSD plots
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
for i = 2: num_freqRanges-1
    
    clear energyDataST energyDataBL NI_population_energy Absolute NI_population_energyRelative
    clear dEnergyData_elecwise
    
    % computing energy population data
    energyDataST = data.analysisDataST{i}; % raw Power from MT converted to logarithm values to the base 10
    energyDataBL = data.analysisData_cBL{i}; % raw Power from MT converted to logarithm values to the base 10
    dEnergyData1 = 10*(energyDataST - energyDataBL); % in decibels
    
    if combineUniqueElectrodeData
        dEnergyData = combineData(dEnergyData1,elecInfo);
    else
        dEnergyData = dEnergyData1;
    end
    
    m_dEnergyData = squeeze(mean(dEnergyData,1));
    median_dEnergyData = squeeze(median(dEnergyData,1));
    
    % computing N.I. population data
    NI_dEnergyData_Cohen =  (dEnergyData(:,1,1)+dEnergyData(:,5,5))./squeeze(dEnergyData(:,1,5));
    %     NI_dEnergyData_Ray =  (2*dEnergyData(:,1,5)./squeeze((dEnergyData(:,1,1))+ (dEnergyData(:,5,5))))-1;
    
    mNI_dEnergyData_Cohen = mean(NI_dEnergyData_Cohen,1);
    median_dEnergyData_Cohen = median(NI_dEnergyData_Cohen,1);
    %     mNI_dEnergyData_Ray = mean(NI_dEnergyData_Ray,1);
    %     median_dEnergyData_Ray = median(NI_dEnergyData_Ray,1);
    %
    
    % Color coded Plots of energyData
    imagesc(median_dEnergyData,'parent',hPlot.hPlot2(i-1,1));
    colorBar_rlvPSD = colorbar(hPlot.hPlot2(i-1,1));
    colorYlabelHandle = get(colorBar_rlvPSD,'Ylabel');
    set(colorYlabelHandle,'String',{'\Delta Power (dB)'},'fontSize',12);
    plotPos = get(hPlot.hPlot2(i-1,1),'Position');
    set(hPlot.hPlot2(i-1,1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.028 plotPos(4)]);
    % caxis(hPlot.hPlot2(2),[0 4]);
    set(hPlot.hPlot2(i-1,1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
%     title(hPlot.hPlot2(i-1,1),['Median NI: ',num2str((median_dEnergyData_Ray))],'fontSize',14,'fontWeight','bold');
    
    if i==3
        set(hPlot.hPlot2(i-1,1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
        xlabel(hPlot.hPlot2(i-1,1),'Contrast of Ori 1(%)','fontSize',12);ylabel(hPlot.hPlot2(i-1,1),'Contrast of Ori 2(%)','fontSize',12);
        
    else
        set(hPlot.hPlot2(i-1,1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',[],'YTickLabel',[]);
    end
    
    
    % NI population histogram
    if i == 2 % gamma
        edges = [-inf 1.6:0.2:5.6 inf];
        edgeCenters = [1.5 1.7:0.2:5.6 5.7];
        n = histcounts(NI_dEnergyData_Cohen,edges);
    elseif i == 3 % high-gamma
        edges = [-inf 0:0.2:2.6 inf];
        edgeCenters = [-0.1 0.1:0.2:2.6 2.7];
        n = histcounts(NI_dEnergyData_Cohen,edges);
    end
    barHandle = bar(hPlot.hPlot2(i-1,2),edgeCenters,n);
    if i == 2
        xlim(hPlot.hPlot2(i-1,2),[0 8]);
    elseif i == 3
        xlim(hPlot.hPlot2(i-1,2),[0 3]);
    end
    
    set(hPlot.hPlot2(i-1,2),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    title(hPlot.hPlot2(i-1,2),['Median NI: ',num2str(round(median_dEnergyData_Cohen,2))],'fontSize',12,'fontWeight','bold');
    if i==3
        xlabel(hPlot.hPlot2(i-1,2),'Normalization Index','fontSize',12);ylabel(hPlot.hPlot2(i-1,2),'Number of electrodes','fontSize',12);
    else
    end
    
    xlim()
    displayRange(hPlot.hPlot2(i-1,2),[1 2],getYLims(hPlot.hPlot2(i-1,2)),'k');
    

    % color scheme of medianChangeinFR and population histogram
    if strcmp(colorScheme,'grayscale')
        colormap(hPlot.hPlot2(i-1,1),'gray');
        set(barHandle,'facecolor',[0.5 0.5 0.5]);
    end
    
    % CRF
    %     % Grating
    versionNum = 2;
    fittingParams = normalizationModelFit(dEnergyData,versionNum);
    cList = [0 6.25 12.5 25 50];%[0 1 2 4 8]/16;
    
    mp = squeeze(mean(fittingParams.dataP,1));
    
    %     mp = squeeze(meanP(i,:,:));
    parMP = getParametersPlaid(mp,versionNum);
    [dp,pmp] = getResponseMatrixPlaid(parMP,mp,versionNum);
    expVar = 1 - (dp/sum((mp(:)-mean(mp(:))).^2));
    
    plot(hPlot.hPlot2(i-1,3),cList,mp(5,:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
    hold (hPlot.hPlot2(i-1,3),'on');
    plot(hPlot.hPlot2(i-1,3),cList,flip(mp(:,1)),'color',[0.8 0.8 0.8],'marker','s','linestyle','none');
    plot(hPlot.hPlot2(i-1,3),cList,[mp(5,1) mp(4,2) mp(3,3) mp(2,4) mp(1,5)],'color','k','marker','v','linestyle','none');
    plot(hPlot.hPlot2(i-1,3),cList,pmp(5,:),'color',[0.5 0.5 0.5]);
    plot(hPlot.hPlot2(i-1,3),cList,flip(pmp(:,1)),'color',[0.8 0.8 0.8]);
        
    plot(hPlot.hPlot2(i-1,3),cList,[pmp(5,1) pmp(4,2) pmp(3,3) pmp(2,4) pmp(1,5)],'color','k');
    title(hPlot.hPlot2(i-1,3),['\alpha=' num2str(parMP(3),3) ', \sigma=' num2str(parMP(4),3) ',ExpVar=' num2str(round(100*expVar)) '%']);
    hold(hPlot.hPlot2(i-1,3),'off');
    
    % % CRF
    if i == 3
        text(0.5,0.8,'Grating 1','color',[0.5 0.5 0.5],'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(i-1,3))
        text(0.5,0.7,'Grating 2','color',[0.8 0.8 0.8],'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(i-1,3))
        text(0.5,0.6,'Plaid','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(i-1,3))
    end
    set(hPlot.hPlot2(i-1,3),'fontSize',12,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    if i==3
        set(hPlot.hPlot2(i-1,3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
    else
        set(hPlot.hPlot2(i-1,3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',[]);
    end
    
    if i == 3
        xlabel(hPlot.hPlot2(i-1,3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(i-1,3),'Change in Power(dB)');
    end
    
    plotPos = get(hPlot.hPlot2(i-1,2),'Position');
    plotPos2 = get(hPlot.hPlot2(i-1,3),'Position');
    set(hPlot.hPlot2(i-1,3),'Position',[plotPos2(1)-(plotPos(3)-plotPos2(3)) plotPos2(2) plotPos(3) plotPos2(4)]);
    
end

rescaleData(hPlot.hPlot1,0,250,getYLims(hPlot.hPlot1),12);
rescaleData(hPlot.hPlot1(1,:),0,250,[-1.5 3.5],12);
rescaleData(hPlot.hPlot1(2,:),0,250,[-4 10],12);
rescaleData(hPlot.hPlot2(:,3),-0.5,50,[-0.2 6],12);
rescaleData(hPlot.hPlot2(1,2),-0.5,5.6,getYLims(hPlot.hPlot2(2)),12);

end
function plotData_SSVEP(hPlot,data,elecInfo,combineUniqueElectrodeData,colorScheme)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;
% num_freqRanges = length(data.analysisDataST)-1; % alpha, gamma, hi-gamma Induced power is being plotted

% mean energy data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x
% freqVals
% energyVsFrequencyDataST = squeeze(mean(data.dataST(:,2,:,:,:),1));
% energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,2,:,:,:),1));
% dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
if combineUniqueElectrodeData
    energyVsFrequencyDataST = squeeze(mean(combineTimeSeriesData(data.dataST(:,2,:,:,:),elecInfo),1));
    energyVsFrequencyDataBL = squeeze(mean(combineTimeSeriesData(data.data_cBL(:,2,:,:,:),elecInfo),1));
    dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
else
    energyVsFrequencyDataST = squeeze(mean(data.dataST(:,2,:,:,:),1));
    energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,2,:,:,:),1));
    dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
end


% PSD plots
if strcmp(colorScheme,'color')
    colors = jet(length(cValsUnique));
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end
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

% Seeting axes properties
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

% energy data: con_Ori2 (rows) x con_Ori2 (columns)

clear energyDataST energyDataBL NI_population_energy Absolute NI_population_energyRelative
energyDataST = data.analysisDataST{4};
energyDataBL = data.analysisData_cBL{4};
dEnergyData1 = 10*(energyDataST - energyDataBL); %across elecs
if combineUniqueElectrodeData
    dEnergyData = combineData(dEnergyData1,elecInfo);
else
    dEnergyData = dEnergyData1;
end

m_dEnergyData = squeeze(mean(dEnergyData,1));
median_dEnergyData = squeeze(median(dEnergyData,1));


% computing N.I. population
NI_dEnergyData_Cohen = (dEnergyData(:,1,1) + dEnergyData(:,5,5))./squeeze(dEnergyData(:,1,5));
mNI_dEnergyData_Cohen = mean(NI_dEnergyData_Cohen,1);
median_dEnergyData_Cohen = median(NI_dEnergyData_Cohen,1);
% NI_dEnergyData_Ray =  (2*dEnergyData(:,1,5)./squeeze((dEnergyData(:,1,1))+ (dEnergyData(:,5,5))))-1;
% mNI_dEnergyData_Ray = mean(NI_dEnergyData_Ray,1);
% median_dEnergyData_Ray = median(NI_dEnergyData_Ray,1);


% Color coded Plots of SSVEP energyData
imagesc(median_dEnergyData,'parent',hPlot.hPlot2(1));
colorBar_rlvPSD = colorbar(hPlot.hPlot2(1));
colorYlabelHandle = get(colorBar_rlvPSD,'Ylabel');
set(colorYlabelHandle,'String','Change in Power(dB)','fontSize',14);
plotPos = get(hPlot.hPlot2(1),'Position');
set(hPlot.hPlot2(1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
% title(hPlot.hPlot2(1),['Median NI: ',num2str(median_dEnergyData_Cohen)],'fontWeight','bold');
set(hPlot.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(1),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(1),'Contrast of Ori 2(%)');

% NI population histogram
edges = [-inf 0:0.2:2.6 inf];
edgeCenters = [-0.2 0.2:0.2:2.6 2.8];
n = histcounts(NI_dEnergyData_Cohen,edges);
barHandle = bar(hPlot.hPlot2(2),edgeCenters,n);
set(hPlot.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
title(hPlot.hPlot2(2),['Median NI: ',num2str(round(median_dEnergyData_Cohen,2))],'fontSize',14,'fontWeight','bold');
if i==3
    %         set(hPlot.hPlot2(i-1,2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
    xlabel(hPlot.hPlot2(2),'Normalization Index','fontSize',12);ylabel(hPlot.hPlot2(2),'Number of electrodes','fontSize',12);
    
else
    %         set(hPlot.hPlot2(i-1,2),'XTickLabel',[],'YTickLabel',[]);
end
if strcmp(colorScheme,'grayscale')
    colormap(hPlot.hPlot2(1),'gray');
    set(barHandle,'facecolor',[0.5 0.5 0.5]);
end
xlim(hPlot.hPlot2(2),[-1 4]);
displayRange(hPlot.hPlot2(2),[1 2],getYLims(hPlot.hPlot2(2)),'k');


% CRF
% errorbar(cValsUnique,dEnergyData(end,:),sem_dEnergyData(end,:),...
%     'Marker','o','LineWidth',2,'color',colors(end,:,:),'parent',hPlot.hPlot2(3))
% hold(hPlot.hPlot2(3),'on');
% errorbar(cValsUnique,flip(dEnergyData(:,1)),sem_dEnergyData(:,1),...
%     'Marker','o','LineWidth',2,'color',colors(1,:,:),'parent',hPlot.hPlot2(3))
% errorbar(cValsUnique,diag(flipud(dEnergyData)),diag(flipud(sem_dEnergyData)),'Marker','v','LineStyle','none','LineWidth',2,'color','k','parent',hPlot.hPlot2(3));
% hold(hPlot.hPlot2(3),'on');
% errorbar(cValsUnique,avg_dEnergyData,sem_avg_dEnergyData,'Marker','o','LineStyle','none','LineWidth',2,'color',[0.5 0.5 0.5],'parent',hPlot.hPlot2(3));
% % errorbar(cValsUnique,sum_dEnergyData,sem_sum_dEnergyData,'Marker','o','LineStyle','--','LineWidth',2,'color',[0.8 0.8 0.8],'parent',hPlot.hPlot2(3));
%
versionNum = 2;
fittingParams = normalizationModelFit(dEnergyData,versionNum);
cList = [0 6.25 12.5 25 50];%[0 1 2 4 8]/16;

mp = squeeze(mean(fittingParams.dataP,1));

%     mp = squeeze(meanP(i,:,:));
parMP = getParametersPlaid(mp,versionNum);
[dp,pmp] = getResponseMatrixPlaid(parMP,mp,versionNum);
expVar = 1 - (dp/sum((mp(:)-mean(mp(:))).^2));

plot(hPlot.hPlot2(3),cList,mp(5,:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
hold (hPlot.hPlot2(3),'on');
plot(hPlot.hPlot2(3),cList,flip(mp(:,1)),'color',[0.8 0.8 0.8],'marker','s','linestyle','none');
plot(hPlot.hPlot2(3),cList,[mp(5,1) mp(4,2) mp(3,3) mp(2,4) mp(1,5)],'color','k','marker','v','linestyle','none');
plot(hPlot.hPlot2(3),cList,pmp(5,:),'color',[0.5 0.5 0.5]);
plot(hPlot.hPlot2(3),cList,flip(pmp(:,1)),'color',[0.8 0.8 0.8]);

plot(hPlot.hPlot2(3),cList,[pmp(5,1) pmp(4,2) pmp(3,3) pmp(2,4) pmp(1,5)],'color','k');
title(hPlot.hPlot2(3),[{['\alpha=' num2str(parMP(3),3) ', \sigma=' num2str(parMP(4),3)]} {['ExpVar=' num2str(round(100*expVar)) '%']}]);
% hold(hPlot.hPlot2(3),'off');




% % CRF
% text(0.5,0.05,'cOri 2: 0%','color',colors(end,:,:),'fontWeight','bold','fontSize',8,'unit','normalized','parent',hPlot.hPlot2(3))
% text(0.5,0.1,'cOri 1: 0%','color',colors(1,:,:),'fontWeight','bold','fontSize',8,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.5,0.25,'Grating 1','color',[0.5 0.5 0.5],'fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.5,0.2,'Grating 2','color',[0.8 0.8 0.8],'fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.5,0.15,'Plaid','color','k','fontWeight','bold','fontSize',10,'unit','normalized','parent',hPlot.hPlot2(3))
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
hold(hPlot.hPlot2(3),'off');

end

function plotExampleElectrodeData(hPlotsFig1,hPlotsFig2,elecNum,elecInfo,firingRateData,energyData,energyDataTF,colorScheme)
cValsUnique = [0 6.25 12.5 25 50];
cValsUnique2 = cValsUnique;
if strcmp(colorScheme,'color')
    colors = jet(length(cValsUnique));
elseif strcmp(colorScheme,'grayscale')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end
cFlipped_Indices = flip(1:length(cValsUnique2));
psthData = squeeze(firingRateData.data(elecNum,:,:,:,:));
spikeRasterData = firingRateData.spikeRasterData(elecNum,1,:,:);

tickLengthPlot = 2*get(hPlotsFig1.hPlot1(1),'TickLength');

combinedData = combineData(firingRateData,energyData,energyDataTF,elecInfo,elecNum);

psthData2 = psthData;
spikeRasterData2 = spikeRasterData;
if elecNum == 29 || elecNum == 40 || elecNum == 169
    psthData2 = flip(flip(permute(psthData,[2 1 3]),1),2);
    spikeRasterData2 = flip(flip(permute(firingRateData.spikeRasterData(elecNum,1,:,:),[1 2 4 3]),4),3);
    combinedData{1} = flip(flip(permute(combinedData{1},[2 1 3]),1),2);
    combinedData{2} = flip(flip(permute(combinedData{2},[2 1 3 4]),1),2);
    combinedData{3} = flip(flip(permute(combinedData{3},[2 1 3]),1),2);
    combinedData{4} = flip(flip(permute(combinedData{4},[2 1 3]),1),2);

end
for c_Ori2 = 1:5
    for c_Ori1 = 1:5
        X = spikeRasterData2{1,1,c_Ori2,c_Ori1};
        axes(hPlotsFig1.hPlot1(c_Ori2,c_Ori1)); %#ok<LAXES>
        hold(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),'on');
        rasterplot(X,1:length(X),'k',2);
        plot(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),firingRateData.timeVals,squeeze(psthData2(c_Ori2,c_Ori1 ,:)),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
        set(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')

    end
end
xlabel(hPlotsFig1.hPlot1(5,1),'Time (s)'); ylabel(hPlotsFig1.hPlot1(5,1),[{'Firing Rate'} {'(spikes/s)'}]);
rescaleData(hPlotsFig1.hPlot1,-0.1,0.5,[-5 0]+ getYLims(hPlotsFig1.hPlot1),14);


for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlotsFig1.hPlot2(1,c_Ori1),firingRateData.timeVals,squeeze(psthData2(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig1.hPlot2(1,c_Ori1),'on');
        plot(hPlotsFig1.hPlot2(2,c_Ori1),firingRateData.timeVals,squeeze(combinedData{1}(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig1.hPlot2(2,c_Ori1),'on');
    end
end

for i = 1:length(cValsUnique)
    set(hPlotsFig1.hPlot2(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
xlabel(hPlotsFig1.hPlot2(1,1),'Time (s)'); ylabel(hPlotsFig1.hPlot2(1,1),[{'Firing Rate'} {'(spikes/s)'}]);
rescaleData(hPlotsFig1.hPlot2,-0.1,0.5,[-5 0]+getYLims(hPlotsFig1.hPlot2),14);


% Figure 2


% dEnergyTF = 10*(energyDataTF.data - energyDataTF.data_cBL);
% energyVsFrequencyDataST = squeeze(mean(energyData.dataST(elecNum,1,:,:,:),1));
% energyVsFrequencyDataBL = squeeze(mean(energyData.data_cBL(elecNum,1,:,:,:),1));
% dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
dEnergyTF = 10*combinedData{2};
% energyVsFrequencyDataST = squeeze(mean(energyData.dataST(elecNum,1,:,:,:),1));
% energyVsFrequencyDataBL = squeeze(mean(energyData.data_cBL(elecNum,1,:,:,:),1));
energyVsFrequencyDataBL = combinedData{4};
dEnergyVsFrequencyData = 10*combinedData{3};



dEnergyTF2 = dEnergyTF;
dEnergyVsFrequencyData2 = dEnergyVsFrequencyData;
% if elecNum == 29 || elecNum == 40 || elecNum == 169
%     dEnergyTF2 = flip(flip(permute(dEnergyTF,[1 2 4 3 5 6]),4),3);
%     dEnergyVsFrequencyData2 = flip(flip(permute(dEnergyVsFrequencyData,[2 1 3]),2),1);
% end

for c_Ori2 = 1:5
    for c_Ori1 = 1:5
        pcolor(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),energyDataTF.timeVals,energyDataTF.freqVals,squeeze(dEnergyTF2(c_Ori2,c_Ori1 ,:,:)));
        caxis(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),[-2 10]);
        shading(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'interp'); %caxis([-5 10]);
        if strcmp(colorScheme,'color')
            colormap(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'jet');
        elseif strcmp(colorScheme,'greyscale')||strcmp(colorScheme,'grayscale')
            colormap(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),flipud(gray));
        end
        set(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')

        if c_Ori2 == 5 && c_Ori1 == 5
            colorBar_tf = colorbar(hPlotsFig2.hPlot1(c_Ori2,c_Ori1));
            colorBar_YLabelHandle = get(colorBar_tf,'Ylabel');
            set(colorBar_YLabelHandle,'String',{'\Delta Power (dB)'},'fontsize',14);
            plotPos = get(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'Position');
            set(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.042 plotPos(4)]);
        end
    end
end
xlabel(hPlotsFig2.hPlot1(5,1),'Time (s)'); ylabel(hPlotsFig2.hPlot1(5,1),'Frequency (Hz)');

rescaleData(hPlotsFig2.hPlot1,-0.1,0.5,getYLims(hPlotsFig2.hPlot1),14);


for c_Ori1 = 1:length(cValsUnique)
    plot(hPlotsFig2.hPlot2(1,c_Ori1),energyData.freqVals,squeeze(energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)-energyVsFrequencyDataBL(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
    hold(hPlotsFig2.hPlot2(1,c_Ori1),'on');
end

for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlotsFig2.hPlot2(1,c_Ori1),energyData.freqVals,squeeze(dEnergyVsFrequencyData2(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig2.hPlot2(1,c_Ori1),'on');
    end
end
for i = 1:length(cValsUnique)
    set(hPlotsFig2.hPlot2(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
xlabel(hPlotsFig2.hPlot2(1,1),'Frequency (Hz)'); ylabel(hPlotsFig2.hPlot2(1,1),{'Change in','Power (dB)'});
rescaleData(hPlotsFig2.hPlot2,0,250,getYLims(hPlotsFig2.hPlot2),14);

end

% combineData
% function xAll = combineData(firingRateData,energyData,energyDataTF,elecInfo,elecNum)
% combineUniqueElectrodes=1;
% N=15; snr=2; d=0.75;
% 
% %%%%%%%%%%%%%%%% Subselect electrodes based on N, snr and d %%%%%%%%%%%%%%%
% goodN = (max([elecInfo.N(1,:);elecInfo.N(2,:)])>N);
% goodSNR = (elecInfo.SNR>snr);
% goodD = (elecInfo.d <= d);
% goodElectrodeIDs = (goodN & goodSNR & goodD);
% 
% %%%%%%%%%%%%%%%%%%% Get information about the electrodes %%%%%%%%%%%%%%%%%%
% allElectrodes = [];
% for i=1:length(elecInfo.elecs)
%     x = elecInfo.elecs{i}{end};
%     if i<=13
%         allElectrodes = cat(2,allElectrodes,x);
%     else
%         allElectrodes = cat(2,allElectrodes,x+100); % We add 100 to each electrode ID in Monkey 2 so that electrode numbers are different from monkey 1
%     end
% end
% 
% %%%%%% For each electrode, find the IDs that need to be averaged %%%%%%%%%%
% clear goodIDList
% if combineUniqueElectrodes
%     
%     uniqueGoodElectrodeList = unique(allElectrodes(goodElectrodeIDs));
%     numGoodElectrodes = length(uniqueGoodElectrodeList);
%     
%     goodIDList = cell(1,numGoodElectrodes);
%     for i=1:numGoodElectrodes
%         goodIDList{i} = find(allElectrodes==uniqueGoodElectrodeList(i));
%     end
% else
%     goodIDs = find(goodElectrodeIDs==1); %#ok<*UNRCH>
%     numGoodElectrodes = length(goodIDs);
%     
%     goodIDList = cell(1,numGoodElectrodes);
%     for i=1:numGoodElectrodes
%         goodIDList{i} = goodIDs(i);
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%% Combine data if necessary %%%%%%%%%%%%%%%%%%%%%%%%%
% % allData = cell(1,4);
% for i= 1:size(goodIDList,2)
%     a = find(elecNum==goodIDList{i}); %#ok<*EFIND>
%     if ~isempty(a)
%         UniqueElecID = i;
%     end
% end
% 
% clear allDataTMP
% % allDataTMP = zeros(numGoodElectrodes,5,5);
% 
% x = squeeze(firingRateData.data(goodIDList{UniqueElecID},1,:,:,:));
% y = squeeze(energyDataTF.data(goodIDList{UniqueElecID},1,:,:,:,:))-squeeze(energyDataTF.data_cBL(goodIDList{UniqueElecID},1,:,:,:,:));
% z = squeeze(energyData.dataST(goodIDList{UniqueElecID},1,:,:,:))-squeeze(energyData.data_cBL(goodIDList{UniqueElecID},1,:,:,:));
% 
% b = squeeze(energyData.data_cBL(goodIDList{UniqueElecID},1,:,:,:));
% 
% if length(goodIDList{UniqueElecID})==1
%     allDataTMP = x;
%     alldTFDataTMP = y;
%     alldPSDDataTMP = z;
%     allblPSDDataTMP = b;
% else
%     xs = zeros(size(x));
%     ys = zeros(size(y));
%     zs = zeros(size(z));
%     bs = zeros(size(b));
%     for j=1:length(goodIDList{UniqueElecID})
%         xs(j,:,:,:) = squeeze(x(j,:,:,:));
%         ys(j,:,:,:,:) = squeeze(y(j,:,:,:,:));
%         zs(j,:,:,:) = squeeze(z(j,:,:,:));
%         bs(j,:,:,:) = squeeze(b(j,:,:,:));
%     end
%     allDataTMP = squeeze(mean(xs,1));
%     alldTFDataTMP = squeeze(mean(ys,1));
%     alldPSDDataTMP = squeeze(mean(zs,1));
%     allblPSDDataTMP = squeeze(mean(bs,1));
% end
% 
% xAll{1} = allDataTMP; xAll{2} = alldTFDataTMP; xAll{3} = alldPSDDataTMP; xAll{4} = allblPSDDataTMP;
% type{1} = 'PSTH';type{2} = 'dTF'; type{3} = 'dPSD';type{3} = 'blPSD'; %#ok<*NASGU>
% 
% % for i=2:4
% %     eDataTMP = squeeze(energyData.dataST{i}) - squeeze(energyData.data_cBL{i});
% %     
% %     allDataTMP = zeros(numGoodElectrodes,5,5);
% %     for j=1:numGoodElectrodes
% %         clear x
% %         x = eDataTMP(goodIDList{j},:,:);
% %         
% %         xs = zeros(size(x));
% %         for k=1:length(goodIDList{j})
% % %             y = squeeze(x(k,:,:));
% % %             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
% %             xs(k,:,:) = squeeze(x(k,:,:));
% %         end
% % 
% %         allDataTMP(j,:,:) = squeeze(mean(xs,1));
% %     end
% %     xAll{i} = allDataTMP;
% % end
% end


% function plotExampleElectrodeData(hPlotsFig6,hPlotsFig7,elecNum,firingRateData,energyDataTF)
%
% psthData = squeeze(firingRateData.data(elecNum,:,:,:,:));
% for c_Ori2 = 1:5
%     for c_Ori1 = 1:5
%         plot(hPlotsFig6.hPlot1(c_Ori2,c_Ori1),firingRateData.timeVals,squeeze(psthData(c_Ori2,c_Ori1 ,:)),'b')
%     end
% end
%
%
% for c_Ori2 = 1:5
%     for c_Ori1 = 1:5
%         pcolor(hPlotsFig7.hPlot1(c_Ori2,c_Ori1),energyDataTF.timeVals,energyDataTF.freqVals,10*squeeze(mean(energyDataTF.BLsubtractedData(elecNum,1,c_Ori2,c_Ori1 ,:,:),1)));
%         shading(hPlotsFig7.hPlot1(c_Ori2,c_Ori1),'interp');
%         colormap(hPlotsFig7.hPlot1(c_Ori2,c_Ori1),'jet')
%         if c_Ori2 == 5 && c_Ori1 == 5
%             colorBar_tf = colorbar(hPlotsFig7.hPlot1(c_Ori2,c_Ori1));
%             colorBar_YLabelHandle = get(colorBar_tf,'Ylabel');
%             set(colorBar_YLabelHandle,'String',{'\Delta Power (dB)'},'fontsize',14);
%             plotPos = get(hPlotsFig7.hPlot1(c_Ori2,c_Ori1),'Position');
%             set(hPlotsFig7.hPlot1(c_Ori2,c_Ori1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.042 plotPos(4)]);
%         end
%     end
% end
%
% end

function data2 = combineTimeSeriesData(data,elecInfo)
% combine Spike Data for unique electrodes
combineUniqueElectrodes=1;
N=15; snr=2; d=0.75;

%%%%%%%%%%%%%%%% Subselect electrodes based on N, snr and d %%%%%%%%%%%%%%%
goodN = (max([elecInfo.N(1,:);elecInfo.N(2,:)])>N);
goodSNR = (elecInfo.SNR>snr);
goodD = (elecInfo.d <= d);
goodElectrodeIDs = (goodN & goodSNR & goodD);

%%%%%%%%%%%%%%%%%%% Get information about the electrodes %%%%%%%%%%%%%%%%%%
allElectrodes = [];
for i=1:length(elecInfo.elecs)
    x = elecInfo.elecs{i}{end};
    if i<=13
        allElectrodes = cat(2,allElectrodes,x);
    else
        allElectrodes = cat(2,allElectrodes,x+100); % We add 100 to each electrode ID in Monkey 2 so that electrode numbers are different from monkey 1
    end
end

%%%%%% For each electrode, find the IDs that need to be averaged %%%%%%%%%%
clear goodIDList
if combineUniqueElectrodes
    
    uniqueGoodElectrodeList = unique(allElectrodes(goodElectrodeIDs));
    numGoodElectrodes = length(uniqueGoodElectrodeList);
    
    goodIDList = cell(1,numGoodElectrodes);
    for i=1:numGoodElectrodes
        goodIDList{i} = find(allElectrodes==uniqueGoodElectrodeList(i));
    end
else
    goodIDs = find(goodElectrodeIDs==1); %#ok<*UNRCH>
    numGoodElectrodes = length(goodIDs);
    
    goodIDList = cell(1,numGoodElectrodes);
    for i=1:numGoodElectrodes
        goodIDList{i} = goodIDs(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%% Combine data if necessary %%%%%%%%%%%%%%%%%%%%%%%%%
% allData = cell(1,4);
% clear allDataTMP
% allDataTMP = zeros(numGoodElectrodes,5,5);
for i=1:numGoodElectrodes
    clear x
    x = squeeze(data(goodIDList{i},1,:,:,:));%squeeze(firingRateData.analysisDataST(goodIDList{i},1,:,:) - firingRateData.analysisData_cBL(goodIDList{i},1,:,:));
    
    if length(goodIDList{i})==1
%         y = x;
%         z = flip(y,1); z = (z+z')/2; z = flip(z,1);
        allDataTMP(i,:,:,:) = x;
    else
        xs = zeros(size(x));
        for j=1:length(goodIDList{i})
%             y = squeeze(x(j,:,:));
%             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
            xs(j,:,:,:) = squeeze(x(j,:,:,:));
        end
        allDataTMP(i,:,:,:) = squeeze(mean(xs,1));
    end
end
data2 = allDataTMP;
end

function data2 = combineData(data,elecInfo)
% combine Spike Data for unique electrodes
combineUniqueElectrodes=1;
N=15; snr=2; d=0.75;

%%%%%%%%%%%%%%%% Subselect electrodes based on N, snr and d %%%%%%%%%%%%%%%
goodN = (max([elecInfo.N(1,:);elecInfo.N(2,:)])>N);
goodSNR = (elecInfo.SNR>snr);
goodD = (elecInfo.d <= d);
goodElectrodeIDs = (goodN & goodSNR & goodD);

%%%%%%%%%%%%%%%%%%% Get information about the electrodes %%%%%%%%%%%%%%%%%%
allElectrodes = [];
for i=1:length(elecInfo.elecs)
    x = elecInfo.elecs{i}{end};
    if i<=13
        allElectrodes = cat(2,allElectrodes,x);
    else
        allElectrodes = cat(2,allElectrodes,x+100); % We add 100 to each electrode ID in Monkey 2 so that electrode numbers are different from monkey 1
    end
end

%%%%%% For each electrode, find the IDs that need to be averaged %%%%%%%%%%
clear goodIDList
if combineUniqueElectrodes
    
    uniqueGoodElectrodeList = unique(allElectrodes(goodElectrodeIDs));
    numGoodElectrodes = length(uniqueGoodElectrodeList);
    
    goodIDList = cell(1,numGoodElectrodes);
    for i=1:numGoodElectrodes
        goodIDList{i} = find(allElectrodes==uniqueGoodElectrodeList(i));
    end
else
    goodIDs = find(goodElectrodeIDs==1); %#ok<*UNRCH>
    numGoodElectrodes = length(goodIDs);
    
    goodIDList = cell(1,numGoodElectrodes);
    for i=1:numGoodElectrodes
        goodIDList{i} = goodIDs(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%% Combine data if necessary %%%%%%%%%%%%%%%%%%%%%%%%%
allDataTMP = zeros(numGoodElectrodes,5,5);
for i=1:numGoodElectrodes
    clear x
    x = squeeze(data(goodIDList{i},:,:));%squeeze(firingRateData.analysisDataST(goodIDList{i},1,:,:) - firingRateData.analysisData_cBL(goodIDList{i},1,:,:));
    
    if length(goodIDList{i})==1
%         y = x;
%         z = flip(y,1); z = (z+z')/2; z = flip(z,1);
        allDataTMP(i,:,:) = x;
    else
        xs = zeros(size(x));
        for j=1:length(goodIDList{i})
%             y = squeeze(x(j,:,:));
%             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
            xs(j,:,:) = squeeze(x(j,:,:));
        end
        allDataTMP(i,:,:) = squeeze(mean(xs,1));
    end
end
data2 = allDataTMP;
end
function fittingParams = normalizationModelFit(data,versionNum)

numElectrodes = size(data,1);

for j=1:numElectrodes
%     % Grating
%     g1 = squeeze(data(j,5,:))';
%     g2 = flip(squeeze(data(j,:,1)),2);
%     
%     g = (g1+g2)/2; % Make symmetric
%     parG = getParametersGrating(g);
%     [dg,pg] = getResponseMatrixGrating(parG,g);
%     eG(j) = 1 - (dg/sum((g-mean(g)).^2)); %#ok<*NASGU>
%     sG(j) = parG(2);
%     dataG(j,:) = g;
    
    %Plaid
    z = squeeze(data(j,:,:));
%     z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
    parP = getParametersPlaid(z,versionNum);
    [dp,pz] = getResponseMatrixPlaid(parP,z,versionNum); %#ok<*ASGLU>
    eP(j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
    aP(j) = parP(2);
    sP(j) = parP(3);
    dataP(j,:,:) = z;
end

fittingParams.eP = eP;
fittingParams.aP = aP;
fittingParams.sP = sP;
fittingParams.dataP = dataP;

end

function fittingParams = normalizationModelFitV2(data)

numElectrodes = size(data,1);

for j=1:numElectrodes
%     % Grating
%     g1 = squeeze(data(j,5,:))';
%     g2 = flip(squeeze(data(j,:,1)),2);
%     
%     g = (g1+g2)/2; % Make symmetric
%     parG = getParametersGrating(g);
%     [dg,pg] = getResponseMatrixGrating(parG,g);
%     eG(j) = 1 - (dg/sum((g-mean(g)).^2)); %#ok<*NASGU>
%     sG(j) = parG(2);
%     dataG(j,:) = g;
    
    %Plaid
    z = squeeze(data(j,:,:));
%     z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
    parP = getParametersPlaid(z);
    [dp,pz] = getResponseMatrixPlaid(parP,z); %#ok<*ASGLU>
    eP(j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
    aP(j) = parP(2);
    sP(j) = parP(3);
    dataP(j,:,:) = z;
end

fittingParams.eP = eP;
fittingParams.aP = aP;
fittingParams.sP = sP;
fittingParams.dataP = dataP;

end

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
function [ElectrodeArrayListAll,allGoodElectrodesStrArray,ElectrodeStringListAll,elecInfo,numElecs] = getElectrodesList(fileNameStringTMP,elecParams,timeRangeForComputation,folderSourceString)

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
allGoodNs = [];
allGoodSNRs = [];
alld_elecs = [];

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
    [tmpElectrodeStringList{i},tmpElectrodeStringArrayList{i},tmpElectrodeArrayList{i},goodNs,goodSNRs,d_elecs,goodElectrodes] = getGoodElectrodesSingleSession(monkeyName,expDate,protocolName,gridType,elecParams,timeRangeForComputation,folderSourceString,versionNum);
    allGoodNs = cat(2,allGoodNs,goodNs);
    allGoodSNRs = cat(2,allGoodSNRs,goodSNRs);
    alld_elecs = cat(2,alld_elecs,d_elecs);
    numElecs = numElecs+length(goodElectrodes);
end
allGoodElectrodesStrArray = tmpElectrodeStringArrayList;
ElectrodeStringListAll = tmpElectrodeStringList;
ElectrodeArrayListAll = tmpElectrodeArrayList;

elecInfo.elecs = ElectrodeArrayListAll;
elecInfo.elecIDs = allGoodElectrodesStrArray;
elecInfo.N = allGoodNs;
elecInfo.SNR = allGoodSNRs;
elecInfo.d = alld_elecs; 

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
        hold(plotHandles(i,j),'on');
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

function plotFigures1and2(monkeyName,folderSourceString,timeRangeForComputation_spikes,timeRangeForComputation_lfp,colorScheme)
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
elecParams.snrCutoff = 2.5;
elecParams.dRange = [0 0.75];
elecParams.unitID = 0;
% elecParams.NICutOff = 1;


tapers_MT = [1 1]; % parameters for MT analysis
removeERPFlag = 0;

% Fixed parameters
% timeRangeForComputationBL_spikes = -0.05+[-diff(timeRangeForComputation_spikes) 0];
% timeRangeForComputationBL_lfp = -0.05+[-diff(timeRangeForComputation_lfp) 0];

folderSourceString_Project = strtok(folderSourceString,'\');

folderSave = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures\Figures1and2');

if ~exist(folderSave,'dir')
    mkdir(folderSave)
end
gridType = 'Microelectrode';
combineUniqueElectrodeData = 0;


freqRanges{1} = [8 12]; % alpha
freqRanges{2} = [30 80]; % gamma
freqRanges{3} = [104 250]; % hi-gamma
freqRanges{4} = [16 16];  % SSVEP



%%%%%%%%%%%%%%%%%%%%%%%%%%%% display properties %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 6: Example Single electrode PSTH data for a single session:
hFigure1 = figure(1);
set(hFigure1,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig1.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.65],0.01,0.01,0); linkaxes(hPlotsFig1.hPlot1);
hPlotsFig1.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.13],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

textH{1} = getPlotHandles(1,1,[0.12 0.97 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.12 0.23 0.01 0.01]);


textString = {'A','B'};
for i = 1:2
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end

% % Figure 7: Example Single electrode TF data for a single session:
hFigure2 = figure(2);
set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig2.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.65],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot1);
hPlotsFig2.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.13],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

textH{1} = getPlotHandles(1,1,[0.12 0.97 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.12 0.23 0.01 0.01]);


textString = {'A','B'};
for i = 1:2
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
% %%%%%%%%%%%%%%%%%%%%% Get Session Details for Monkey(s) %%%%%%%%%%%%%%%%%%%
fileNameStringListAll = getFileNameStringList(monkeyName,gridType);
%%%%%%%%%%%%%%%%%%%%% Get data for all Good Electrodes %%%%%%%%%%%%%%%%%%%%
disp('all Electrodes:')
elecParams.oriSelectiveFlag = 0;
elecParams.getSpikeElectrodesFlag = 1;
[electrodeList_All,numElecs] = getElectrodesList(fileNameStringListAll,elecParams,timeRangeForComputation_spikes,folderSourceString);
disp([num2str(numElecs) ' Good Electrodes'])

fileSave1 = fullfile(folderSave,[monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) ...
    '_allElecs_spikesT' num2str(round(1000*timeRangeForComputation_spikes(1))) '_' num2str(round(1000*timeRangeForComputation_spikes(2))) ...
    '_lfpT' num2str(round(1000*timeRangeForComputation_lfp(1))) '_' num2str(round(1000*timeRangeForComputation_lfp(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) '.mat']);


timeRangeParameters(1).stRange = timeRangeForComputation_spikes;
timeRangeParameters(1).blRange = -0.05+[-diff(timeRangeForComputation_spikes) 0];
timeRangeParameters(1).erpRange = [0.05 0.2];

timeRangeParameters(2).stRange = timeRangeForComputation_lfp;
timeRangeParameters(2).blRange = -0.05+[-diff(timeRangeForComputation_lfp) 0];
timeRangeParameters(2).erpRange = [0.05 0.2];

if exist(fileSave1,'file')
    disp(['Loading file ' fileSave1]);
    load(fileSave1);
else
    % get Data all Session for monkey(s) for all Electrodes
    %     timeRangeParameters.blRange = timeRangeForComputationBL;
    [erpData,firingRateData,fftData,energyData,energyDataTF,oriTuningData,NI_Data,~]  = ...
        getDataFigures1and2(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,elecParams,removeERPFlag); %#ok<ASGLU>
    save(fileSave1,'erpData','firingRateData','fftData','energyData','energyDataTF','oriTuningData','NI_Data')
    
end

elecNum = 23;
plotExampleElectrodeData(hPlotsFig1,hPlotsFig2,elecNum,firingRateData,energyData,energyDataTF,colorScheme)
if strcmp(colorScheme,'color')
    folderSaveFigs = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\Figures\Color');
elseif strcmp(colorScheme,'greyscale')||strcmp(colorScheme,'grayscale')
    folderSaveFigs = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\Figures\GrayScale');
end

FigName1 = fullfile(folderSaveFigs,['Figure 1_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_elec' num2str(elecNum)...
    '_T' num2str(round(1000*timeRangeParameters(1).stRange(1))) '_' num2str(round(1000*timeRangeParameters(1).stRange(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) ]);
saveas(hFigure1,[FigName1 '.fig'])
saveas(hFigure1,[FigName1,'.tif'])

% Figure 7

FigName2 = fullfile(folderSaveFigs,['Figure 2_' monkeyName '_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_elec' num2str(elecNum)...
    '_T' num2str(round(1000*timeRangeParameters(2).stRange(1))) '_' num2str(round(1000*timeRangeParameters(2).stRange(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) ]);
saveas(hFigure2,[FigName2 '.fig'])
saveas(hFigure2,[FigName2,'.tif'])
end

function plotExampleElectrodeData(hPlotsFig1,hPlotsFig2,elecNum,firingRateData,energyData,energyDataTF,colorScheme)
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

if elecNum == 9
    psthData2 = psthData;
    spikeRasterData2 = spikeRasterData;
elseif elecNum == 17 || elecNum == 22 || elecNum == 23 || elecNum == 51 || elecNum == 78
    psthData2 = flip(flip(permute(psthData,[2 1 3]),1),2);
    spikeRasterData2 = flip(flip(permute(firingRateData.spikeRasterData(elecNum,1,:,:),[1 2 4 3]),4),3);
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
    end
end

for i = 1:length(cValsUnique)
    set(hPlotsFig1.hPlot2(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
xlabel(hPlotsFig1.hPlot2(1,1),'Time (s)'); ylabel(hPlotsFig1.hPlot2(1,1),[{'Firing Rate'} {'(spikes/s)'}]);
rescaleData(hPlotsFig1.hPlot2,-0.1,0.5,[-5 0]+getYLims(hPlotsFig1.hPlot2),14);

% Figure 2


dEnergyTF = 10*(energyDataTF.data - energyDataTF.data_cBL);
energyVsFrequencyDataST = squeeze(mean(energyData.dataST(elecNum,1,:,:,:),1));
energyVsFrequencyDataBL = squeeze(mean(energyData.data_cBL(elecNum,1,:,:,:),1));
dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);

if elecNum == 9
    dEnergyTF2 = dEnergyTF;
    dEnergyVsFrequencyData2 = dEnergyVsFrequencyData;
elseif elecNum == 17 || elecNum == 22 || elecNum == 23 || elecNum == 51 || elecNum == 78
    dEnergyTF2 = flip(flip(permute(dEnergyTF,[1 2 4 3 5 6]),4),3);
    dEnergyVsFrequencyData2 = flip(flip(permute(dEnergyVsFrequencyData,[2 1 3]),2),1);
end

for c_Ori2 = 1:5
    for c_Ori1 = 1:5
        pcolor(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),energyDataTF.timeVals,energyDataTF.freqVals,squeeze(mean(dEnergyTF2(elecNum,1,c_Ori2,c_Ori1 ,:,:),1)));
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
function [ElectrodeArrayListAll,numElecs] = getElectrodesList(fileNameStringTMP,elecParams,timeRangeForComputation,folderSourceString)

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
    [tmpElectrodeStringList{i},tmpElectrodeArrayList{i},goodElectrodes] = getGoodElectrodesSingleSession(monkeyName,expDate,protocolName,gridType,elecParams,timeRangeForComputation,folderSourceString,versionNum);
    numElecs = numElecs+length(goodElectrodes);
end

ElectrodeArrayListAll = tmpElectrodeArrayList;
end

% Draw lines for timeTange or FreqRange
% function displayRange(plotHandles,range,yLims,colorName)
% [nX,nY] = size(plotHandles);
% %yLims = getYLims(plotHandles);
% 
% yVals = yLims(1):(yLims(2)-yLims(1))/100:yLims(2);
% xVals1 = range(1) + zeros(1,length(yVals));
% xVals2 = range(2) + zeros(1,length(yVals));
% 
% for i=1:nX
%     for j=1:nY
%         hold(plotHandles(i,j),'on');
%         plot(plotHandles(i,j),xVals1,yVals,'color',colorName);
%         plot(plotHandles(i,j),xVals2,yVals,'color',colorName);
%     end
% end
% end

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
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
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
% set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
% set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
% set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
% set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
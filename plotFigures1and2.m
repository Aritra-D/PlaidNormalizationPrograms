function plotFigures1and2(monkeyName,folderSourceString,elecNum,timeRangeForComputation_spikes,timeRangeForComputation_lfp,colorScheme)
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
freqRanges{2} = [32 80]; % gamma
freqRanges{3} = [104 248]; % hi-gamma
freqRanges{4} = [16 16];  % SSVEP



%%%%%%%%%%%%%%%%%%%%%%%%%%%% display properties %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 6: Example Single electrode PSTH data for a single session:
hFigure1 = figure(1);
set(hFigure1,'units','normalized','outerposition',[0 0 1 1])
% hPlotsFig1.hPlot1 = getPlotHandles(5,5,[0.2 0.38 0.55 0.55],0.01,0.01,0); linkaxes(hPlotsFig1.hPlot1);
% hPlotsFig1.hPlot2 = getPlotHandles(2,5,[0.2 0.08 0.55 0.22],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

hPlotsFig1.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.6],0.01,0.01,0); linkaxes(hPlotsFig1.hPlot1);
hPlotsFig1.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.1125],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

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
hPlotsFig2.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.6],0.01,0.01,0); linkaxes(hPlotsFig2.hPlot1);
hPlotsFig2.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.1125],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot2);

% hPlotsFig2.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.65],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot1);
% hPlotsFig2.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.13],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

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
% [electrodeList_All,numElecs] = getElectrodesList(fileNameStringListAll,elecParams,timeRangeForComputation_spikes,folderSourceString);
[electrodeList_All,elecIDs,electrodeListString_All,NList_All,SNRList_All,dList_All,numElecs] = getElectrodesList(fileNameStringListAll,elecParams,timeRangeForComputation_spikes,folderSourceString); %#ok<*ASGLU>
elecInfo.elecs = electrodeList_All;
elecInfo.elecIDs = elecIDs;
elecInfo.N = NList_All;
elecInfo.SNR = SNRList_All;
elecInfo.d = dList_All;  %#ok<*STRNU>

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
        getDataFigures1and2(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,elecParams,removeERPFlag);
    save(fileSave1,'elecInfo','erpData','firingRateData','fftData','energyData','energyDataTF','oriTuningData','NI_Data')
    
end

% elecNum = 29;


plotExampleElectrodeData(hPlotsFig1,hPlotsFig2,elecNum,elecInfo,firingRateData,energyData,energyDataTF,colorScheme)
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

print(hFigure1,[FigName1,'HQ.tif'],'-dtiff','-r300')
print(hFigure2,[FigName2,'HQ.tif'],'-dtiff','-r300')
end

function plotExampleElectrodeData(hPlotsFig1,hPlotsFig2,elecNum,elecInfo,firingRateData,energyData,energyDataTF,colorScheme)

freqRanges{1} = [8 12]; % alpha
freqRanges{2} = [32 80]; % gamma
freqRanges{3} = [104 248]; % hi-gamma
freqRanges{4} = [16 16];  % SSVEP

cValsUnique = [0 6.25 12.5 25 50];
cValsUnique2 = cValsUnique;
if strcmp(colorScheme,'color')
    colors = jet(length(cValsUnique));
elseif strcmp(colorScheme,'grayscale')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end
cFlipped_Indices = flip(1:length(cValsUnique2));

tickLengthPlot = 2*get(hPlotsFig1.hPlot1(1),'TickLength');

combinedData = combineData(firingRateData,energyData,energyDataTF,elecInfo,elecNum);

if ~isempty(elecNum)
    psthData = squeeze(firingRateData.data(elecNum,:,:,:,:));
    spikeRasterData = firingRateData.spikeRasterData(elecNum,1,:,:);
    if elecNum == 29 || elecNum == 40 || elecNum == 169
        psthData2 = flip(flip(permute(psthData,[2 1 3]),1),2);
        spikeRasterData2 = flip(flip(permute(firingRateData.spikeRasterData(elecNum,1,:,:),[1 2 4 3]),4),3);
        combinedData{1} = flip(flip(permute(combinedData{1},[2 1 3]),1),2);
        combinedData{2} = flip(flip(permute(combinedData{2},[2 1 3 4]),1),2);
        combinedData{3} = flip(flip(permute(combinedData{3},[2 1 3]),1),2);
        combinedData{4} = flip(flip(permute(combinedData{4},[2 1 3]),1),2);
        
    end
    
    rasterScale = round(max(psthData2(:)))+1;
    for c_Ori2 = 1:5
        for c_Ori1 = 1:5
            clear X
            X = spikeRasterData2{1,1,c_Ori2,c_Ori1};
            axes(hPlotsFig1.hPlot1(c_Ori2,c_Ori1)); %#ok<LAXES>
            hold(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),'on');
            subplot(hPlotsFig1.hPlot1(c_Ori2,c_Ori1))
            tickLen = rasterScale/length(X);
            rasterplot(X,(tickLen:tickLen:rasterScale)-tickLen/2,'k',tickLen);
            plot(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),firingRateData.timeVals,squeeze(psthData2(c_Ori2,c_Ori1,:)),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
            set(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        end
    end
    
else
    psthData = squeeze(mean(firingRateData.data,1));
    for c_Ori2 = 1:5
        for c_Ori1 = 1:5
            plot(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),firingRateData.timeVals,squeeze(psthData(c_Ori2,c_Ori1 ,:)),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
            set(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        end
    end
end
for cOri1 = 1:5
    title(hPlotsFig1.hPlot1(1,cOri1),[num2str(cValsUnique(cOri1)) ' %'])
end
cValsUniqueFlipped = flip(cValsUnique);
figure(1)
for cOri2 = 1:5
    textH{cOri2} = getPlotHandles(1,1,[0.81 0.82-(cOri2-1)*0.115 0.01 0.01]);
end

textString = {'50 %','25 %','12.5 %','6.25 %','0 %'};

for cOri2 = 1:5
    set(textH{cOri2},'Visible','Off')
    text(0.35,1.15,textString{cOri2},'unit','normalized','fontsize',16,'fontweight','bold','parent',textH{cOri2});
end
if ~isempty(elecNum)
    textString1 =['Contrast of Grating 1 (Orientation: 112.5 ' char(176) ')'];  textString2 = ['Contrast of Grating 2 (Orientation: 22.5 ' char(176) ')'];
    textH1 = getPlotHandles(1,1,[0.36 0.97 0.01 0.01]);
    textH2 = getPlotHandles(1,1,[0.88 0.85 0.01 0.01]);
    set(textH1,'Visible','Off');set(textH2,'Visible','Off');
    text(0.35,1.15,textString1,'unit','normalized','fontsize',16,'fontweight','bold','parent',textH1);
    text(0.35,1.15,textString2,'unit','normalized','fontsize',16,'fontweight','bold','rotation',270,'parent',textH2);
end
rescaleData(hPlotsFig1.hPlot1,-0.1,0.5,[-5 0]+ getYLims(hPlotsFig1.hPlot1),14);


for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        if ~isempty(elecNum)
            plot(hPlotsFig1.hPlot2(1,c_Ori1),firingRateData.timeVals,squeeze(psthData2(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        else
            plot(hPlotsFig1.hPlot2(1,c_Ori1),firingRateData.timeVals,squeeze(psthData(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        end
        hold(hPlotsFig1.hPlot2(1,c_Ori1),'on');
        %         plot(hPlotsFig1.hPlot2(2,c_Ori1),firingRateData.timeVals,squeeze(combinedData{1}(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        %         hold(hPlotsFig1.hPlot2(2,c_Ori1),'on');
    end
end

for i = 1:length(cValsUnique)
    set(hPlotsFig1.hPlot2(1,i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    %     set(hPlotsFig1.hPlot2(2,i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    
end
xlabel(hPlotsFig1.hPlot1(5,1),'Time (s)'); ylabel(hPlotsFig1.hPlot2(1,1),[{'Firing Rate'} {'(spikes/s)'}]);

xlabel(hPlotsFig1.hPlot2(1,1),'Time (s)'); ylabel(hPlotsFig1.hPlot2(1,1),[{'Firing Rate'} {'(spikes/s)'}]);
rescaleData(hPlotsFig1.hPlot2,-0.1,0.5,[-5 0]+getYLims(hPlotsFig1.hPlot2),14);
rescaleData(hPlotsFig1.hPlot1,-0.1,0.5,[-5 0]+ getYLims(hPlotsFig1.hPlot1),14);


% Figure 2

if strcmp(colorScheme,'color')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end
% dEnergyTF = 10*(energyDataTF.data - energyDataTF.data_cBL);
% energyVsFrequencyDataST = squeeze(mean(energyData.dataST(elecNum,1,:,:,:),1));
% energyVsFrequencyDataBL = squeeze(mean(energyData.data_cBL(elecNum,1,:,:,:),1));
% dEnergyVsFrequencyData = 10*(energyVsFrequencyDataST-energyVsFrequencyDataBL);
dEnergyTF = 10*combinedData{2};
% energyVsFrequencyDataST = squeeze(mean(energyData.dataST(elecNum,1,:,:,:),1));
% energyVsFrequencyDataBL = squeeze(mean(energyData.data_cBL(elecNum,1,:,:,:),1));
energyVsFrequencyDataBL = combinedData{4};
dEnergyVsFrequencyData = 10*combinedData{3};


if ~isempty(elecNum)
    dEnergyTF2 = dEnergyTF;
    dEnergyVsFrequencyData2 = dEnergyVsFrequencyData;
    energyVsFrequencyDataBL2 = energyVsFrequencyDataBL;
else
    dEnergyTF2 = squeeze(mean(dEnergyTF,1));
    dEnergyVsFrequencyData2 = squeeze(mean(dEnergyVsFrequencyData,1));
    energyVsFrequencyDataBL2 = squeeze(mean(energyVsFrequencyDataBL,1));
    
end
% if elecNum == 29 || elecNum == 40 || elecNum == 169
%     dEnergyTF2 = flip(flip(permute(dEnergyTF,[1 2 4 3 5 6]),4),3);
%     dEnergyVsFrequencyData2 = flip(flip(permute(dEnergyVsFrequencyData,[2 1 3]),2),1);
% end

for c_Ori2 = 1:5
    for c_Ori1 = 1:5
        pcolor(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),energyDataTF.timeVals,energyDataTF.freqVals,squeeze(dEnergyTF2(c_Ori2,c_Ori1 ,:,:)));
        caxis(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),[-2 10]);
        shading(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'interp'); %caxis([-5 10]);
        xLims = get(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'XLim');
        p1 = line(xLims,[freqRanges{2}(1) freqRanges{2}(1)],'parent',hPlotsFig2.hPlot1(c_Ori2,c_Ori1));
        p2 = line(xLims,[freqRanges{2}(2) freqRanges{2}(2)],'parent',hPlotsFig2.hPlot1(c_Ori2,c_Ori1));
        set(p1,'color','k','LineWidth',2);
        set(p2,'color','k','LineWidth',2);
        
        p3 = line(xLims,[freqRanges{3}(1) freqRanges{3}(1)],'parent',hPlotsFig2.hPlot1(c_Ori2,c_Ori1));
        p4 = line(xLims,[freqRanges{3}(2) freqRanges{3}(2)],'parent',hPlotsFig2.hPlot1(c_Ori2,c_Ori1));
        set(p3,'color',[0.45 0.45 0.45],'LineStyle','--','LineWidth',2);
        set(p4,'color',[0.45 0.45 0.45],'LineStyle','--','LineWidth',2);
        
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
    plot(hPlotsFig2.hPlot2(1,c_Ori1),energyData.freqVals,squeeze(energyVsFrequencyDataBL2(cFlipped_Indices(c_Ori1),c_Ori1,:)-energyVsFrequencyDataBL2(cFlipped_Indices(c_Ori1),c_Ori1,:)),'color','k','LineWidth',2);
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

for cOri1 = 1:5
    title(hPlotsFig2.hPlot1(1,cOri1),[num2str(cValsUnique(cOri1)) ' %'])
end
cValsUniqueFlipped = flip(cValsUnique);

figure(2)
for cOri2 = 1:5
    textH{cOri2} = getPlotHandles(1,1,[0.88 0.83-(cOri2-1)*0.115 0.01 0.01]);
end

textString = {'50 %','25 %','12.5 %','6.25 %','0 %'};

for cOri2 = 1:5
    set(textH{cOri2},'Visible','Off')
    text(0.35,1.15,textString{cOri2},'unit','normalized','fontsize',16,'fontweight','bold','parent',textH{cOri2});
end
if ~isempty(elecNum)
    textString1 =['Contrast of Grating 1 (Orientation: 112.5 ' char(176) ')'];  textString2 = ['Contrast of Grating 2 (Orientation: 22.5 ' char(176) ')'];
    figure(2)
    textH1 = getPlotHandles(1,1,[0.38 0.97 0.01 0.01]);
    textH2 = getPlotHandles(1,1,[0.95 0.83 0.01 0.01]);
    set(textH1,'Visible','Off');set(textH2,'Visible','Off');
    text(0.35,1.15,textString1,'unit','normalized','fontsize',16,'fontweight','bold','parent',textH1);
    text(0.35,1.15,textString2,'unit','normalized','fontsize',16,'fontweight','bold','rotation',270,'parent',textH2);
end
% xLims = get(hPlotsFig2.hPlot(1,1),'XLim');


rescaleData(hPlotsFig2.hPlot2,0,250,getYLims(hPlotsFig2.hPlot2),14);
for i=1:5
    displayRange(hPlotsFig2.hPlot2(1,i),[freqRanges{2}(1) freqRanges{2}(2)],getYLims(hPlotsFig2.hPlot2),'k');
    displayRange(hPlotsFig2.hPlot2(1,i),[freqRanges{3}(1) freqRanges{3}(2)],getYLims(hPlotsFig2.hPlot2),[0.5 0.5 0.5]);
end

end

% combineData
function xAll = combineData(firingRateData,energyData,energyDataTF,elecInfo,elecNum)
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

if ~isempty(elecNum)
    % allData = cell(1,4);
    for i= 1:size(goodIDList,2)
        a = find(elecNum==goodIDList{i}); %#ok<*EFIND>
        if ~isempty(a)
            UniqueElecID = i;
        end
    end
    
    clear allDataTMP alldTFDataTMP alldPSDDataTMP allblPSDDataTMP
    % allDataTMP = zeros(numGoodElectrodes,5,5);
    
    x = squeeze(firingRateData.data(goodIDList{UniqueElecID},1,:,:,:));
    y = squeeze(energyDataTF.data(goodIDList{UniqueElecID},1,:,:,:,:))-squeeze(energyDataTF.data_cBL(goodIDList{UniqueElecID},1,:,:,:,:));
    z = squeeze(energyData.dataST(goodIDList{UniqueElecID},1,:,:,:))-squeeze(energyData.data_cBL(goodIDList{UniqueElecID},1,:,:,:));
    
    b = squeeze(energyData.data_cBL(goodIDList{UniqueElecID},1,:,:,:));
    
    if length(goodIDList{UniqueElecID})==1
        allDataTMP = x;
        alldTFDataTMP = y;
        alldPSDDataTMP = z;
        allblPSDDataTMP = b;
    else
        xs = zeros(size(x));
        ys = zeros(size(y));
        zs = zeros(size(z));
        bs = zeros(size(b));
        for j=1:length(goodIDList{UniqueElecID})
            xs(j,:,:,:) = squeeze(x(j,:,:,:));
            ys(j,:,:,:,:) = squeeze(y(j,:,:,:,:));
            zs(j,:,:,:) = squeeze(z(j,:,:,:));
            bs(j,:,:,:) = squeeze(b(j,:,:,:));
        end
        allDataTMP = squeeze(mean(xs,1));
        alldTFDataTMP = squeeze(mean(ys,1));
        alldPSDDataTMP = squeeze(mean(zs,1));
        allblPSDDataTMP = squeeze(mean(bs,1));
    end
    
else
    clear allDataTMP alldTFDataTMP alldPSDDataTMP allblPSDDataTMP
    
    for i=1:numGoodElectrodes
        x = squeeze(firingRateData.data(goodIDList{i},1,:,:,:));
        y = squeeze(energyDataTF.data(goodIDList{i},1,:,:,:,:))-squeeze(energyDataTF.data_cBL(goodIDList{i},1,:,:,:,:));
        z = squeeze(energyData.dataST(goodIDList{i},1,:,:,:))-squeeze(energyData.data_cBL(goodIDList{i},1,:,:,:));
        
        b = squeeze(energyData.data_cBL(goodIDList{i},1,:,:,:));
        
        if length(goodIDList{i})==1
            allDataTMP(i,:,:,:) = x;
            alldTFDataTMP(i,:,:,:,:) = y;
            alldPSDDataTMP(i,:,:,:) = z;
            allblPSDDataTMP(i,:,:,:) = b;
        else
            xs = zeros(size(x));
            ys = zeros(size(y));
            zs = zeros(size(z));
            bs = zeros(size(b));
            for j=1:length(goodIDList{i})
                xs(j,:,:,:) = squeeze(x(j,:,:,:));
                ys(j,:,:,:,:) = squeeze(y(j,:,:,:,:));
                zs(j,:,:,:) = squeeze(z(j,:,:,:));
                bs(j,:,:,:) = squeeze(b(j,:,:,:));
            end
            allDataTMP(i,:,:,:) = squeeze(mean(xs,1));
            alldTFDataTMP(i,:,:,:,:) = squeeze(mean(ys,1));
            alldPSDDataTMP(i,:,:,:) = squeeze(mean(zs,1));
            allblPSDDataTMP(i,:,:,:) = squeeze(mean(bs,1));
        end
    end
end

xAll{1} = allDataTMP; xAll{2} = alldTFDataTMP; xAll{3} = alldPSDDataTMP; xAll{4} = allblPSDDataTMP;
type{1} = 'PSTH';type{2} = 'dTF'; type{3} = 'dPSD';type{3} = 'blPSD'; %#ok<*NASGU>

% for i=2:4
%     eDataTMP = squeeze(energyData.dataST{i}) - squeeze(energyData.data_cBL{i});
%
%     allDataTMP = zeros(numGoodElectrodes,5,5);
%     for j=1:numGoodElectrodes
%         clear x
%         x = eDataTMP(goodIDList{j},:,:);
%
%         xs = zeros(size(x));
%         for k=1:length(goodIDList{j})
% %             y = squeeze(x(k,:,:));
% %             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
%             xs(k,:,:) = squeeze(x(k,:,:));
%         end
%
%         allDataTMP(j,:,:) = squeeze(mean(xs,1));
%     end
%     xAll{i} = allDataTMP;
% end
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
function [ElectrodeArrayListAll,allGoodElectrodesStrArray,ElectrodeStringListAll,allGoodNs,allGoodSNRs,alld_elecs,numElecs] = getElectrodesList(fileNameStringTMP,elecParams,timeRangeForComputation,folderSourceString)

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
end

% function [ElectrodeArrayListAll,numElecs] = getElectrodesList(fileNameStringTMP,elecParams,timeRangeForComputation,folderSourceString)
%
% % [~,tmpElectrodeArrayList,~] = getGoodElectrodesDetails(fileNameStringTMP,oriSelectiveFlag,folderSourceString);
%
% gridType = 'microelectrode';
%
% numSessions = length(fileNameStringTMP);
% tmpElectrodeStringList = cell(1,numSessions);
% tmpElectrodeArrayList = cell(1,numSessions);
% numElecs = 0;
%
% Monkey1_ExpDates = dataInformationPlaidNorm('alpaH',gridType,0);
% Monkey1_SessionNum = length(Monkey1_ExpDates);
% % Monkey2_ExpDates = dataInformationPlaidNorm('kesariH',gridType,0);
% % Monkey2_SessionNum = length(Monkey2_ExpDates);
%
% for i = 1:numSessions
%     clear monkeyName
%     if strcmp(fileNameStringTMP{i}(1:5),'alpaH')
%         monkeyName = 'alpaH';
%         expDate = fileNameStringTMP{i}(6:11);
%         protocolName = fileNameStringTMP{i}(12:end);
%     elseif strcmp(fileNameStringTMP{i}(1:7),'kesariH')
%         monkeyName = 'kesariH';
%         expDate = fileNameStringTMP{i}(8:13);
%         protocolName = fileNameStringTMP{i}(14:end);
%     end
%     if i == 1
%         disp(['MonkeyName: ' ,monkeyName])
%     elseif i == Monkey1_SessionNum+1 % 13 Sessions are from alpaH; 9 Sessions from kesariH;
%         disp(['MonkeyName: ' ,monkeyName])
%     end
%     versionNum = 2;
%     [tmpElectrodeStringList{i},tmpElectrodeArrayList{i},goodElectrodes] = getGoodElectrodesSingleSession(monkeyName,expDate,protocolName,gridType,elecParams,timeRangeForComputation,folderSourceString,versionNum);
%     numElecs = numElecs+length(goodElectrodes);
% end
%
% ElectrodeArrayListAll = tmpElectrodeArrayList;
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
        plot(plotHandles(i,j),xVals1,yVals,'color',colorName,'LineWidth',1.5);
        plot(plotHandles(i,j),xVals2,yVals,'color',colorName,'LineWidth',1.5);
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
function plotFigures(monkeyName,folderSourceString)
if ~exist('folderSourceString','var') 
   folderSourceString = 'M:\Data\PlaidNorm\';
end
close all; % closes any open figure to avoid any overlaying issues

% Variable Parameters
spikeCutoff = 20;
snrCutoff = 2;
timeRangeForComputation = [0.15 0.4]; % expressed in second
timeRangeForComputationBL = [-diff(timeRangeForComputation) 0];
dRange = [0 0.75];
tapers_MT = [1 1]; % parameters for MT analysis
removeERPFlag = 1;

% Fixed parameters
folderSave = fullfile(folderSourceString,'Projects\PlaidNormalizationProject\savedData_Figures');
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
hPlotsFig1.hPlot2 = getPlotHandles(1,3,[0.15 0.2 0.7 0.3],0.1,0.05,0);

% Figure 2 (Spike data from all Elecs)
hFigure2 = figure(2);
set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig2.hPlot1 = getPlotHandles(1,5,[0.15 0.65 0.7 0.2],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot1);
hPlotsFig2.hPlot2 = getPlotHandles(1,3,[0.15 0.2 0.7 0.3],0.1,0.05,0);
% 
% % Figure 3 (PSDs (A); deltaPSDs (B) along increasing contrast of Ori 1 
% % with contrasts of Ori 2 presented in different colors;
% % absolute color map 5x5; relative color map 5x5; 
% % CRF (no Norm, Pref,Pref+Null, Avg, Null) 
% % for alpha (C), gamma (D) and high-gamma (E) for ERP-subtracted data
hFigure3 = figure(3);
set(hFigure3,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig3.hPlot1 = getPlotHandles(2,5,[0.25 0.65 0.5 0.3],0.01,0.01,1); linkaxes(hPlotsFig3.hPlot1);
hPlotsFig3.hPlot2 = getPlotHandles(3,3,[0.25 0.05 0.5 0.55],0.1,0.01,0);
% 
% % Figure 4 (PSDs (A); deltaPSDs (B) along increasing contrast of Ori 1 
% % with contrasts of Ori 2 presented in different colors;
% % absolute color map 5x5; relative color map 5x5; 
% % CRF (no Norm, Pref,Pref+Null, Avg, Null) 
% % for SSVEP (C) for non-ERP subtracted data
hFigure4 = figure(4);
set(hFigure4,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig4.hPlot1 = getPlotHandles(2,5,[0.15 0.5 0.7 0.4],0.01,0.01,1); linkaxes(hPlotsFig4.hPlot1);
hPlotsFig4.hPlot2 = getPlotHandles(1,3,[0.15 0.1 0.7 0.3],0.1,0.05,0);
% 
% 
% % Figure 5: Comparison of fitted parameters for FR, alpha, gamma, high
% % gamma, SSVEP--- Three plots for each NI, sigma and alpha fitted
% % parameters
% hFigure5 = figure(5);
% set(hFigure5,'units','normalized','outerposition',[0 0 1 1])
% hPlotsFig5.hPlot1 = getPlotHandles(5,3,[0.3 0.1 0.4 0.8],0.06,0.01,1); linkaxes(hPlotsFig4.hPlot1);

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
else
%     oriSelectiveFlag = 0;
    % get Data all Session for particular monkey or both combined for all
    % Electrodes
    [erpData,firingRateData,fftData,energyData,~,NI_Data,~] = ...
    getData(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); %#ok<ASGLU>
    save(fileSave1,'erpData','firingRateData','fftData','energyData','NI_Data')
end

% Put plot Functions for figures 1,3
plotData_spikes(hPlotsFig1,firingRateData,0) % spikes for static gratings, Fig 1
rescaleData(hPlotsFig1.hPlot1,-0.1,0.5,getYLims(hPlotsFig1.hPlot1));
rescaleData(hPlotsFig1.hPlot2(3),0,50,getYLims(hPlotsFig1.hPlot2(3)));

plotData_energy(hPlotsFig3,energyData) % alpha, gamma, hi-gamma for static gratings, Fig 3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Get data for Orientation selective Good Electrodes %%%%%%%%%%%
disp('Orientation-Tuned Electrodes:')
oriSelectiveFlag = 1; 
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
else
%     oriSelectiveFlag = 1;
    % get Data all Session for particular monkey or both combined for 
    % ori-tuned Electrodes
    [erpData,firingRateData,fftData,energyData,~,NI_Data,~] = ...
    getData(folderSourceString,fileNameStringListAll,electrodeList_OriTuned,timeRangeParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); %#ok<ASGLU>
    save(fileSave2,'erpData','firingRateData','fftData','energyData','NI_Data')
end

% Plotting Functions
plotData_spikes(hPlotsFig2,firingRateData,0) % spikes for static gratings from ori-selective electrodes, Fig 1
rescaleData(hPlotsFig2.hPlot1,-0.1,0.5,getYLims(hPlotsFig2.hPlot1));
rescaleData(hPlotsFig2.hPlot2(3),0,50,getYLims(hPlotsFig2.hPlot2(3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

removeERPFlag = 0; % Evoked Response
fileSave3 = fullfile(folderSave,[monkeyName '_N' num2str(spikeCutoff) '_S' num2str(snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(dRange(1)) '_' num2str(dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(combineUniqueElectrodeData) ...
    '_gse' num2str(getSpikeElectrodesFlag) '_gridType_' gridType '_UnitID' num2str(unitID) '.mat']); 
if exist(fileSave1,'file')
    disp(['Loading file ' fileSave3]);
    load(fileSave3);
else
    oriSelectiveFlag = 0;
    % get Data all Session for particular monkey or both combined for all
    % Electrodes
    [erpData,firingRateData,fftData,energyData,~,NI_Data,~] = ...
    getData(folderSourceString,fileNameStringListAll,electrodeList_All,timeRangeParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); %#ok<ASGLU>
    save(fileSave1,'erpData','firingRateData','fftData','energyData','NI_Data')
end

plotData_energy(hPlotsFig4,energyData) % SSVEP Evoked, Fig 4; 


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,electrodeArray] = ...
    getData(folderSourceString,fileNameStringTMP,ElectrodeListTMP,dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod)

numDatasets = length(fileNameStringTMP);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,electrodeArray]...
= getDataSingleSession(folderSourceString,fileNameStringTMP{1},...
ElectrodeListTMP{1},dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); 

if length(fileNameStringTMP)>1
    for i=2:numDatasets
        if isempty(ElectrodeListTMP{i}{end})
           continue
        end
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [erpDataTMP,firingRateDataTMP,fftDataTMP,energyDataTMP,~,NI_DataTMP,electrodeArrayTMP] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},...
            ElectrodeListTMP{i},dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod);
        
        erpData.data = cat(1,erpData.data,erpDataTMP.data);
        erpData.analysisDataBL = cat(1,erpData.analysisDataBL,erpDataTMP.analysisDataBL);
        erpData.analysisData_cBL = cat(1,erpData.analysisData_cBL,erpDataTMP.analysisData_cBL);
        erpData.analysisDataST = cat(1,erpData.analysisDataST,erpDataTMP.analysisDataST);
        
        firingRateData.data = cat(1,firingRateData.data,firingRateDataTMP.data);
        firingRateData.analysisDataBL = cat(1,firingRateData.analysisDataBL,firingRateDataTMP.analysisDataBL);
        firingRateData.analysisData_cBL = cat(1,firingRateData.analysisData_cBL,firingRateDataTMP.analysisData_cBL);
        firingRateData.analysisDataST = cat(1,firingRateData.analysisDataST,firingRateDataTMP.analysisDataST);
        
        fftData.dataBL = cat(1,fftData.dataBL,fftDataTMP.dataBL);
        fftData.data_cBL = cat(1,fftData.data_cBL,fftDataTMP.data_cBL);
        fftData.dataST = cat(1,fftData.dataST,fftDataTMP.dataST);
        for j = 1:length(fftData.analysisDataBL)
            fftData.analysisDataBL{j} = cat(1,fftData.analysisDataBL{j},fftDataTMP.analysisDataBL{j});
            fftData.analysisData_cBL{j} = cat(1,fftData.analysisData_cBL{j},fftDataTMP.analysisData_cBL{j});
            fftData.analysisDataST{j} = cat(1,fftData.analysisDataST{j},fftDataTMP.analysisDataST{j});
        end
        
        energyData.dataBL = cat(1,energyData.dataBL,energyDataTMP.dataBL);
        energyData.data_cBL = cat(1,energyData.data_cBL,energyDataTMP.data_cBL);
        energyData.dataST = cat(1,energyData.dataST,energyDataTMP.dataST);
        for j =1:length(energyData.analysisDataBL)
            energyData.analysisDataBL{j} = cat(1,energyData.analysisDataBL{j},energyDataTMP.analysisDataBL{j});
            energyData.analysisData_cBL{j} = cat(1,energyData.analysisData_cBL{j},energyDataTMP.analysisData_cBL{j});
            energyData.analysisDataST{j} = cat(1,energyData.analysisDataST{j},energyDataTMP.analysisDataST{j});
        end
        
        NI_Data.erp = cat(1, NI_Data.erp, NI_DataTMP.erp);
        NI_Data.firingRate = cat(1, NI_Data.firingRate, NI_DataTMP.firingRate);
        for j =1:length(NI_Data.fft)
            NI_Data.fft{j} = cat(2, NI_Data.fft{j}, NI_DataTMP.fft{j});
            NI_Data.energy{j} = cat(2, NI_Data.energy{j}, NI_DataTMP.energy{j});
        end

        % Combining OriData across sessions may be required! Right now, there is no requirement!        
        electrodeArray = cat(2,electrodeArray,electrodeArrayTMP);
    end
end
end

function [erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,electrodeArray] = ...
    getDataSingleSession(folderSourceString,fileNameStringTMP,...
    ElectrodeListTMP,dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod)

gridType = 'microelectrode';
if strcmp(fileNameStringTMP(1:5),'alpaH')       
    monkeyName = 'alpaH'; 
    expDate = fileNameStringTMP(6:11); 
    protocolName = fileNameStringTMP(12:end);
    oriTuning_protocolName = ['GRF_00' num2str(str2double(protocolName(5:end))-1)]; % The protocol Number is just the immediate precedent of the main protocol 
elseif strcmp(fileNameStringTMP(1:7),'kesariH')
    monkeyName = 'kesariH';
    expDate = fileNameStringTMP(8:13); 
    protocolName = fileNameStringTMP(14:end);
    oriTuning_protocolName = ['GRF_00' num2str(str2double(protocolName(5:end))-1)]; % The protocol Number is just the immediate precedent of the main protocol 
end


folderName = fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,protocolName);
tuningProtocol_folderName =  fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,oriTuning_protocolName);                           

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');
folderSave = fullfile(folderSourceString,'Projects\PlaidNormalizationProject\savedDataV2');
folderSave_oriTuning = fullfile(tuningProtocol_folderName,'savedData');

if ~exist(folderSave,'dir')
    mkdir(folderSave);
end
if ~exist(folderSave_oriTuning,'dir')
    mkdir(folderSave_oriTuning);
end

% Load Orientation Tuning dataFile for ori Tuning protocol for all elecs
oriTuningDataFile = fullfile(folderSave_oriTuning,['oriTuningData_' num2str(1000*dataParameters.stRange(1)) 'ms_' num2str(1000*dataParameters.stRange(2)) 'ms.mat']);
if exist(oriTuningDataFile,'file')
    disp(['Loading file ' oriTuningDataFile]);
    load(oriTuningDataFile);
else
    % Get OrientationTuning Data
    [computationVals,PO,OS] = savePrefOriAndOriSelectivitySpikes(monkeyName,expDate,oriTuning_protocolName,folderSourceString,gridType);
end

oriTuningData.PO = PO(ElectrodeListTMP{end});
oriTuningData.OS = OS(ElectrodeListTMP{end});
oriTuningData.FR = computationVals(ElectrodeListTMP{end},:);

if oriSelectiveFlag 
    fileToSave = fullfile(folderSave,[fileNameStringTMP '_OriTunedElecData_StimPeriod_' num2str(1000*dataParameters.stRange(1)) '_' num2str(1000*dataParameters.stRange(2)) 'ms_tapers' num2str(tapers_MT(1)) '_' num2str(tapers_MT(2)) '_' strtok(LFPdataProcessingMethod) '.mat']);
else
    fileToSave = fullfile(folderSave,[fileNameStringTMP '_allElecData_StimPeriod_' num2str(1000*dataParameters.stRange(1)) '_' num2str(1000*dataParameters.stRange(2)) 'ms_tapers' num2str(tapers_MT(1)) '_' num2str(tapers_MT(2)) '_' strtok(LFPdataProcessingMethod) '.mat']);
end

if exist(fileToSave,'file')
    disp(['Loading file ' fileToSave]);
    load(fileToSave)
else
    % Get Combinations
    [parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,~,~,oValsUnique,cValsUnique,tValsUnique, ...
    aValsUnique2,eValsUnique2,~,~,oValsUnique2,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract);

    if aValsUnique ~= aValsUnique2 || eValsUnique ~= eValsUnique2
        error('Azimuths and/or elevations do not match!');
    end                                          
    a=1; e=1; s=1; f=1; o=1; 
    if tValsUnique ~= tValsUnique2
        error('Azimuths and/or elevations do not match!');
    else
        tList = 1:length(tValsUnique);
    end

    % electrode info
    ElectrodeList = ElectrodeListTMP{end};
    electrodeArray = ElectrodeList;

    % TimeVals info
    [~,timeVals,~,~] = loadlfpInfo(folderLFP);

    % Set up fft
    Fs = round(1/(timeVals(2)-timeVals(1)));
    range = dataParameters.blRange;
    rangePos = round(diff(range)*Fs);
    erpRangePos = round(diff(dataParameters.erpRange)*Fs);
    blPos = find(timeVals>=dataParameters.blRange(1),1)+ (1:rangePos);
    stPos = find(timeVals>=dataParameters.stRange(1),1)+ (1:rangePos);
    erpPos = find(timeVals>=dataParameters.erpRange(1),1)+ (1:erpRangePos);
    freqVals = 0:1/diff(range):Fs-1/diff(range);
    numFreqs = length(freqRanges);
    unitID = 0;

    % Set up params for MT
    params.tapers   = tapers_MT;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 250];
    params.trialave = 1;

%     cList_Ori1 = 1:length(cValsUnique); 
    cListFlipped_Ori2 = flip(1:length(cValsUnique2)); 

    % Main Loop (Stores data in elec x tempFreq x Contrast of Ori 2 x
    % Contrast of Ori 1 x dataPoints)

    for iElec = 1:length(ElectrodeList)
        % Get LFP data
        clear analogData
        load(fullfile(folderLFP,['elec' num2str(ElectrodeList(iElec)) '.mat']));
        % Get Spike data
        clear spikeData
        load(fullfile(folderSpikes,['elec' num2str(ElectrodeList(iElec)) '_SID' num2str(unitID) '.mat']));

        % Get bad trials
        badTrialFile = fullfile(folderSegment,'badTrials.mat');
        if ~exist(badTrialFile,'file')
            disp('Bad trial file does not exist...');
            badTrials=[];
        else
            badTrials = loadBadTrials(badTrialFile);
            if iElec == 1
                disp([num2str(length(badTrials)) ' bad trials']);
            end
        end

        disp(['Processing electrode (' num2str(iElec) ') : ' num2str(ElectrodeList(iElec))]);

        for t = 1:length(tList)
            for c_Ori2 = 1:length(cValsUnique2)
                for c_Ori1 = 1:length(cValsUnique)

                    clear goodPos fftBL fftST erp dataBL dataST
                    goodPos = parameterCombinations{a,e,s,f,o,c_Ori1,tList(t)};
                    goodPos2 = parameterCombinations2{a,e,s,f,o,cListFlipped_Ori2(c_Ori2),tList(t)};
                    goodPos = intersect(goodPos,goodPos2);
                    goodPos = setdiff(goodPos,badTrials);

                    if isempty(goodPos)
                        disp('No entries for this combination..');
                    else
%                         disp(['pos=(' num2str(cListFlipped_Ori2(c_Ori2)) ',' num2str(c_Ori1) ') ,n=' num2str(length(goodPos))]);
                        N(iElec,t,c_Ori2,c_Ori1) = length(goodPos);
                        if round(diff(dataParameters.blRange)*Fs) ~= round(diff(dataParameters.stRange)*Fs)
                            disp('baseline and stimulus ranges are not the same');
                        else
                            
                           if t == 1 % ERP and firing rate data processed only for static stimuli 
                               % Event-related potential
                               erp = mean(analogData(goodPos,:),1); %#ok<NODEF>
                               erpDataTMP(iElec,t,c_Ori2,c_Ori1,:) = erp;
                               RMSvalsBL(iElec,t,c_Ori2,c_Ori1) = rms(erp(blPos));
                               RMSvalsERP(iElec,t,c_Ori2,c_Ori1) = rms(erp(erpPos));

                               % PSTH & firing rate 
                               [psthData(iElec,t,c_Ori2,c_Ori1,:),xsFR] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
                               firingRatesBL(iElec,t,c_Ori2,c_Ori1) = mean(getSpikeCounts(spikeData(goodPos),dataParameters.blRange))/diff(dataParameters.blRange);
                               firingRatesST(iElec,t,c_Ori2,c_Ori1) = mean(getSpikeCounts(spikeData(goodPos),dataParameters.stRange))/diff(dataParameters.stRange);
                           end
                           
                           % fft data is processed for both static and
                           % flickering stimuli
                           if strcmp(strtok(LFPdataProcessingMethod),'Evoked')
                               fftBL = squeeze(mean(abs(fft(analogData(goodPos,blPos),[],2))));
                               fftST = squeeze(mean(abs(fft(analogData(goodPos,stPos),[],2))));
                           elseif strcmp(strtok(LFPdataProcessingMethod),'Induced')
                               fftBL = squeeze(mean(abs(fft(removeERP(analogData(goodPos,blPos)),[],2))));
                               fftST = squeeze(mean(abs(fft(removeERP(analogData(goodPos,stPos)),[],2))));
                           end
                           fftDataBL(iElec,t,c_Ori2,c_Ori1,:) = conv2Log(fftBL); %#ok<*AGROW>
                           fftDataST(iElec,t,c_Ori2,c_Ori1,:) = conv2Log(fftST);
                           
                           % Power Estimation by MT method
                           % for both static and flickering stimuli
                           if strcmp(strtok(LFPdataProcessingMethod),'Evoked')
                               dataBL = analogData(goodPos,blPos)';
                               dataST = analogData(goodPos,stPos)';
                           elseif strcmp(strtok(LFPdataProcessingMethod),'Induced')
                               dataBL = removeERP(analogData(goodPos,blPos))';
                               dataST = removeERP(analogData(goodPos,stPos))';                           
                           end
                           [tmpEBL,freqValsBL] = mtspectrumc(dataBL,params);
                           [tmpEST,freqValsST] = mtspectrumc(dataST,params);
                           
                           if isequal(freqValsBL,freqValsST)
                               freqValsMT = freqValsST;
                           end
                           mEnergyVsFreqBL(iElec,t,c_Ori2,c_Ori1,:) = conv2Log(tmpEBL);
                           mEnergyVsFreqST(iElec,t,c_Ori2,c_Ori1,:) = conv2Log(tmpEST);

                           % computing analysis Data for particular
                           % frequency band
                           
                           if t == 1
                               for i=1:numFreqs-1
                                   fftAmpBL{i}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(fftBL(:),freqVals,freqRanges{i}));
                                   fftAmpST{i}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(fftST(:),freqVals,freqRanges{i}));
                                   energyValsBL{i}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),freqValsMT,freqRanges{i}));
                                   energyValsST{i}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),freqValsMT,freqRanges{i}));
                               end
                               
                           elseif t == 2 %% fft and energy data for only SSVEP frequency
                               fftAmpBL{numFreqs}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(fftBL(:),freqVals,freqRanges{numFreqs}));
                               fftAmpST{numFreqs}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(fftST(:),freqVals,freqRanges{numFreqs}));
                               energyValsBL{numFreqs}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),freqValsMT,freqRanges{numFreqs}));
                               energyValsST{numFreqs}(iElec,c_Ori2,c_Ori1,:) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),freqValsMT,freqRanges{numFreqs}));
                           end
                        end
                    end
                end
            end
        end
    end
    
    % display Stimulus repeats in a Table for Static Stimuli & Flickering
    % Stimuli
    for t = 1:length(tList)
        clear allElecStimReps
        allElecStimReps = N(:,t,:,:);
        count = 0;
        for iElec = 1:size(allElecStimReps,1)
            if isequal(squeeze(allElecStimReps(iElec,:,:,:)),squeeze(allElecStimReps(1,:,:,:)))
                count = count+1;
            else
                error('Stimulus repeats are not equal for different electrodes')
            end
        end
        
        if count == size(allElecStimReps,1)
            if t == 1
                disp('Stim Repeats for Static Stimuli, con_Ori2 (rows) x con_Ori1 (columns)')
            elseif t == 2  
                disp('Stim Repeats for Flickering Stimuli, con_Ori2 (rows) x con_Ori1 (columns)')
            end
            disp(squeeze(allElecStimReps(1,:,:,:)))
        end
    end
    
    
    % Time-Domain data
    erpData.data = erpDataTMP;
    erpData.analysisDataBL = RMSvalsBL;
    erpData.analysisData_cBL = getCommonBaseline(RMSvalsBL); % gets common (mean baseline) for all con_Ori2 x con_Ori1 contrast conditions 
    erpData.analysisDataST = RMSvalsERP;
    erpData.timeVals = timeVals;
    erpData.N = N;

    firingRateData.data = psthData;
    firingRateData.analysisDataBL = firingRatesBL;
    firingRateData.analysisData_cBL = getCommonBaseline(firingRatesBL); % gets common (mean baseline) for all con_Ori2 x con_Ori1 contrast conditions 
    firingRateData.analysisDataST = firingRatesST;
    firingRateData.timeVals = xsFR;
    firingRateData.N = N;

    % Freq-Domain data
    fftData.dataBL = fftDataBL;
    fftData.data_cBL = getCommonBaseline(fftDataBL);% gets common (mean baseline) for all 5 x 5 contrast conditions 
    fftData.dataST = fftDataST;
    fftData.analysisDataBL = fftAmpBL;
    fftData.analysisData_cBL = getCommonBaseline(fftAmpBL);% gets fft amplitude of common (mean baseline) for all 5 x 5 contrast conditions for selected freqRanges
    fftData.analysisDataST = fftAmpST;
    fftData.freqVals = freqVals;
    fftData.N = N;

    energyData.dataBL = mEnergyVsFreqBL;
    energyData.data_cBL = getCommonBaseline(mEnergyVsFreqBL);% gets common (mean baseline) for all 5 x 5 contrast conditions 
    energyData.dataST = mEnergyVsFreqST;
    energyData.analysisDataBL = energyValsBL;
    energyData.analysisData_cBL = getCommonBaseline(energyValsBL);% gets energyData of  common (mean baseline) for all 5 x 5 contrast conditions for selected freqRanges
    energyData.analysisDataST = energyValsST;
    energyData.freqVals = freqValsMT;
    energyData.N = N;

    % Segregation into Preferred-null axis is done when analysis is being
    % done for orientation selective electrodes
    if oriSelectiveFlag
        elecs_neededtoFlipped = find(abs(oriTuningData.PO-oValsUnique2)<abs(oriTuningData.PO-oValsUnique));
        if ~isempty(elecs_neededtoFlipped)
            erpData = segregate_Pref_Null_data(erpData,elecs_neededtoFlipped);
            firingRateData = segregate_Pref_Null_data(firingRateData,elecs_neededtoFlipped);
            fftData = segregate_Pref_Null_data(fftData,elecs_neededtoFlipped);
            energyData = segregate_Pref_Null_data(energyData,elecs_neededtoFlipped);
        end
    end

    
    % Get Normalization Indices
    NI_Data.erp = getNI(erpData);
    NI_Data.firingRate = getNI(firingRateData);
    NI_Data.fft = getNI(fftData);
    NI_Data.energy = getNI(energyData);

    % Save Data for particular session
    save(fileToSave,'erpData','firingRateData','fftData','energyData','oriTuningData','NI_Data','electrodeArray');
end
end



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
caxis(hPlot.hPlot2(1),[0 20]);
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
caxis(hPlot.hPlot2(2),[0 20]);
set(hPlot.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot.hPlot2(2),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(2),'Contrast of Ori 2(%)');


% CRF
errorbar(cValsUnique,spikeRateDataST(end,:),sem_spikeRate(end,:),...
    'Marker','o','LineWidth',2,'color',colors(end,:,:),'parent',hPlot.hPlot2(3))
hold(hPlot.hPlot2(3),'on');
errorbar(cValsUnique,diag(flipud(spikeRateDataST)),diag(flipud(sem_spikeRate)),'Marker','o','LineWidth',2,'color','k','parent',hPlot.hPlot2(3));
hold(hPlot.hPlot2(3),'off');
text(0.5,0.2,'cOri 2: 0%','color',colors(end,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
text(0.5,0.1,'cOri 1 = cOri 2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
set(hPlot.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot.hPlot2(3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlot.hPlot2(3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(3),'Absolute Spike Rate (spike/s)');

end
function plotData_energy(hPlot,data)

cValsUnique = [0 12.5 25 50 100]/2;
cValsUnique2 = [0 12.5 25 50 100]/2;
num_freqRanges = length(data.analysisDataST)-1; % alpha, gamma, hi-gamma Induced power is being plotted

% mean energy data across electrodes: con_Ori2 (rows) x con_Ori2 (columns) x
% freqVals
energyVsFrequencyDataST = squeeze(mean(data.dataST(:,1,:,:,:),1));
energyVsFrequencyDataBL = squeeze(mean(data.data_cBL(:,1,:,:,:),1));
dEnergyVsFrequencyData = energyVsFrequencyDataST-energyVsFrequencyDataBL;

% PSD plots
colors = jet(length(cValsUnique));
cFlipped_Indices = flip(1:length(cValsUnique2)); 

for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlot1.hPlot1(1,c_Ori1),data.freqVals,squeeze(energyVsFrequencyDataST(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot1.hPlot1(1,c_Ori1),'on');
        plot(hPlot1.hPlot1(2,c_Ori1),data.freqVals,squeeze(dEnergyVsFrequencyData(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlot1.hPlot1(2,c_Ori1),'on');
    end
end


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

    % remove Outlier elecs (add as a function)
%     OutlierVals = [-15 15];
%     NI_population_outlier = find(NI_population_spikeRateRelative<OutlierVals(1) | NI_population_spikeRateRelative>OutlierVals(2));
%     NI_population_outlierVals = NI_population_spikeRateRelative(NI_population_outlier);
%     NI_population_spikeRateRelative = NI_population_spikeRateRelative(setdiff(1:length(NI_population_spikeRateRelative),NI_population_outlier));
%     fprintf(['Deleting Electrode number: ',num2str(NI_population_outlier) ' \nfor NI calculation because NI value(s) '...
%         num2str(NI_population_outlierVals) '\nfalls outside range ' num2str(OutlierVals(1)) ' < NI values < ' num2str(OutlierVals(2)) '\n'] )



    % Color coded Plots of energyData

    imagesc(energyDataST,'parent',hPlot1.hPlot2(i,1));
    imagesc(dEnergyData,'parent',hPlot1.hPlot2(i,2));

    
    % CRF
    errorbar(cValsUnique,dEnergyData(end,:),sem_dEnergyData(end,:),...
        'Marker','o','LineWidth',2,'color',colors(end,:,:),'parent',hPlot1.hPlot2(i,3))
    hold(hPlot1.hPlot2(i,3),'on');
    errorbar(cValsUnique,diag(flipud(dEnergyData)),diag(flipud(sem_dEnergyData)),'Marker','o','LineWidth',2,'color','k','parent',hPlot1.hPlot2(i,3));
    hold(hPlot1.hPlot2(i,3),'off');
end


% set(hPlot.hPlot1(1),'XLim',[-0.1 0.5]);
% set(hPlot.hPlot1(1),'YLim',[0 90]);
% tickLengthPlot = 2*get(hPlot.hPlot1(1),'TickLength');
% xlabel(hPlot.hPlot1(1),'Time (s)')
% ylabel(hPlot.hPlot1(1),'Spike rate(spike/s)')
% displayRange(hPlot.hPlot1,[0.2 0.4],getYLims(hPlot.hPlot1),'k');
% 
% 
% for i = 1:length(cValsUnique)
%     set(hPlot.hPlot1(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
% end
% 
% 
% 
% colorBar_absSpikeRate = colorbar(hPlot.hPlot2(1)); 
% colorYlabelHandle = get(colorBar_absSpikeRate,'Ylabel');
% set(colorYlabelHandle,'String','Absolute Spike Rate (spikes/s)','fontSize',14);
% plotPos = get(hPlot.hPlot2(1),'Position');
% set(hPlot.hPlot2(1),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
% title(hPlot.hPlot2(1),['Mean NI: ',num2str(round(mean(NI_population_spikeRateAbsolute),2))],'fontWeight','bold');
% caxis(hPlot.hPlot2(1),[0 20]);
% set(hPlot.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
% set(hPlot.hPlot2(1),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
% xlabel(hPlot.hPlot2(1),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(1),'Contrast of Ori 2(%)');
% 
% 
% 
% colorBar_rlvSpikeRate = colorbar(hPlot.hPlot2(2)); 
% colorYlabelHandle = get(colorBar_rlvSpikeRate,'Ylabel');
% set(colorYlabelHandle,'String','Change in Spike Rate (spikes/s)','fontSize',14);
% plotPos = get(hPlot.hPlot2(2),'Position');
% set(hPlot.hPlot2(2),'Position',[plotPos(1) plotPos(2) plotPos(3)+0.02 plotPos(4)]);
% title(hPlot.hPlot2(2),['Mean NI: ',num2str(round(mean(NI_population_spikeRateRelative),2))],'fontWeight','bold');
% caxis(hPlot.hPlot2(2),[0 20]);
% set(hPlot.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
% set(hPlot.hPlot2(2),'XTick',1:length(cValsUnique),'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
% xlabel(hPlot.hPlot2(2),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(2),'Contrast of Ori 2(%)');
% 
% 
% 
% text(0.5,0.2,'cOri 2: 0%','color',colors(end,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
% text(0.5,0.1,'cOri 1 = cOri 2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlot.hPlot2(3))
% set(hPlot.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
% set(hPlot.hPlot2(3),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
% xlabel(hPlot.hPlot2(3),'Contrast of Ori 1(%)');ylabel(hPlot.hPlot2(3),'Absolute Spike Rate (spike/s)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Accessory Functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load LFP Info
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored); %#ok<NODEF>
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

% Get parameter combinations
function [parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,...
    cValsUnique,tValsUnique,aValsUnique2,eValsUnique2,sValsUnique2,...
    fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;            end

end

% Get Bad Trials
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
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
        fileNameStringListAll{pos} = fileNameStringList{i}{j};
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
    elseif i == 13
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

% Get Induced LFP data by subtracting trialaveraged ERP data from trialwise LFP Data
function Y = removeERP(X)
Y = X-repmat(mean(X,1),size(X,1),1);
end

% Get MeanEnergy for different frequency bands
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange)

posToAverage = intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2)));
eValue   = mean(mEnergy(posToAverage));
end

% Normalize data for ERP and Spike data
function normData = normalizeData(x)
for iElec = 1:size(x.data,1)
    for t = 1:size(x.data,2)
        normData.data(iElec,t,:,:,:) = x.data(iElec,t,:,:,:)./max(max(max(abs(x.data(iElec,t,:,:,:)))));
        normData.analysisDataBL(iElec,t,:,:) = x.analysisDataBL(iElec,t,:,:)./max(max(abs(x.analysisDataBL(iElec,t,:,:))));
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

% get Common Baseline across all 5 (Ori 1) x 5 (Ori 2) contrast conditions
function data_BL = getCommonBaseline(data_BL)

if iscell(data_BL)
    size_data_BL = numel(size(data_BL{1}));
    num_con_Ori2 = size(data_BL{1},3);
    num_con_Ori1 = size(data_BL{1},4);
else
    size_data_BL = numel(size(data_BL));
    num_con_Ori2 = size(data_BL,3);
    num_con_Ori1 = size(data_BL,4);
end

if size_data_BL == 4 % baseline for analysis data (elec x TF x Num_Contrast_Ori1 x Num_Contrast_Ori2 x (analysisDataBL))
    if iscell(data_BL)
        for iElec = 1:size(data_BL{1},1)
            for iTF = 1:size(data_BL{1},2)
                for k = 1:length(data_BL)
                    data_BL{k}(iElec,iTF,:,:) = repmat(mean(mean(squeeze(data_BL{k}(iElec,iTF,:,:)),2),1),[num_con_Ori2 num_con_Ori1]);
                end
            end
        end
    else
        for iElec = 1:size(data_BL,1)
            for iTF = 1:size(data_BL,2)
               data_BL(iElec,iTF,:,:) = repmat(mean(mean(squeeze(data_BL(iElec,iTF,:,:)),2),1),[num_con_Ori2 num_con_Ori1]);
            end
        end
    end
elseif size_data_BL == 5 % baseline for timeSeries/PSD data (elec x TF x Num_Contrast_Ori1 x Num_Contrast_Ori2 x time/FreqVals x (dataBL))
    for iElec = 1:size(data_BL,1)
        for iTF = 1:size(data_BL,2)
           data_BL(iElec,iTF,:,:,:) = reshape(repmat(squeeze(mean(mean(squeeze(data_BL(iElec,iTF,:,:,:)),2),1))',[num_con_Ori2*num_con_Ori1 1]),[num_con_Ori2 num_con_Ori1 size(data_BL,5)]);
        end
    end    
end
end

% Segregate data into Preferred (positive x-axis -- increasing contrast) and Null (positive y axis -- increasing contrast)
function data = segregate_Pref_Null_data(data,elecs_neededtoFlipped)
    for iElec = 1:length(elecs_neededtoFlipped)
        for iTF = 1:2
            disp ([num2str(iElec), ' ' num2str(iTF)])
            if numel(fieldnames(data)) == 6
                if iTF == 1
                    data.data(elecs_neededtoFlipped(iElec),:,:,:,:) = flip(flip(permute(squeeze(data.data(elecs_neededtoFlipped(iElec),:,:,:,:)),[2 1 3]),1),2);
                end
            elseif numel(fieldnames(data)) == 8
                data.dataBL(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.dataBL(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
                data.data_cBL(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.data_cBL(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
                data.dataST(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.dataST(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
            end
        end
    end
    
    if ~iscell(data.analysisDataBL) && ~iscell(data.analysisDataST)
        data.analysisDataBL(elecs_neededtoFlipped(iElec),:,:) = flip(flip(permute(squeeze(data.analysisDataBL(elecs_neededtoFlipped(iElec),:,:)),[2 1]),1),2);
        data.analysisData_cBL(elecs_neededtoFlipped(iElec),:,:) = flip(flip(permute(squeeze(data.analysisData_cBL(elecs_neededtoFlipped(iElec),:,:)),[2 1]),1),2);
        data.analysisDataST(elecs_neededtoFlipped(iElec),:,:) = flip(flip(permute(squeeze(data.analysisDataST(elecs_neededtoFlipped(iElec),:,:)),[2 1]),1),2);
    elseif iscell(data.analysisDataBL) && iscell(data.analysisDataST)
        for iCell = 1:length(data.analysisDataBL)
            data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),:,:) = flip(flip(permute(squeeze(data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),:,:)),[2 1]),1),2);
            data.analysisData_cBL{iCell}(elecs_neededtoFlipped(iElec),:,:) = flip(flip(permute(squeeze(data.analysisData_cBL{iCell}(elecs_neededtoFlipped(iElec),:,:)),[2 1]),1),2);
            data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),:,:) = flip(flip(permute(squeeze(data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),:,:)),[2 1]),1),2);
        end
    end
end

% computing NI for different neural measures
function [NI_absolute,NI_relative] = getNI(data)
if iscell(data.analysisDataST)
    for iElec = 1:size(data.analysisDataST{1},1)
        for k = 1:length(data.analysisDataST) 
            clear responseMatrix_elec
            responseMatrix_elec = squeeze(data.analysisDataST{k}(iElec,:,:));
            diff_responseMatrix_elec = squeeze(data.analysisDataST{k}(iElec,:,:))-squeeze(data.analysisDataBL{k}(iElec,:,:));
            NI_absolute{k}(iElec) = 2*responseMatrix_elec(1,5)/(responseMatrix_elec(1,1)+responseMatrix_elec(5,5))-1;
            NI_relative{k}(iElec) = 2*diff_responseMatrix_elec(1,5)/(diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))-1;
        end
    end
else
    for iElec = 1:size(data.analysisDataST,1)
        for iTF = 1:size(data.analysisDataST,2)
            clear responseMatrix_elec
            responseMatrix_elec = squeeze(data.analysisDataST(iElec,iTF,:,:));
            diff_responseMatrix_elec = squeeze(data.analysisDataST(iElec,iTF,:,:))-squeeze(data.analysisDataBL(iElec,iTF,:,:));
            NI_absolute(iElec,iTF) = 2*responseMatrix_elec(1,5)/(responseMatrix_elec(1,1)+responseMatrix_elec(5,5))-1;
            NI_relative(iElec,iTF) = 2*diff_responseMatrix_elec(1,5)/(diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))-1;
        end
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
function rescaleData(plotHandles,xMin,xMax,yLims)

[numRows,numCols] = size(plotHandles);
labelSize=14;
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

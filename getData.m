function[erpData,firingRateData,fftData,energyData,energyDataTF,oriTuningData,NI_Data,electrodeArray,N] = ...
    getData(folderSourceString,fileNameStringTMP,ElectrodeListTMP,dataParameters,tapers_MT,freqRanges,elecParams,removeERPFlag)

numDatasets = length(fileNameStringTMP);
dataset_init = 1;

while isempty(ElectrodeListTMP{dataset_init}{end})
    dataset_init= dataset_init+1;
end

disp(['Working on dataset ' num2str(dataset_init) 'of ' num2str(numDatasets)]);
[erpData,firingRateData,fftData,energyData,energyDataTF,oriTuningData,NI_Data,electrodeArray]...
    = getDataSingleSession(folderSourceString,fileNameStringTMP{dataset_init},...
    ElectrodeListTMP{dataset_init},dataParameters,tapers_MT,freqRanges,elecParams,removeERPFlag);
for j=1:2
    N{1,j} = squeeze(firingRateData.N(1,j,:,:));
end

if length(fileNameStringTMP)>1
    for i= dataset_init+1:numDatasets
        if isempty(ElectrodeListTMP{i}{end})
            continue
        end
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [erpDataTMP,firingRateDataTMP,fftDataTMP,energyDataTMP,energyDataTFTMP,oriTuningDataTMP,NI_DataTMP,electrodeArrayTMP] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},...
            ElectrodeListTMP{i},dataParameters,tapers_MT,freqRanges,elecParams,removeERPFlag);
        
        erpData.data = cat(1,erpData.data,erpDataTMP.data);
        erpData.analysisDataBL = cat(1,erpData.analysisDataBL,erpDataTMP.analysisDataBL);
        erpData.analysisData_cBL = cat(1,erpData.analysisData_cBL,erpDataTMP.analysisData_cBL);
        erpData.analysisDataST = cat(1,erpData.analysisDataST,erpDataTMP.analysisDataST);
        
        firingRateData.data = cat(1,firingRateData.data,firingRateDataTMP.data);
        firingRateData.spikeRasterData = cat(1,firingRateData.spikeRasterData,firingRateDataTMP.spikeRasterData);
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
        
        energyDataTF.data = cat(1,energyDataTF.data,energyDataTFTMP.data);
        energyDataTF.data_cBL = cat(1,energyDataTF.data_cBL,energyDataTFTMP.data_cBL);

        
        NI_Data.ERP_RMS_ST_Ray = cat(1, NI_Data.ERP_RMS_ST_Ray, NI_DataTMP.ERP_RMS_ST_Ray);
        NI_Data.ERP_dRMS_Ray = cat(1, NI_Data.ERP_dRMS_Ray, NI_DataTMP.ERP_dRMS_Ray);
        NI_Data.ERP_RMS_ST_Cohen = cat(1, NI_Data.ERP_RMS_ST_Cohen, NI_DataTMP.ERP_RMS_ST_Cohen);
        NI_Data.ERP_dRMS_Cohen = cat(1, NI_Data.ERP_dRMS_Cohen, NI_DataTMP.ERP_dRMS_Cohen);
        NI_Data.firingRate_ST_Ray = cat(1, NI_Data.firingRate_ST_Ray, NI_DataTMP.firingRate_ST_Ray);
        NI_Data.dfiringRate_Ray = cat(1, NI_Data.dfiringRate_Ray, NI_DataTMP.dfiringRate_Ray);
        NI_Data.firingRate_ST_Cohen = cat(1, NI_Data.firingRate_ST_Cohen, NI_DataTMP.firingRate_ST_Cohen);
        NI_Data.dfiringRate_Cohen = cat(1, NI_Data.dfiringRate_Cohen, NI_DataTMP.dfiringRate_Cohen);
        for j =1:length(NI_Data.fft_ST_Ray)
            NI_Data.fft_ST_Ray{j} = cat(2, NI_Data.fft_ST_Ray{j}, NI_DataTMP.fft_ST_Ray{j});
            NI_Data.dfft_Ray{j} = cat(2, NI_Data.dfft_Ray{j}, NI_DataTMP.dfft_Ray{j});
            NI_Data.fft_ST_Cohen{j} = cat(2, NI_Data.fft_ST_Cohen{j}, NI_DataTMP.fft_ST_Cohen{j});
            NI_Data.dfft_Cohen{j} = cat(2, NI_Data.dfft_Cohen{j}, NI_DataTMP.dfft_Cohen{j});
            NI_Data.energy_ST_Ray{j} = cat(2, NI_Data.energy_ST_Ray{j}, NI_DataTMP.energy_ST_Ray{j});
            NI_Data.denergy_Ray{j} = cat(2, NI_Data.denergy_Ray{j}, NI_DataTMP.denergy_Ray{j});
            NI_Data.energy_ST_Cohen{j} = cat(2, NI_Data.energy_ST_Cohen{j}, NI_DataTMP.energy_ST_Cohen{j});
            NI_Data.denergy_Cohen{j} = cat(2, NI_Data.denergy_Cohen{j}, NI_DataTMP.denergy_Cohen{j});
        end
        
        oriTuningData.PO = cat(2,oriTuningData.PO,oriTuningDataTMP.PO);
        oriTuningData.OS = cat(2,oriTuningData.OS,oriTuningDataTMP.OS);
        oriTuningData.FR = cat(1,oriTuningData.FR,oriTuningDataTMP.FR);
        oriTuningData.OriPairFR = cat(1,oriTuningData.OriPairFR,oriTuningDataTMP.OriPairFR);
        
        electrodeArray = cat(2,electrodeArray,electrodeArrayTMP);
        
        for j=1:2
            N{i,j} = squeeze(firingRateDataTMP.N(1,j,:,:)); 
        end
    end
end
N = N';
end

function [erpData,firingRateData,fftData,energyData,energyDataTF,oriTuningData,NI_Data,electrodeArray] = ...
    getDataSingleSession(folderSourceString,fileNameStringTMP,...
    ElectrodeListTMP,dataParameters,tapers_MT,freqRanges,elecParams,removeERPFlag)

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
folderExtract_oriTuning = fullfile(tuningProtocol_folderName,'extractedData');
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');
folderSave = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\savedDataV2');
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
    [computationVals,PO,OS] = savePrefOriAndOriSelectivitySpikes(monkeyName,expDate,oriTuning_protocolName,folderSourceString,gridType,dataParameters.stRange);
end

oriTuningData.PO = PO(ElectrodeListTMP{end});
oriTuningData.OS = OS(ElectrodeListTMP{end});
oriTuningData.FR = computationVals(ElectrodeListTMP{end},:);

[~,~,~,~,~,oValsUnique_Grating,~,~] = loadOriTuningParameterCombinations(folderExtract_oriTuning);
[~,~,~,~,~,~,oValsUnique_Plaid,~,~,~,~,~,~,oValsUnique2_Plaid,~,~] = loadParameterCombinations(folderExtract);

oriPairIndex = [find(oValsUnique_Plaid == oValsUnique_Grating),find(oValsUnique2_Plaid == oValsUnique_Grating)];
oriTuningData.OriPairFR = oriTuningData.FR(:,oriPairIndex);

if elecParams.oriSelectiveFlag
    fileToSave = fullfile(folderSave,[fileNameStringTMP ...
        '_N' num2str(elecParams.spikeCutoff) ...
        '_S' num2str(elecParams.snrCutoff)  '_oriTunedElecData_T_' ...
        num2str(1000*dataParameters.stRange(1))...
        '_' num2str(1000*dataParameters.stRange(2)) ...
        '_d'  num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2)) ...
        '_tapers'  num2str(tapers_MT(2)) ...
        '_removeERP' num2str(removeERPFlag) ...
        '_gse' num2str(elecParams.getSpikeElectrodesFlag)...
        '_Microelectrode_UnitID' num2str(elecParams.unitID) '.mat']);
else
    fileToSave = fullfile(folderSave,[fileNameStringTMP ...
        '_N' num2str(elecParams.spikeCutoff) ...
        '_S' num2str(elecParams.snrCutoff)  '_allElecData_T_' ...
        num2str(1000*dataParameters.stRange(1))...
        '_' num2str(1000*dataParameters.stRange(2)) ...
        '_d'  num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2)) ...
        '_tapers'  num2str(tapers_MT(2)) ...
        '_removeERP' num2str(removeERPFlag) ...
        '_gse' num2str(elecParams.getSpikeElectrodesFlag)...
        '_Microelectrode_UnitID' num2str(elecParams.unitID) '.mat']);
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
    
    % Set up movingWindow parameters for time-frequency plot
    winSize = 0.1;
    winStep = 0.025;
    movingwin = [winSize winStep];
    
    cListFlipped_Ori2 = flip(1:length(cValsUnique2)); % helps in plotting the responses from low to high contrast
    
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
                    clear energy_tf energyBL_tf
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
                                spikeRasterData{iElec,t,c_Ori2,c_Ori1} = spikeData(goodPos);
                                firingRatesBL(iElec,t,c_Ori2,c_Ori1) = mean(getSpikeCounts(spikeData(goodPos),dataParameters.blRange))/diff(dataParameters.blRange);
                                firingRatesST(iElec,t,c_Ori2,c_Ori1) = mean(getSpikeCounts(spikeData(goodPos),dataParameters.stRange))/diff(dataParameters.stRange);
                            end
                            
                            % fft data is processed for both static and
                            % flickering stimuli
                            if removeERPFlag == 0
                                fftBL = squeeze(mean(abs(fft(analogData(goodPos,blPos),[],2))));
                                fftST = squeeze(mean(abs(fft(analogData(goodPos,stPos),[],2))));
                            elseif removeERPFlag == 1
                                fftBL = squeeze(mean(abs(fft(removeERP(analogData(goodPos,blPos)),[],2))));
                                fftST = squeeze(mean(abs(fft(removeERP(analogData(goodPos,stPos)),[],2))));
                            end
                            fftDataBL(iElec,t,c_Ori2,c_Ori1,:) = conv2Log(fftBL); %#ok<*AGROW>
                            fftDataST(iElec,t,c_Ori2,c_Ori1,:) = conv2Log(fftST);
                            
                            % Power Estimation by MT method
                            % for both static and flickering stimuli
                            if removeERPFlag == 0
                                dataBL = analogData(goodPos,blPos)';
                                dataST = analogData(goodPos,stPos)';
                            elseif removeERPFlag == 1
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
                                    fftAmpBL{i}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(fftBL(:),freqVals,freqRanges{i}));
                                    fftAmpST{i}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(fftST(:),freqVals,freqRanges{i}));
                                    energyValsBL{i}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),freqValsMT,freqRanges{i}));
                                    energyValsST{i}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),freqValsMT,freqRanges{i}));
                                end
                                
                            elseif t == 2 %% fft and energy data for only SSVEP frequency
                                fftAmpBL{numFreqs}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(fftBL(:),freqVals,freqRanges{numFreqs}));
                                fftAmpST{numFreqs}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(fftST(:),freqVals,freqRanges{numFreqs}));
                                energyValsBL{numFreqs}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),freqValsMT,freqRanges{numFreqs}));
                                energyValsST{numFreqs}(iElec,c_Ori2,c_Ori1) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),freqValsMT,freqRanges{numFreqs}));
                            end
                            
                            
                            % computing time-frequency spectrum by multi-taper method
                            [tmpE_tf,tmpT_tf,freqVals_tf] = mtspecgramc(analogData(goodPos,:)',movingwin,params);
                            
                            timeVals_tf_= tmpT_tf + timeVals(1);
                            energy_tf = conv2Log(tmpE_tf)';
                            energyBL_tf = mean(energy_tf(:,timeVals_tf_>=dataParameters.blRange(1)& timeVals_tf_<=dataParameters.blRange(2)),2);
                            
                            mEnergy_tf(iElec,t,c_Ori2,c_Ori1,:,:) = energy_tf;
                            mEnergyBL_tf(iElec,t,c_Ori2,c_Ori1,:,:) = repmat(energyBL_tf,1,length(timeVals_tf_));

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
    firingRateData.spikeRasterData = spikeRasterData;
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
    
    % Time-Frequency data
    energyDataTF.data = mEnergy_tf;
    energyDataTF.data_cBL = getCommonBaseline(mEnergyBL_tf);
    energyDataTF.timeVals = timeVals_tf_;
    energyDataTF.freqVals = freqVals_tf;
    energyDataTF.N = N;
    
    % Segregation into Preferred-null axis is done when analysis is being
    % done for orientation selective electrodes
    if elecParams.oriSelectiveFlag
        elecs_neededtoFlipped = find(abs(oriTuningData.PO-oValsUnique2)<abs(oriTuningData.PO-oValsUnique));
        if ~isempty(elecs_neededtoFlipped)
            erpData = segregate_Pref_Null_data(erpData,elecs_neededtoFlipped);
            firingRateData = segregate_Pref_Null_data(firingRateData,elecs_neededtoFlipped);
            fftData = segregate_Pref_Null_data(fftData,elecs_neededtoFlipped);
            energyData = segregate_Pref_Null_data(energyData,elecs_neededtoFlipped);
        end
    end
    
    
    % Get Normalization Indices
    [NI_Data.ERP_RMS_ST_Ray, NI_Data.ERP_dRMS_Ray,NI_Data.ERP_RMS_ST_Cohen, NI_Data.ERP_dRMS_Cohen] = getNI(erpData);
    [NI_Data.firingRate_ST_Ray,NI_Data.dfiringRate_Ray,NI_Data.firingRate_ST_Cohen,NI_Data.dfiringRate_Cohen] = getNI(firingRateData);
    [NI_Data.fft_ST_Ray,NI_Data.dfft_Ray,NI_Data.fft_ST_Cohen,NI_Data.dfft_Cohen] = getNI(fftData);
    [NI_Data.energy_ST_Ray,NI_Data.denergy_Ray,NI_Data.energy_ST_Cohen,NI_Data.denergy_Cohen] = getNI(energyData);
    
    % Save Data for particular session
    save(fileToSave,'erpData','firingRateData','fftData','energyData','energyDataTF','oriTuningData','NI_Data','electrodeArray');
end
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

% Get Tuning protocol Parameter combinations
function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = ...
    loadOriTuningParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;            end

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

% Get Induced LFP data by subtracting trialaveraged ERP data from trialwise LFP Data
function Y = removeERP(X)
Y = X-repmat(mean(X,1),size(X,1),1);
end

% Get MeanEnergy for different frequency bands
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange)

posToAverage = intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2)));
eValue   = mean(mEnergy(posToAverage));
end

% get Common Baseline across all 5 (Ori 2) x 5 (Ori 1) contrast conditions
function data_BL = getCommonBaseline(data_BL)

if iscell(data_BL)
    size_data_BL = size(data_BL,2);
    num_con_Ori2 = size(data_BL{1},2);
    num_con_Ori1 = size(data_BL{1},3);
else
    size_data_BL = numel(size(data_BL));
    num_con_Ori2 = size(data_BL,3);
    num_con_Ori1 = size(data_BL,4);
end

if size_data_BL == 4 % baseline for analysis data (elec x Num_Contrast_Ori2 x Num_Contrast_Ori1 x (analysisDataBL))
    if iscell(data_BL)
        for iElec = 1:size(data_BL{1},1)
            for k = 1:length(data_BL)
                data_BL{k}(iElec,:,:) = repmat(mean(mean(squeeze(data_BL{k}(iElec,:,:)),2),1),[num_con_Ori2 num_con_Ori1]);
            end
        end
    else
        for iElec = 1:size(data_BL,1)
            for iTF = 1:size(data_BL,2)
                data_BL(iElec,iTF,:,:) = repmat(mean(mean(squeeze(data_BL(iElec,iTF,:,:)),2),1),[num_con_Ori2 num_con_Ori1]);
            end
        end
    end
elseif size_data_BL == 5 % baseline for timeSeries/PSD data (elec x TF x Num_Contrast_Ori2 x Num_Contrast_Ori1 x time/FreqVals x (dataBL))
    for iElec = 1:size(data_BL,1)
        for iTF = 1:size(data_BL,2)
            data_BL(iElec,iTF,:,:,:) = reshape(repmat(squeeze(mean(mean(squeeze(data_BL(iElec,iTF,:,:,:)),2),1))',[num_con_Ori2*num_con_Ori1 1]),[num_con_Ori2 num_con_Ori1 size(data_BL,5)]);
        end
    end
elseif size_data_BL == 6 % baseline for time-Frequency data (elec x TF x Num_Contrast_Ori2 x Num_Contrast_Ori1 x time/FreqVals x (dataBL_tf))
    for iElec = 1:size(data_BL,1)
        for iTF = 1:size(data_BL,2)
            data_BL(iElec,iTF,:,:,:,:) = repmat(mean(mean(squeeze(data_BL(iElec,iTF,:,:,:,:)),2),1),num_con_Ori2);
        end
    end    
end
end

% Segregate data into Preferred (positive x-axis -- increasing contrast) and Null (positive y axis -- increasing contrast)
function data = segregate_Pref_Null_data(data,elecs_neededtoFlipped)
for iElec = 1:length(elecs_neededtoFlipped)
    for iTF = 1:2
        disp (['elec: ' num2str(iElec), ', TF: ' num2str(iTF)])
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
    
    
    if ~iscell(data.analysisDataBL) && ~iscell(data.analysisDataST)
        data.analysisDataBL(elecs_neededtoFlipped(iElec),:,:,:) = flip(flip(permute(squeeze(data.analysisDataBL(elecs_neededtoFlipped(iElec),:,:,:)),[2 1]),1),2);
        data.analysisData_cBL(elecs_neededtoFlipped(iElec),:,:,:) = flip(flip(permute(squeeze(data.analysisData_cBL(elecs_neededtoFlipped(iElec),:,:,:)),[2 1]),1),2);
        data.analysisDataST(elecs_neededtoFlipped(iElec),:,:,:) = flip(flip(permute(squeeze(data.analysisDataST(elecs_neededtoFlipped(iElec),:,:,:)),[2 1]),1),2);
    elseif iscell(data.analysisDataBL) && iscell(data.analysisDataST)
        for iCell = 1:length(data.analysisDataBL)
            data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),:,:,:) = flip(flip(permute(squeeze(data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),:,:,:)),[2 1]),1),2);
            data.analysisData_cBL{iCell}(elecs_neededtoFlipped(iElec),:,:,:) = flip(flip(permute(squeeze(data.analysisData_cBL{iCell}(elecs_neededtoFlipped(iElec),:,:,:)),[2 1]),1),2);
            data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),:,:,:) = flip(flip(permute(squeeze(data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),:,:,:)),[2 1]),1),2);
        end
    end
end
end

% computing NI for different neural measures
function [NI_absolute_Ray,NI_relative_Ray,NI_absolute_Cohen,NI_relative_Cohen] = getNI(data)
if iscell(data.analysisDataST)
    for iElec = 1:size(data.analysisDataST{1},1)
        for k = 1:length(data.analysisDataST)
            clear responseMatrix_elec
            responseMatrix_elec = squeeze(data.analysisDataST{k}(iElec,:,:));
            diff_responseMatrix_elec = squeeze(data.analysisDataST{k}(iElec,:,:))-squeeze(data.analysisData_cBL{k}(iElec,:,:));
            NI_absolute_Ray{k}(iElec) = 2*responseMatrix_elec(1,5)/(responseMatrix_elec(1,1)+responseMatrix_elec(5,5))-1;
            NI_relative_Ray{k}(iElec) = 2*diff_responseMatrix_elec(1,5)/(diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))-1;
            NI_absolute_Cohen{k}(iElec) = (responseMatrix_elec(1,1)+ responseMatrix_elec(5,5))/responseMatrix_elec(1,5);
            NI_relative_Cohen{k}(iElec) = (diff_responseMatrix_elec(1,1)+ diff_responseMatrix_elec(5,5))/diff_responseMatrix_elec(1,5);
        end
    end
else
    for iElec = 1:size(data.analysisDataST,1)
        for iTF = 1:size(data.analysisDataST,2)
            clear responseMatrix_elec
            responseMatrix_elec = squeeze(data.analysisDataST(iElec,iTF,:,:));
            diff_responseMatrix_elec = squeeze(data.analysisDataST(iElec,iTF,:,:))-squeeze(data.analysisData_cBL(iElec,iTF,:,:));
            NI_absolute_Ray(iElec,iTF) = 2*responseMatrix_elec(1,5)/(responseMatrix_elec(1,1)+responseMatrix_elec(5,5))-1;
            NI_relative_Ray(iElec,iTF) = 2*diff_responseMatrix_elec(1,5)/(diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))-1;
            NI_absolute_Cohen(iElec,iTF) = (responseMatrix_elec(1,1)+responseMatrix_elec(5,5))/responseMatrix_elec(1,5);
            NI_relative_Cohen(iElec,iTF) = (diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))/diff_responseMatrix_elec(1,5);
        end
    end
end
end



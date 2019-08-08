% getNeuronQualityData returns a number of parameters characterizing the
% quality of isolation, snr, firing rate etc.


function saveNeuronSNRandRatesUnsortedPlaidNorm(monkeyName,stimulusPeriod,versionNum)

gridType = 'Microelectrode';
[expDates,protocolNames,~,~,datafolderSourceString] = dataInformationPlaidNorm(monkeyName,gridType);
numDays = length(expDates);

folderSourceString = strtok(datafolderSourceString,'\');
if versionNum == 1
    folderSave = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\snrAndRatesPlaidNorm');
elseif versionNum == 2
    folderSave = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\snrAndRatesPlaidNormV2');
end

if ~exist(folderSave,'dir')
    mkdir(folderSave);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:numDays
    expDate = expDates{i};
    protocolName = protocolNames{i};
    disp([expDate protocolName]);
    
    saveSpikeCountsAllElectrodes(monkeyName,expDate,protocolName,folderSourceString,gridType,folderSave,stimulusPeriod,versionNum);
    saveSNRFromSegments(monkeyName,expDate,protocolName,folderSourceString,folderSave);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveSNRFromSegments(monkeyName,expDate,protocolName,folderSourceString,folderSave)

folderSegment = [folderSourceString 'data\PlaidNorm\data\' monkeyName '\' 'Microelectrode\' expDate '\' protocolName '\segmentedData\Segments\'];
load([folderSegment 'segmentInfo.mat']);
segmentChannelsStored = sort(segmentChannelsStored);

disp('Computing SNR')
count=1;

for i=1:length(segmentChannelsStored)
    %disp(i);
    load([folderSegment 'elec' num2str(segmentChannelsStored(i))]);
    
    uniqueIDs = unique(unitID);
    for j=1:length(uniqueIDs)
        [snr(i),meanSpike(i,:),N(i),signal(i),noise(i)] = getSNR(segmentData(:,unitID==uniqueIDs(j))); %#ok<*NODEF>
        disp(['ElecNum: ' num2str(segmentChannelsStored(i)) ', SNR Count: ' num2str(length(snr))])
        segmentElectrodeList(i) = segmentChannelsStored(i);
        segmentUnitList(i)      =  uniqueIDs(j);  %#ok<*AGROW,*NASGU>
        %     count=count+1;
    end
    %     for j=1:length(uniqueIDs)
    %         [snr(count),meanSpike(count,:),N(count),signal(count),noise(count)] = getSNR(segmentData(:,unitID==uniqueIDs(j))); %#ok<*NODEF>
    %         disp(['ElecNum: ' num2str(segmentChannelsStored(i)) ', SNR Count: ' num2str(length(snr))])
    %         segmentElectrodeList(count) = segmentChannelsStored(i);
    %         segmentUnitList(count)      =  uniqueIDs(j);  %#ok<*AGROW,*NASGU>
    %         count=count+1;
    %     end
end

% Save
save([appendIfNotPresent(folderSave,'\') monkeyName expDate protocolName 'unsortedSNR.mat'] ,'snr','meanSpike','N','signal','noise','segmentElectrodeList','segmentUnitList');
end
function saveSpikeCountsAllElectrodes(monkeyName,expDate,protocolName,folderSourceString,gridType,folderSave,stimulusPeriod,versionNum)

validStimNum = 1;
baselinePeriod = -0.05+[-diff(stimulusPeriod) 0];

% Get parameter combinations
load(fullfile(folderSourceString,'data\PlaidNorm\data\',monkeyName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'));
% Get bad trials
badTrialFile = fullfile(folderSourceString,'data\PlaidNorm\data\',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat');
if exist(badTrialFile,'file')
    load(fullfile(folderSourceString,'data\PlaidNorm\data\',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'));
else
    [~,badTrials]= findBadTrialsWithLFPv4(monkeyName,expDate,protocolName,folderSourceString,gridType,[45,58,67,68,77],1,6,1000,-2000,[],1,[-0.5 0.5],1); % CheckTheseElectrodes alpa-[8 25 57 66 67] kesari - [45,58,67,68,77]
end
% Get Spike Data
load(fullfile(folderSourceString,'data\PlaidNorm\data\',monkeyName,gridType,expDate,protocolName,'segmentedData','Spikes','spikeInfo.mat'));
neuralChannelsStored = sort(neuralChannelsStored);

numElectrodes = length(neuralChannelsStored);
disp('Computing Spike Counts')

cListFlipped_Ori1 = flip(1:length(cValsUnique));
cListFlipped_Ori2 = flip(1:length(cValsUnique2));

for i=1:numElectrodes
    clear spikeData
    load(fullfile(folderSourceString,'data\PlaidNorm\data\',monkeyName,gridType,expDate,protocolName,'segmentedData','Spikes',['elec' num2str(neuralChannelsStored(i)) '_SID' num2str(SourceUnitID(i)) '.mat']));
    disp(['ElecNum: ' num2str(neuralChannelsStored(i))])
    
    if versionNum == 1
        
        for c_Ori1=1:length(cValsUnique)
            for c_Ori2=1:length(cValsUnique2)
                clear goodPos goodPos2
                goodPos = parameterCombinations{1,1,1,1,1,cListFlipped_Ori1(c_Ori1),1}; %#ok<USENS>
                goodPos2 = parameterCombinations2{1,1,1,1,1,c_Ori2,1}; %#ok<USENS>
                goodPos = intersect(goodPos,goodPos2);
                goodPos = setdiff(goodPos,badTrials);
                
                nStim(i,c_Ori1,c_Ori2) = sum(getSpikeCounts(spikeData(goodPos),stimulusPeriod));
                frStim(i,c_Ori1,c_Ori2) = (nStim(i,c_Ori1,c_Ori2)/length(goodPos))/diff(stimulusPeriod);
            end
        end
        
    elseif versionNum == 2
        
        for c_Ori2=1:length(cValsUnique2)
            for c_Ori1=1:length(cValsUnique)
                clear goodPos goodPos2
                goodPos = parameterCombinations{1,1,1,1,1,c_Ori1,1};
                goodPos2 = parameterCombinations2{1,1,1,1,1,cListFlipped_Ori2(c_Ori2),1};
                goodPos = intersect(goodPos,goodPos2);
                goodPos = setdiff(goodPos,badTrials);
                
                nStim(i,c_Ori2,c_Ori1) = sum(getSpikeCounts(spikeData(goodPos),stimulusPeriod));
                frStim(i,c_Ori2,c_Ori1) = (nStim(i,c_Ori2,c_Ori1)/length(goodPos))/diff(stimulusPeriod);
            end
        end
        
    end
    % For baseline, we combine all contrasts
    clear goodPos goodPos2
    goodPos = parameterCombinations{1,1,1,1,1,6,1};
    goodPos2 = parameterCombinations2{1,1,1,1,1,6,1};
    goodPos = intersect(goodPos,goodPos2);
    goodPos = setdiff(goodPos,badTrials);
    nBL(i)  = sum(getSpikeCounts(spikeData(goodPos),baselinePeriod));
    frBL(i) = (nBL(i)/length(goodPos))/diff(baselinePeriod);
    
end

save(fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSpikeCountsAndRatesPlaidNorm' num2str(round(1000*stimulusPeriod(1))) '_' num2str(round(1000*stimulusPeriod(2))) '.mat']),'nStim','frStim','nBL','frBL','neuralChannelsStored','SourceUnitID');

end

% getNeuronQualityData returns a number of parameters characterizing the
% quality of isolation, snr, firing rate etc.


function saveNeuronSNRandRatesUnsortedPlaidNorm(monkeyName,stimulusPeriod)

gridType = 'Microelectrode';
[expDates,protocolNames,~,folderSourceString] = dataInformationPlaidNorm(monkeyName,gridType);

numDays = length(expDates);
folderSave = 'E:\Projects\PlaidNormalization Project\snrAndRatesPlaidNorm';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:numDays
    expDate = expDates{i};
    protocolName = protocolNames{i};
    disp([expDate protocolName]);

    saveSpikeCountsAllElectrodes(monkeyName,expDate,protocolName,folderSourceString,gridType,folderSave,stimulusPeriod);
    saveSNRFromSegments(monkeyName,expDate,protocolName,folderSourceString,folderSave);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveSNRFromSegments(monkeyName,expDate,protocolName,folderSourceString,folderSave)

folderSegment = [folderSourceString 'data\' monkeyName '\' 'Microelectrode\' expDate '\' protocolName '\segmentedData\Segments\'];
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
function saveSpikeCountsAllElectrodes(monkeyName,expDate,protocolName,folderSourceString,gridType,folderSave,stimulusPeriod)

validStimNum = 1;
baselinePeriod = -0.05+[-diff(stimulusPeriod) 0];

% Get parameter combinations
load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'extractedData','parameterCombinations.mat'));
% Get bad trials
badTrialFile = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat');
if exist(badTrialFile,'file')
    load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','badTrials.mat'));
else
    [~,badTrials]= findBadTrialsWithLFPv4(monkeyName,expDate,protocolName,folderSourceString,gridType,[45,58,67,68,77],1,6,1000,-2000,[],1,[-0.5 0.5],1); % CheckTheseElectrodes alpa-[8 25 57 66 67] kesari - [45,58,67,68,77]
end
% Get Spike Data
load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','Spikes','spikeInfo.mat'));
neuralChannelsStored = sort(neuralChannelsStored);

numElectrodes = length(neuralChannelsStored);
nStim = zeros(numElectrodes,5,5);
frStim = zeros(numElectrodes,5,5);
nBL = zeros(1,numElectrodes);
frBL = zeros(1,numElectrodes);

disp('Computing Spike Counts')

for i=1:numElectrodes
    clear spikeData
    load(fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName,'segmentedData','Spikes',['elec' num2str(neuralChannelsStored(i)) '_SID' num2str(SourceUnitID(i)) '.mat']));
    disp(['ElecNum: ' num2str(neuralChannelsStored(i))])

    for cNum1=1:5
        for cNum2=1:5
            clear goodPos
            goodPos = parameterCombinations{1,1,1,1,1,cNum1,1}; %#ok<USENS>
            goodPos2 = parameterCombinations2{1,1,1,1,1,cNum2,1}; %#ok<USENS>
            goodPos = intersect(goodPos,goodPos2);
            goodPos = setdiff(goodPos,badTrials);

            nStim(i,cNum1,cNum2) = sum(getSpikeCounts(spikeData(goodPos),stimulusPeriod));
            frStim(i,cNum1,cNum2) = (nStim(i,cNum1,cNum2)/length(goodPos))/diff(stimulusPeriod);
        end
    end

    % For baseline, we combine all contrasts
    
        clear goodPos
        goodPos = parameterCombinations{1,1,1,1,1,6,1}; 
        goodPos2 = parameterCombinations2{1,1,1,1,1,6,1}; 
        goodPos = intersect(goodPos,goodPos2);
        goodPos = setdiff(goodPos,badTrials);
        nBL(i)  = sum(getSpikeCounts(spikeData(goodPos),baselinePeriod));
        frBL(i) = (nBL(i)/length(goodPos))/diff(baselinePeriod);
    
end

save(fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSpikeCountsAndRatesPlaidNorm' num2str(round(1000*stimulusPeriod(1))) '_' num2str(round(1000*stimulusPeriod(2))) '.mat']),'nStim','frStim','nBL','frBL','neuralChannelsStored','SourceUnitID');

end

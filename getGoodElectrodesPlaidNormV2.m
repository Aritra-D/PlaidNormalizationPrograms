% Electrode selection. Originally it was part of displayAttentionDataV2,
% but was made a separate program to get the list of good electrodes for
% energy computation.

function [allGoodElectrodes,allDaysToCombine,allGoodElectrodesStr,goodElectrodes,goodDays] = getGoodElectrodesPlaidNormV2(fileNameStringTMP,gridType,dRange,combineUniqueElectrodeData,getSpikeElectrodesFlag,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList)
if ~exist('gridType','var');                         gridType = 'microelectrode';           end
if ~exist('combineUniqueElectrodeData','var');       combineUniqueElectrodeData = 0;        end
if ~exist('getSpikeElectrodesFlag', 'var');          getSpikeElectrodesFlag = 1;            end
if ~exist('unitID', 'var');                          unitID = 0;                            end
if ~exist('spikeCutoff', 'var');                     spikeCutoff = 20;                      end
if ~exist('snrCutoff', 'var');                       snrCutoff = 2;                         end
if ~exist('timeRangeFRComputation', 'var');          timeRangeFRComputation = [0.15 .4];    end
if ~exist('contrastIndexList', 'var');               contrastIndexList = 1:5;               end
if ~exist('dRange', 'var');                          dRange = [0 0.75];                      end

% Getting Day & Protocol Info from fileNameString
if strcmp(fileNameStringTMP{1}(1:5),'alpaH') 
    monkeyName = 'alpaH'; 
    SelectedexpDate = fileNameStringTMP{1}(6:11); 
%     SelectedprotocolName = fileNameStringTMP{1}(12:end); 
    combineElecsAllDays = 0;
    
    if length(fileNameStringTMP) == 12
        combineElecsAllDays = 1;
    end
        
        
elseif strcmp(fileNameStringTMP{1}(1:7),'kesariH') 
    monkeyName = 'kesariH'; 
    SelectedexpDate = fileNameStringTMP{1}(8:13); 
%     SelectedprotocolName = fileNameStringTMP{1}(14:end); 
    combineElecsAllDays = 0;
    
    if length(fileNameStringTMP) == 10
        combineElecsAllDays = 1;
    end
    
elseif length(fileNameStringTMP) == 22
    
end


    
% Selection criteria should be compatible with Ray and Maunsell, 2010,
% Neuron and 2011, JNS, both of which used the attention data.

impedanceCutoff=2500;
allUsefulElectrodes  = [];
allDayIndices        = [];

OrientationTuningFlag = 0; % Plaid Protocols
[expDates,protocolNames,positionList,folderSourceString] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag);
numDays = length(expDates);

% Get a list of bad electrodes, Following Ray and Maunsell, 2011, JNS, a common list of bad electrodes was used

badElectrodeListAllDays = [];
for i=1:numDays
    expDate = expDates{i};

    % badElectrodes
    load(fullfile(folderSourceString, 'data',monkeyName,gridType,expDate, 'impedanceValues.mat'));
    badElectrodeListAllDays = cat(2,badElectrodeListAllDays,find(impedanceValues>impedanceCutoff));
end

badElectrodeList = unique(badElectrodeListAllDays);

% RF Info
load([monkeyName gridType 'RFData.mat']);
electrodeList = highRMSElectrodes(find(highRMSElectrodes<=81)); % Taking only the microelectrode set; ignoring highRMS ECoG elecs 81-90

if combineElecsAllDays == 0
    DayNum = find(strcmp(expDates,SelectedexpDate));
    
        a=positionList(1); e=positionList(2);

    clear d
    for j=1:length(electrodeList)
        azi = rfStats(electrodeList(j)).meanAzi;
        ele = rfStats(electrodeList(j)).meanEle;
        
        d(j) = sqrt(sum((azi-a)^2+(ele-e)^2)); %#ok<*AGROW>
    end  
    usefulElectrodes = setdiff(electrodeList(intersect(find(d>=dRange(1)),find(d<dRange(2)))),badElectrodeList);
        % If spikes are required, check which electrodes have enough spikes
    if getSpikeElectrodesFlag
        usefulElectrodes = intersect(usefulElectrodes,getRateAndSNRInfo(monkeyName,expDates{DayNum},protocolNames{DayNum},unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList));
    end
    
    
    if ~isempty(usefulElectrodes)
        disp(['day' num2str(DayNum) ': electrodes - ' num2str(usefulElectrodes)]);
    else
        disp('No Good Electrode is found for this day! Please Choose another Session!')
    end
    
    goodElectrodes = usefulElectrodes;
    goodDays       = DayNum;
    
    allGoodElectrodes = goodElectrodes;
    allGoodElectrodesStr = '';
    
    for i=1:length(allGoodElectrodes) % For each unique electrode
        thisElectrode = allGoodElectrodes(i);
        allDaysToCombine = goodDays;
        allGoodElectrodesStr = [allGoodElectrodesStr 'elec' num2str(thisElectrode) '_' expDates{goodDays} '_' protocolNames{goodDays} '|'];
    end
    
    allGoodElectrodesStr = [allGoodElectrodesStr 'all'];

elseif combineElecsAllDays == 1

    
    
    
    for i=1:numDays
        a=positionList(1); e=positionList(2);

        clear d
        for j=1:length(electrodeList)
            azi = rfStats(electrodeList(j)).meanAzi;
            ele = rfStats(electrodeList(j)).meanEle;

            d(j) = sqrt(sum((azi-a)^2+(ele-e)^2)); %#ok<*AGROW>
        end

        usefulElectrodes = setdiff(electrodeList(intersect(find(d>=dRange(1)),find(d<dRange(2)))),badElectrodeList);

        % If spikes are required, check which electrodes have enough spikes
        if getSpikeElectrodesFlag
            usefulElectrodes = intersect(usefulElectrodes,getRateAndSNRInfo(monkeyName,expDates{i},protocolNames{i},unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList));
        end

        if ~isempty(usefulElectrodes)

            disp(['day' num2str(i) ': electrodes - ' num2str(usefulElectrodes)]);
            allUsefulElectrodes = cat(2,allUsefulElectrodes,usefulElectrodes);
            allDayIndices       = cat(2,allDayIndices,i+zeros(1,length(usefulElectrodes)));
        end
    end


    goodElectrodes = allUsefulElectrodes;
    goodDays       = allDayIndices;

    disp([num2str(length(goodElectrodes)) ' good electrodes, ' num2str(length(unique(goodElectrodes))) ' unique.']);



    clear allGoodElectrodes allDaysToCombine allGoodElectrodesStr 
    if combineUniqueElectrodeData
        allUniqueGoodElectrodes  = unique(goodElectrodes);

        allGoodElectrodesStr = '';

        for i=1:length(allUniqueGoodElectrodes) % For each unique electrode
            thisElectrode = allUniqueGoodElectrodes(i);

            daysToAverage = goodDays(goodElectrodes==thisElectrode);

            allDaysToCombine{i} = daysToAverage;
            allGoodElectrodesStr = [allGoodElectrodesStr 'elec' num2str(thisElectrode) ', N=' num2str(length(daysToAverage)) '|'];
        end
        allGoodElectrodesStr = [allGoodElectrodesStr 'all'];
        allGoodElectrodes    = allUniqueGoodElectrodes;

    else
        allGoodElectrodes = goodElectrodes;
        allGoodElectrodesStr = '';

        for i=1:length(allGoodElectrodes) % For each unique electrode
            thisElectrode = allGoodElectrodes(i);
            allDaysToCombine{i} = goodDays(i);
            allGoodElectrodesStr = [allGoodElectrodesStr 'elec' num2str(thisElectrode) '_' expDates{goodDays(i)} '_' protocolNames{goodDays(i)} '|'];
        end
        allGoodElectrodesStr = [allGoodElectrodesStr 'all'];
    end
end

end
function goodElectrodes = getRateAndSNRInfo(monkeyName,expDate,protocolName,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList)

% The travelingWavesProject required spike electrodes, so we stored the SNR
% and number of spikes for the unattended condition between 0.2 to 0.4
% seconds in folderSave. See saveNeuronSNRandRatesSorted in
% travelingWavesProject for more information

% Good electrodes: >spikeCutoff spikes between 150-400 ms & SNR>snrCutoff
% Returns all electrodes that meet these criteria

% Modification - We allow different intervals for computation of firing
% rates.
folderSave = 'E:\Projects\PlaidNormalization Project\snrAndRatesPlaidNorm';

fileSaveFR = fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSpikeCountsAndRatesPlaidNorm' num2str(round(1000*timeRangeFRComputation(1))) '_' num2str(round(1000*timeRangeFRComputation(2))) '.mat']);
if exist(fileSaveFR,'file')
    x=load(fileSaveFR);
else
end
y=load(fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSNR.mat']));

goodPos = (x.SourceUnitID==unitID)&(min(min(x.nStim(:,contrastIndexList,contrastIndexList),[],2),[],3)'>spikeCutoff)&((y.snr>snrCutoff)==1);
goodElectrodes = x.neuralChannelsStored(goodPos);
end
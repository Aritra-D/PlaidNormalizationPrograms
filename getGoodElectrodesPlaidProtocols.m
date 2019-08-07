% Electrode selection. This displays all the good Electrodes for selected
% criteria for Plaid Normalization Project

function [allGoodElectrodes,allDaysToCombine,allGoodElectrodesStr,goodElectrodes,goodDays] = getGoodElectrodesPlaidProtocols(monkeyName,versionNum,gridType,dRange,combineUniqueElectrodeData,getSpikeElectrodesFlag,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList)

% if ~exist('folderSourceString','var');               folderSourceString = 'E:\';            end
if ~exist('gridType','var');                         gridType = 'microelectrode';           end
if ~exist('combineUniqueElectrodeData','var');       combineUniqueElectrodeData = 0;        end
if ~exist('getSpikeElectrodesFlag', 'var');          getSpikeElectrodesFlag = 1;            end
if ~exist('unitID', 'var');                          unitID = 0;                            end
if ~exist('spikeCutoff', 'var');                     spikeCutoff = 20;                      end
if ~exist('snrCutoff', 'var');                       snrCutoff = 2;                         end
if ~exist('timeRangeFRComputation', 'var');          timeRangeFRComputation = [0.15 0.4];    end
if ~exist('contrastIndexList', 'var');               contrastIndexList = {[1,1],[5,5]};     end
if ~exist('dRange', 'var');                          dRange = [0 0.75];                     end

% Selection criteria should be compatible with Ray and Maunsell, 2010,
% Neuron and 2011, JNS, both of which used the attention data.

impedanceCutoff=2500;
allUsefulElectrodes  = [];
allDayIndices        = [];

[expDates,protocolNames,positionList,~,dataFolderSourceString] = dataInformationPlaidNorm(monkeyName,gridType,0); % OrientationTuningFlag set to zero

numDays = length(expDates);

% Get a list of bad electrodes, Following Ray and Maunsell, 2011, JNS, a common list of bad electrodes was used

badElectrodeListAllDays = [];
for i=1:numDays
    expDate = expDates{i};
    
    % badElectrodes
    load(fullfile(dataFolderSourceString,'data',monkeyName,gridType,expDate, 'impedanceValues.mat'));
    badElectrodeListAllDays = cat(2,badElectrodeListAllDays,find(impedanceValues>impedanceCutoff));
end

badElectrodeList = unique(badElectrodeListAllDays);

% RF Info
load([monkeyName gridType 'RFData.mat']);
electrodeList = highRMSElectrodes(find(highRMSElectrodes<=81));

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
        usefulElectrodes = intersect(usefulElectrodes,getRateAndSNRInfo(monkeyName,expDates{i},protocolNames{i},dataFolderSourceString,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList,versionNum));
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
function goodElectrodes = getRateAndSNRInfo(monkeyName,expDate,protocolName,dataFolderSourceString,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList,versionNum)

% Good electrodes: >spikeCutoff spikes between 150-400 ms & SNR>snrCutoff
% Returns all electrodes that meet these criteria

% Modification - We allow different intervals for computation of firing
% rates.

ProjectFolderSourceString = strtok(dataFolderSourceString,'\');
if versionNum == 1
    folderSave = fullfile(ProjectFolderSourceString,'Projects\Aritra_PlaidNormalizationProject\snrAndRatesPlaidNorm\');
elseif versionNum == 2
    folderSave = fullfile(ProjectFolderSourceString,'Projects\Aritra_PlaidNormalizationProject\snrAndRatesPlaidNormV2\');
end

fileSaveFR = fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSpikeCountsAndRatesPlaidNorm' num2str(round(1000*timeRangeFRComputation(1))) '_' num2str(round(1000*timeRangeFRComputation(2))) '.mat']);
if exist(fileSaveFR,'file')
    x=load(fileSaveFR);
else
end
y=load(fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSNR.mat']));

goodPos1 = (x.SourceUnitID==unitID)&...
    ((x.nStim(:,contrastIndexList{1}(1),contrastIndexList{1}(2))'>spikeCutoff)|...
    (x.nStim(:,contrastIndexList{2}(1),contrastIndexList{2}(2))'>spikeCutoff))&...
    ((y.snr>snrCutoff)==1);

goodPos2 = (x.SourceUnitID==unitID)&...
    ((x.frStim(:,1,1)'+ x.frStim(:,5,5)')<= 2.*x.frStim(:,1,5)')&...
    (y.snr>snrCutoff==1);

goodPos = intersect(goodPos1,goodPos2);

goodElectrodes = x.neuralChannelsStored(goodPos);
end
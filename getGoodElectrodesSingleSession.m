% Electrode selection. Originally it was part of displayAttentionDataV2,
% but was made a separate program to get the list of good electrodes for
% energy computation.

function [allGoodElectrodesArray,allGoodElectrodesStr,goodElectrodes] = getGoodElectrodesSingleSession(monkeyName,gridType,SessionNum,folderSourceString,dRange,getSpikeElectrodesFlag,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList)

if ~exist('gridType','var');                         gridType = 'microelectrode';           end
if ~exist('getSpikeElectrodesFlag', 'var');          getSpikeElectrodesFlag = 1;            end
if ~exist('unitID', 'var');                          unitID = 0;                            end
if ~exist('spikeCutoff', 'var');                     spikeCutoff = 20;                      end
if ~exist('snrCutoff', 'var');                       snrCutoff = 2;                         end
if ~exist('timeRangeFRComputation', 'var');          timeRangeFRComputation = [0.15 .4];    end
if ~exist('contrastIndexList', 'var');               contrastIndexList = 1:5;               end
if ~exist('dRange', 'var');                          dRange = [0 0.75];                     end

    
% Selection criteria should be compatible with Ray and Maunsell, 2010,
% Neuron and 2011, JNS, both of which used the attention data.

impedanceCutoff=2500;
OrientationTuningFlag = 0; % Plaid Protocols
[expDates,protocolNames,positionList] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag);
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

% Session Selection
clear expDate
expDate = expDates{SessionNum};
protocolName = protocolNames{SessionNum};


% RF Info
load([monkeyName gridType 'RFData.mat']);
electrodeList = highRMSElectrodes(find(highRMSElectrodes<=81)); % Taking only the microelectrode set; ignoring highRMS ECoG elecs 81-90

% Stimulus Center Position 
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
    usefulElectrodes = intersect(usefulElectrodes,getRateAndSNRInfo(monkeyName,expDate,protocolName,folderSourceString,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList));
end

if ~isempty(usefulElectrodes)
    disp(['day' num2str(SessionNum) ': electrodes - ' num2str(usefulElectrodes)]);
else
    disp(['day' num2str(SessionNum) ': electrodes - No Electrodes' ]);
end

goodElectrodes = usefulElectrodes;
pos = 1;

for iElec =1:length(goodElectrodes)
    goodElectrodesListArray{pos} = goodElectrodes(iElec);
    pos = pos+1;
end
goodElectrodesListArray{pos} = goodElectrodes;

allGoodElectrodesArray = goodElectrodesListArray;
allGoodElectrodesStr = '';

for i=1:length(goodElectrodes) % For each unique electrode
    thisElectrode = goodElectrodes(i);
    allGoodElectrodesStr = [allGoodElectrodesStr 'elec' num2str(thisElectrode) '_' expDate '_' protocolName '|'];
end

if ~isempty(goodElectrodes)
allGoodElectrodesStr = [allGoodElectrodesStr 'all'];
else
    allGoodElectrodesStr = 'No Electrode Found!';
end

end
function goodElectrodes = getRateAndSNRInfo(monkeyName,expDate,protocolName,folderSourceString,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList)

% Good electrodes: >spikeCutoff spikes between 150-400 ms & SNR>snrCutoff
% Returns all electrodes that meet these criteria
% Modification - We allow different intervals for computation of firing
% rates.

folderSave = fullfile(folderSourceString,'Projects\PlaidNormalizationProject\snrAndRatesPlaidNorm');%'E:\Projects\PlaidNormalizationProject\snrAndRatesPlaidNorm';
fileSaveFR = fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSpikeCountsAndRatesPlaidNorm' num2str(round(1000*timeRangeFRComputation(1))) '_' num2str(round(1000*timeRangeFRComputation(2))) '.mat']);
if exist(fileSaveFR,'file')
    x=load(fileSaveFR);
else
end

y=load(fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSNR.mat']));

goodPos = (x.SourceUnitID==unitID)&(min(min(x.nStim(:,contrastIndexList,contrastIndexList),[],2),[],3)'>spikeCutoff)&((y.snr>snrCutoff)==1);
goodElectrodes = x.neuralChannelsStored(goodPos);
end
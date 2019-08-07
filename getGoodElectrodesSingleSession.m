% Electrode selection. Originally it was part of displayAttentionDataV2,
% but was made a separate program to get the list of good electrodes for
% energy computation.

function [allGoodElectrodesStr,allGoodElectrodesArray,goodElectrodes] ...
    = getGoodElectrodesSingleSession(monkeyName,expDate,protocolName,...
    gridType,elecParams,timeRangeFRComputation,folderSourceString,versionNum,...
    contrastIndexList)

if ~isempty(elecParams)
    spikeCutoff = elecParams.spikeCutoff;
    snrCutoff = elecParams.snrCutoff;
    dRange = elecParams.dRange;
    getSpikeElectrodesFlag = elecParams.getSpikeElectrodesFlag;
    unitID = elecParams.unitID;
    getOriSelectiveFlag = elecParams.oriSelectiveFlag;
else
    if ~exist('spikeCutoff', 'var');                     spikeCutoff = 20;                      end
    if ~exist('snrCutoff', 'var');                       snrCutoff = 2;                       end
    if ~exist('dRange', 'var');                          dRange = [0 0.75];                     end
    if ~exist('getSpikeElectrodesFlag', 'var');          getSpikeElectrodesFlag = 1;            end
    if ~exist('unitID', 'var');                          unitID = 0;                            end
    if ~exist('getOriSelectiveFlag', 'var');             getOriSelectiveFlag = 0;               end
end

if ~exist('gridType','var');                         gridType = 'microelectrode';           end
if ~exist('contrastIndexList', 'var');               contrastIndexList = {[1,1],[5,5]};     end

% Selection criteria is selected as per Ray and Maunsell, 2010,
% Neuron and 2011, JNS, both of which used the attention data.
impedanceCutoff=2500;
[expDates,~,positionList,oriList,~] = dataInformationPlaidNorm(monkeyName,gridType,0); % Plaid Protocols
numDays = length(expDates);

% Get a list of bad electrodes, Following Ray and Maunsell, 2011, JNS, a common list of bad electrodes was used

badElectrodeListAllDays = [];
for i=1:numDays
    exp_Date = expDates{i};
    
    % badElectrodes
    load(fullfile(folderSourceString, 'data',monkeyName,gridType,exp_Date, 'impedanceValues.mat'));
    badElectrodeListAllDays = cat(2,badElectrodeListAllDays,find(impedanceValues>impedanceCutoff));
end

badElectrodeList = unique(badElectrodeListAllDays);

% RF Info
load([monkeyName gridType 'RFData.mat']);
electrodeList = highRMSElectrodes(find(highRMSElectrodes<=81)); % Taking only the microelectrode set; ignoring highRMS ECoG elecs 81-90 if any

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
    usefulElectrodes = intersect(usefulElectrodes,getRateAndSNRInfo(monkeyName,expDate,protocolName,folderSourceString,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList,versionNum));
end

SessionNum = find(strcmp(expDate,expDates));

if getOriSelectiveFlag
    oriTuning_protocolName = ['GRF_00' num2str(str2double(protocolName(5:end))-1)]; % The protocol Number is just the immediate precedent of the main protocol
    folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,oriTuning_protocolName);
    folderSave = fullfile(folderName,'savedData');
    if ~exist(folderSave,'dir')
        mkdir(folderSave);
    end
    fileToSave = fullfile(folderSave,['oriTuningData_' num2str(1000*timeRangeFRComputation(1)) 'ms_' num2str(1000*timeRangeFRComputation(2)) 'ms.mat']);
    if exist(fileToSave,'file')
        load(fileToSave)
    else
        [~,PO, ~] = savePrefOriAndOriSelectivitySpikes(monkeyName,expDate,oriTuning_protocolName,folderSourceString,gridType);
    end
    oriCutOff = 15;
    if oriList{SessionNum}(1) == 0
        oriSelectiveElectrodes = find((PO<=oriList{SessionNum}(1)+oriCutOff & PO>=oriList{SessionNum}(1)+180-oriCutOff) | ...
            (PO<=oriList{SessionNum}(2)+oriCutOff & PO>=oriList{SessionNum}(2)-oriCutOff));
    else
        oriSelectiveElectrodes = find((PO<=oriList{SessionNum}(1)+oriCutOff & PO>=oriList{SessionNum}(1)-oriCutOff ) | ...
            (PO<=oriList{SessionNum}(2)+oriCutOff & PO>=oriList{SessionNum}(2)-oriCutOff));
    end
    
    usefulElectrodes = intersect(usefulElectrodes,oriSelectiveElectrodes);
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
function goodElectrodes = getRateAndSNRInfo(monkeyName,expDate,protocolName,folderSourceString,unitID,spikeCutoff,snrCutoff,timeRangeFRComputation,contrastIndexList,versionNum)

% Good electrodes: spikes>spikeCutoff between 150-400 ms for individual gratings at 50% contrast & SNR>snrCutoff
% Returns all electrodes that meet these criteria
% Modification - We allow different intervals for computation of firing
% rates.
if versionNum == 1
    folderSave = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\snrAndRatesPlaidNorm');
elseif versionNum == 2
    folderSave = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\snrAndRatesPlaidNormV2');
end
fileSaveFR = fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSpikeCountsAndRatesPlaidNorm' num2str(round(1000*timeRangeFRComputation(1))) '_' num2str(round(1000*timeRangeFRComputation(2))) '.mat']);
if exist(fileSaveFR,'file')
    x=load(fileSaveFR);
else
end

y=load(fullfile(folderSave,[monkeyName expDate protocolName 'unsortedSNR.mat']));

goodPos = (x.SourceUnitID==unitID)&...
    ((x.nStim(:,contrastIndexList{1}(1),contrastIndexList{1}(2))'>spikeCutoff)|...
    (x.nStim(:,contrastIndexList{2}(1),contrastIndexList{2}(2))'>spikeCutoff))&...
    (y.snr>snrCutoff==1);

% goodPos2 = (x.SourceUnitID==unitID)&...
%           ((x.frStim(:,1,1)'+ x.frStim(:,5,5)')<2*x.frStim(:,1,5)')&...
%           (y.snr>snrCutoff==1);

% goodPos = intersect(goodPos1,goodPos2);
goodElectrodes = x.neuralChannelsStored(goodPos);
end
function [compVals,PO,OS] = savePrefOriAndOriSelectivitySpikesAndLFP(monkeyName,expDate,protocolName,folderSourceString,gridType,electrodeList,timeForComputation)

dbstop if error
folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderSpikes = fullfile(folderSegment,'Spikes');
folderLFP = fullfile(folderSegment,'LFP');
folderBadTrials = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\badTrials');%fullfile(folderSegment,'badTrials.mat');

folderSave = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\oriTuningData');
if ~exist(folderSave,'dir')
    mkdir(folderSave);
end

% % load Spike Information
% [neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes); %#ok<ASGLU>
%
% load LFP Information
[~,timeVals,~,~] = loadlfpInfo(folderLFP); %#ok<*STOUT>

% freqRanges{1} = [8 12]; % alpha
freqRanges{1} = [32 80]; % gamma
freqRanges{2} = [104 248]; % hi-gamma


% Get Parameter Combinations
[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);
% Defining the orientation values; Lablib rounds of the 0.5 deg value in the oValsUnique.
% Including them back to use in legends and axis labels.
for i=1:length(oValsUnique)
    if mod(i,2)==0
        oValsUnique(i)=oValsUnique(i)+0.5;
    else
    end
end
aLen = length(aValsUnique);
eLen = length(eValsUnique);
sLen = length(sValsUnique);
fLen = length(fValsUnique);
oLen = length(oValsUnique);
cLen = length(cValsUnique);
tLen = length(tValsUnique);

unitID = 0;

% Get bad trials
badTrialFile = fullfile(folderBadTrials,[monkeyName expDate protocolName '_badTrials.mat']);
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    disp(['Loading' badTrialFile])
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end


% Set Up Timing Params for MT
Fs = round(1/(timeVals(2)-timeVals(1)));
rangePos = round(diff(timeForComputation)*Fs);
stPos = find(timeVals>=timeForComputation(1),1)+ (1:rangePos);
% Set up params for MT
tapers_MT = [1 1];
params.tapers   = tapers_MT;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 250];
params.trialave = 1;

% PO = zeros(1,size(electrodeList,2));
% OS = zeros(1,size(electrodeList,2));
spkData = zeros(length(electrodeList),oLen);
% mEnergyVsFreqST = zeros(length(electrodeList),oLen);
energyValsST = cell(1,length(freqRanges));
compVals = cell(1,length(freqRanges)+1);

% Electrode Loop
for iElec = 1:length(electrodeList)
    if iElec ==1
        disp(['monkeyName: ' monkeyName, ', expDate:' expDate, ', protocolName: ' protocolName])
    end
    % Get Spike Data
    clear signal spikeData
    spkFile = fullfile(folderSpikes,['elec' num2str(electrodeList(iElec)) '_SID' num2str(unitID) '.mat']);
    disp(['elec: ' num2str(iElec) '--' num2str(electrodeList(iElec))])

    if exist(spkFile,'file')
        load(spkFile); %#ok<*LOAD>
    else
        error('Spike data not found')
    end
    
    % Get LFP Data
    clear analogData
    lfpFile = fullfile(folderLFP,['elec' num2str(electrodeList(iElec)) '.mat']);
    if exist(lfpFile,'file')
        load(lfpFile); %#ok<*LOAD>
    else
        error('LFP data not found')
    end
    
    % Main Loop
    for iOri = 1:oLen
        clear goodPos
        goodPos = parameterCombinations{aLen,eLen,sLen,fLen,iOri,cLen,tLen};
        goodPos = setdiff(goodPos,badTrials);
        
        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            % compute spikes
            clear firingRate
            firingRate = mean(getSpikeCounts(spikeData(goodPos),timeForComputation))/diff(timeForComputation);  
            spkData(iElec,iOri) = firingRate;
            
            % compute LFP
            clear dataST
            dataST = analogData(goodPos,stPos)';
            [tmpEST,freqValsST] = mtspectrumc(dataST,params);
%             mEnergyVsFreqST(iElec,iOri) = energyST;
            
            for i=1:length(freqRanges)
                energyValsST{i}(iElec,iOri) = getMeanEnergyForAnalysis(tmpEST(:),freqValsST,freqRanges{i});
            end
            
        end
    end
    
    compVals{1} = spkData;
    compVals{2} = energyValsST{1};
    compVals{3} = energyValsST{2};
    
    [PO{1}(iElec),OS{1}(iElec)] = getOrientationTuning(spkData(iElec,:),oValsUnique); %#ok<*AGROW>
    
    for i=1:length(freqRanges)
        [PO{i+1}(iElec),OS{i+1}(iElec)] = getOrientationTuning(energyValsST{i}(iElec,:),oValsUnique);
    end
    
%     PO.Spikes = PO_Spikes;
%     PO.Energy = PO_Energy;
%     OS.Spikes = OS_Spikes;
%     OS.Energy = OS_Energy;
    
end


fileToSave = fullfile(folderSave,[monkeyName,'_',gridType,'_',expDate,'_',protocolName,'_oriTuningData_' num2str(1000*timeForComputation(1)) 'ms_' num2str(1000*timeForComputation(2)) 'ms.mat']);

save(fileToSave,'compVals','PO','OS');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Accessory Functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These functions are mostly curated from displaySingleChannelGRF under
% Programs\ProgramsMAP\Common Programs\display

% Load LFP Info
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored); %#ok<NODEF>
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

% load Data
% function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
% fileName = fullfile(folderSpikes,'spikeInfo.mat');
% if exist(fileName,'file')
%     load(fileName);
%     [neuralChannelsStored,I]=sort(neuralChannelsStored); %#ok<NODEF>
%     SourceUnitID=SourceUnitID(I); %#ok<NODEF>
% else
%     neuralChannelsStored=[];
%     SourceUnitID=[];
% end
% end

function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;             end
if ~exist('cValsUnique','var');    cValsUnique=100;                     end
if ~exist('tValsUnique','var');    tValsUnique=0;                       end
end
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end

% Get MeanEnergy for across frequency bins of a desired frequency band
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange)

posToAverage = intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2)));
eValue   = mean(mEnergy(posToAverage));
end

function [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique)
num=0;
den=0;

for j=1:length(oValsUnique)
    num = num+computationVals(j)*sind(2*oValsUnique(j));
    den = den+computationVals(j)*cosd(2*oValsUnique(j));
end

prefOrientation = round(90*atan2(num,den)/pi);
orientationSelectivity = abs(den+1i*num)/sum(computationVals);

if prefOrientation<0
    prefOrientation = prefOrientation+180;
end
end


function [computationVals,PO,OS] = getPrefOriAndOriSelectivitySpikes(subjectName,expDate,protocolName,folderSourceString,gridType)

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderSpikes = fullfile(folderSegment,'Spikes');

% load Spike Information
[neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes); %#ok<ASGLU>

% Get Parameter Combinations
[parameterCombinations,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract);

aLen = length(aValsUnique);
eLen = length(eValsUnique);
sLen = length(sValsUnique);
fLen = length(fValsUnique);
oLen = length(oValsUnique);
cLen = length(cValsUnique);
tLen = length(tValsUnique);

% Signal Timing Parameters
% signalRange = [-0.1 0.5];
% fftRange = [0 250];
% blRange = [-0.2 0];
% stRange = [0.2 0.4];
% For orientation selectivity analysis
timeForComputation = [150 400]/1000; % s

% neuralChannelString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs);
unitID = 0;
% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end
    
    PO = zeros(1,size(neuralChannelsStored,2));
    OS = zeros(1,size(neuralChannelsStored,2));
    computationVals = zeros(length(neuralChannelsStored),oLen);
    % Electrode Loop
    for iElec = 1:length(neuralChannelsStored)
    % Get Spike Data
    clear signal spikeData
    fileName = fullfile(folderSpikes,['elec' num2str(iElec) '_SID' num2str(unitID) '.mat']);
        if exist(fileName,'file')
        load(fileName);
        else
            continue
        end

        % Main Loop

        for iOri = 1:oLen
            clear goodPos
            goodPos = parameterCombinations{aLen,eLen,sLen,fLen,iOri,cLen,tLen};
            goodPos = setdiff(goodPos,badTrials);
                
            if isempty(goodPos)
                disp('No entries for this combination..')
            else
                disp(['pos=' num2str(iOri) ',n=' num2str(length(goodPos))]);
                % compute spikes
                clear firingRate
                firingRate = mean(getSpikeCounts(spikeData(goodPos),timeForComputation))/diff(timeForComputation);
                computationVals(iElec,iOri) = firingRate;
            end
        end
        if iElec == 1 || iElec == length(neuralChannelsStored)
        disp(['Orientation selectivity values calculated between ' num2str(timeForComputation(1)) '-' num2str(timeForComputation(2))  ' s']);
        end
%         disp(['orientation values: ' num2str(computationVals)]);
        [PO(iElec),OS(iElec)] = getOrientationTuning(computationVals(iElec,:),oValsUnique);
%         disp(['elec: ' num2str(iElec) ', prefOri: ' num2str(round(PO(iElec))) ', sel: ' num2str(OS(iElec))]);
    end
end

% load Data
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = fullfile(folderSpikes,'spikeInfo.mat');
if exist(fileName,'file')
    load(fileName);
    [neuralChannelsStored,I]=sort(neuralChannelsStored); %#ok<NODEF>
    SourceUnitID=SourceUnitID(I); %#ok<NODEF>
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end

function [parameterCombinations,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique] = loadParameterCombinations(folderExtract) %#ok<STOUT>

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;             end
if ~exist('cValsUnique','var');    cValsUnique=100;                     end
if ~exist('tValsUnique','var');    tValsUnique=0;                       end
end
function badTrials = loadBadTrials(badTrialFile) %#ok<STOUT>
load(badTrialFile);
end

% function outString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs)
% outString='';
% for i=1:length(neuralChannelsStored)
%     outString = cat(2,outString,[num2str(neuralChannelsStored(i)) ', SID ' num2str(SourceUnitIDs(i)) '|']);
% end 
% end

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


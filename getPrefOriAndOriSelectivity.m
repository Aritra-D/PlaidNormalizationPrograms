function [PO,OS] = getPrefOriAndOriSelectivity(subjectName,expDate,protocolName,folderSourceString,gridType)

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
% folderSpikes = fullfile(folderSegment,'Spikes');

% load LFP Information
[analogChannelsStored,timeVals,~,analogInputNums] = loadlfpInfo(folderLFP);
% [neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);

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
blRange = [-0.2 0];
stRange = [0.2 0.4];
% For orientation selectivity analysis
timeForComputation = [40 100]/1000; % s

[~,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);


% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end
    
    PO = zeros(1,size(analogChannelsStored,2));
    OS = zeros(1,size(analogChannelsStored,2));
    
    % Electrode Loop
    for iElec = 1:length(analogChannelsStored)
    % Get LFP Data
    clear signal analogData
    load(fullfile(folderLFP,analogChannelStringArray{iElec}));
    

        % Main Loop
        computationVals = zeros(1,oLen);
        for iOri = 1:oLen
            clear goodPos
            goodPos = parameterCombinations{aLen,eLen,sLen,fLen,iOri,cLen,tLen};
            goodPos = setdiff(goodPos,badTrials);
                
            if isempty(goodPos)
                disp('No entries for this combination..')
            else
%                 disp(['pos=' num2str(iOri) ',n=' num2str(length(goodPos))]);

                Fs = round(1/(timeVals(2)-timeVals(1)));
                if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
                    disp('baseline and stimulus ranges are not the same');
                else
%                     range = blRange;
%                     rangePos = round(diff(range)*Fs);
%                     blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
%                     stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
%                     xs = 0:1/diff(range):Fs-1/diff(range);
                end

                xsComputation = intersect(find(timeVals>=timeForComputation(1)),find(timeVals<timeForComputation(2)));
%                 freqComputation = intersect(find(xs>=freqForComputation(1)),find(xs<=freqForComputation(2)));

                % compute ERP
                clear erp
                erp = mean(analogData(goodPos,:),1);
                computationVals(iOri) = abs(min(erp(xsComputation)));
                
                
            end
        end
        if iElec == 1 || iElec == length(analogChannelsStored)
        disp(['Orientation selectivity values calculated between ' num2str(timeForComputation(1)) '-' num2str(timeForComputation(2))  ' s']);
        end
%         disp(['orientation values: ' num2str(computationVals)]);
        [PO(iElec),OS(iElec)] = getOrientationTuning(computationVals,oValsUnique);
        disp(['elec: ' num2str(iElec) ', prefOri: ' num2str(round(PO(iElec))) ', sel: ' num2str(OS(iElec))]);
    end
end

% load Data
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored);
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

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

function [outString,outArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums)
outString='';
count=1;
for i=1:length(analogChannelsStored)
    outArray{count} = ['elec' num2str(analogChannelsStored(i))]; %#ok<AGROW>
    outString = cat(2,outString,[outArray{count} '|']);
    count=count+1;
end
if ~isempty(analogInputNums)
    for i=1:length(analogInputNums)
        outArray{count} = ['ainp' num2str(analogInputNums(i))];
        outString = cat(2,outString,[outArray{count} '|']);
        count=count+1;
    end
end
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


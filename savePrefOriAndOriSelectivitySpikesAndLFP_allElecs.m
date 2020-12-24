function [compVals,PO,OS] = savePrefOriAndOriSelectivitySpikesAndLFP_allElecs(monkeyName,gridType,folderSourceString,timeForComputation)

dbstop if error

OrientationTuningFlag = 1;
[expDates_oriTuning, protocolNames_oriTuning,~,~,~] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag);

% freqRanges{1} = [8 12]; % alpha
freqRanges{1} = [32 80]; % gamma
freqRanges{2} = [104 248]; % hi-gamma
    
for iDay = 1:length(expDates_oriTuning)
    clear PO OS computation Vals 
    
    folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDates_oriTuning{iDay},protocolNames_oriTuning{iDay});
   
    % Get folders
    folderExtract = fullfile(folderName,'extractedData');
    folderSegment = fullfile(folderName,'segmentedData');
    folderSpikes = fullfile(folderSegment,'Spikes');
    folderLFP = fullfile(folderSegment,'LFP');
    folderBadTrials = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\badTrials');%fullfile(folderSegment,'badTrials.mat');
    
    
    folderSave = fullfile(strtok(folderSourceString,'\'),'Projects\Aritra_PlaidNormalizationProject\oriTuningData_allElecs');
    if ~exist(folderSave,'dir')
        mkdir(folderSave);
    end
    
    % % load Spike Information
    [neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes); %#ok<ASGLU>
    
    % load LFP Information
    [analogChannelsStored,timeVals,~,~] = loadlfpInfo(folderLFP); %#ok<*STOUT>
    
    electrodeList{1} = neuralChannelsStored; electrodeList{2} = analogChannelsStored;
    

    
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
    badTrialFile = fullfile(folderBadTrials,[monkeyName expDates_oriTuning{iDay} protocolNames_oriTuning{iDay} '_badTrials.mat']);
    if ~exist(badTrialFile,'file')
        disp('Bad trial file does not exist...');
        badTrials=[];
    else
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
    
    
    spkData = zeros(length(electrodeList{1}),oLen);
    energyValsST = cell(1,length(freqRanges));
    compVals = cell(1,length(freqRanges)+1);
    
    % Electrode Loop
    for iElecType=1:length(electrodeList)
        for iElec = 1:length(electrodeList{iElecType})
            if iElec ==1
                disp(['monkeyName: ' monkeyName, ', expDate:' expDates_oriTuning{iDay}, ', protocolName: '  protocolNames_oriTuning{iDay}])
            end
            
            if iElecType==1
                % Get Spike Data
                clear signal spikeData
                spkFile = fullfile(folderSpikes,['elec' num2str(neuralChannelsStored(iElec)) '_SID' num2str(unitID) '.mat']);
%                 disp(['elec: ' num2str(iElec) '--' num2str(neuralChannelsStored(iElec))])
                
                if exist(spkFile,'file')
                    load(spkFile); %#ok<*LOAD>
                else
                    error('Spike data not found')
                end
                
            elseif iElecType == 2
                % Get LFP Data
                clear analogData
                lfpFile = fullfile(folderLFP,['elec' num2str(analogChannelsStored(iElec)) '.mat']);
                if exist(lfpFile,'file')
                    load(lfpFile); %#ok<*LOAD>
                else
                    error('LFP data not found')
                end
            end
            
            % Main Loop
            for iOri = 1:oLen
                clear goodPos
                goodPos = parameterCombinations{aLen,eLen,sLen,fLen,iOri,cLen,tLen};
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    disp('No entries for this combination..')
                else
                    if iElecType==1
                        % compute spikes
                        clear firingRate
                        firingRate = mean(getSpikeCounts(spikeData(goodPos),timeForComputation))/diff(timeForComputation);
                        spkData(iElec,iOri) = firingRate;
                        
                    elseif iElecType == 2
                        % compute LFP
                        clear dataST
                        dataST = analogData(goodPos,stPos)';
                        [tmpEST,freqValsST] = mtspectrumc(dataST,params);
                        for i=1:length(freqRanges)
                            energyValsST{i}(iElec,iOri) = getMeanEnergyForAnalysis(tmpEST(:),freqValsST,freqRanges{i});
                        end
                    end
                end
            end
            
            if iElecType ==1
            compVals{1} = spkData;
            [PO{1}(iElec),OS{1}(iElec)] = getOrientationTuning(spkData(iElec,:),oValsUnique); %#ok<*AGROW>

            elseif iElecType ==2
                compVals{2} = energyValsST{1};
                compVals{3} = energyValsST{2};
                for i=1:length(freqRanges)
                    [PO{i+1}(iElec),OS{i+1}(iElec)] = getOrientationTuning(energyValsST{i}(iElec,:),oValsUnique);
                end
            end
        end
    end
    
    fileToSave = fullfile(folderSave,[monkeyName,'_',gridType,'_',expDates_oriTuning{iDay},'_', protocolNames_oriTuning{iDay},'_oriTuningData_' num2str(1000*timeForComputation(1)) 'ms_' num2str(1000*timeForComputation(2)) 'ms_GammaRange' num2str(freqRanges{1}(1)) '_' num2str(freqRanges{1}(2)) '_' 'Hz.mat']);
    save(fileToSave,'compVals','PO','OS','neuralChannelsStored','analogChannelsStored');
end


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

% load Spike Info
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


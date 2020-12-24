function combineProtocols_EyeData(monkeyName,expDate,protocolNames,gridType,folderSourceString,skipStims)
% Combine two (adjacent) protocols recorded on the same day, having the same parameters.

if ~exist('skipStims','var'); skipStims = 0; end

folderName1 = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolNames{1});
folderName2 = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolNames{2});

%% Do some sanity checking %%

% check for the same parameter combinations
pc1 = load(fullfile(folderName1,'extractedData','parameterCombinations.mat'));
pc2 = load(fullfile(folderName2,'extractedData','parameterCombinations.mat'));
if ~isequal(size(pc1.parameterCombinations),size(pc2.parameterCombinations)) || ...
        ~isequal(pc1.aValsUnique,pc2.aValsUnique) || ~isequal(pc1.eValsUnique,pc2.eValsUnique) || ...
        ~isequal(pc1.sValsUnique,pc2.sValsUnique) || ~isequal(pc1.fValsUnique,pc2.fValsUnique) || ...
        ~isequal(pc1.oValsUnique,pc2.oValsUnique) || ~isequal(pc1.cValsUnique,pc2.cValsUnique) || ...
        ~isequal(pc1.tValsUnique,pc2.tValsUnique)
    error('Cannot combine protocols, parameter mismatch!');
end

% check for the same parameter combinations (DTF/Plaids)
if isfield(pc1,'parameterCombinations2') && isfield(pc2,'parameterCombinations2')
    if ~isequal(size(pc1.parameterCombinations2),size(pc2.parameterCombinations2)) || ...
            ~isequal(pc1.aValsUnique2,pc2.aValsUnique2) || ~isequal(pc1.eValsUnique2,pc2.eValsUnique2) || ...
            ~isequal(pc1.sValsUnique2,pc2.sValsUnique2) || ~isequal(pc1.fValsUnique2,pc2.fValsUnique2) || ...
            ~isequal(pc1.oValsUnique2,pc2.oValsUnique2) || ~isequal(pc1.cValsUnique2,pc2.cValsUnique2) || ...
            ~isequal(pc1.tValsUnique2,pc2.tValsUnique2)
        error('Cannot combine protocols, parameter mismatch!');
    end
elseif (isfield(pc1,'parameterCombinations2') && ~isfield(pc2,'parameterCombinations2')) || ...
        (~isfield(pc1,'parameterCombinations2') && isfield(pc2,'parameterCombinations2'))
    error('Cannot combine protocols, parameter mismatch!');
end

% make sure the same LFP electrodes are available in both protocols
lfpElec1 = load(fullfile(folderName1,'segmentedData','LFP','lfpInfo.mat'));
lfpElec2 = load(fullfile(folderName2,'segmentedData','LFP','lfpInfo.mat'));
if ~isequal(lfpElec1.analogInputNums,lfpElec2.analogInputNums) || ...
        ~isequal(lfpElec1.timeVals,lfpElec2.timeVals) || ...
        ~isequal(lfpElec1.electrodesStored,lfpElec2.electrodesStored) || ...
        ~isequal(lfpElec1.analogChannelsStored,lfpElec2.analogChannelsStored)
    error('Cannot combine protocols, LFP info mismatch!');
end

% make sure the same Spiking electrodes are available in both protocols
spkElec1 = load(fullfile(folderName1,'segmentedData','Spikes','spikeInfo.mat'));
spkElec2 = load(fullfile(folderName2,'segmentedData','Spikes','spikeInfo.mat'));
if ~isequal(spkElec1.neuralChannelsStored,spkElec2.neuralChannelsStored) || ...
        ~isequal(spkElec1.SourceUnitID,spkElec2.SourceUnitID)
    error('Cannot combine protocols, Spike info mismatch!');
end

% make sure the same Segment electrodes are available in both protocols
sgmtElec1 = load(fullfile(folderName1,'segmentedData','Segments','segmentInfo.mat'));
sgmtElec2 = load(fullfile(folderName2,'segmentedData','Segments','segmentInfo.mat'));
if ~isequal(sgmtElec1.segmentChannelsStored,sgmtElec2.segmentChannelsStored)
    error('Cannot combine protocols, Segment info mismatch!');
end

% everything ok, get offset for trial numbers in second set
st1 = load(fullfile(folderName1,'extractedData','goodStimNums.mat'));
numStims1 = length(st1.goodStimNums)-skipStims;
stimsToSkip = st1.goodStimNums(numStims1+1:end);

% create new directory to save combined data
folderNameCombined = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,strjoin(protocolNames,'_'));
makeDirectory(folderNameCombined);
makeDirectory(fullfile(folderNameCombined,'extractedData'));
makeDirectory(fullfile(folderNameCombined,'segmentedData','eyeData'));
% makeDirectory(fullfile(folderNameCombined,'segmentedData','LFP'));
% makeDirectory(fullfile(folderNameCombined,'segmentedData','Spikes'));

%% Parameter Combinations %%

disp('Combining parameters');
for aPos = 1:size(pc1.parameterCombinations,1)
    for ePos = 1:size(pc1.parameterCombinations,2)
        for sPos = 1:size(pc1.parameterCombinations,3)
            for fPos = 1:size(pc1.parameterCombinations,4)
                for oPos = 1:size(pc1.parameterCombinations,5)
                    for cPos = 1:size(pc1.parameterCombinations,6)
                        for tPos = 1:size(pc1.parameterCombinations,7)
                            pc1.parameterCombinations{aPos,ePos,sPos,fPos,oPos,cPos,tPos} = union(setdiff(pc1.parameterCombinations{aPos,ePos,sPos,fPos,oPos,cPos,tPos},stimsToSkip), ...
                                pc2.parameterCombinations{aPos,ePos,sPos,fPos,oPos,cPos,tPos}+numStims1);
                        end
                    end
                end
            end
        end
    end
end
save(fullfile(folderNameCombined,'extractedData','parameterCombinations'),'-struct','pc1');
if isfield(pc1,'parameterCombinations2') && isfield(pc2,'parameterCombinations2')
    for aPos = 1:size(pc1.parameterCombinations2,1)
        for ePos = 1:size(pc1.parameterCombinations2,2)
            for sPos = 1:size(pc1.parameterCombinations2,3)
                for fPos = 1:size(pc1.parameterCombinations2,4)
                    for oPos = 1:size(pc1.parameterCombinations2,5)
                        for cPos = 1:size(pc1.parameterCombinations2,6)
                            for tPos = 1:size(pc1.parameterCombinations2,7)
                                pc1.parameterCombinations2{aPos,ePos,sPos,fPos,oPos,cPos,tPos} = union(setdiff(pc1.parameterCombinations2{aPos,ePos,sPos,fPos,oPos,cPos,tPos},stimsToSkip), ...
                                    pc2.parameterCombinations2{aPos,ePos,sPos,fPos,oPos,cPos,tPos}+numStims1);
                            end
                        end
                    end
                end
            end
        end
    end
    save(fullfile(folderNameCombined,'extractedData','parameterCombinations'),'-struct','pc1'); % overwrite
end

%% Eye Data

% Extracted Eye Data
disp('Combining Extracted Eye Data')

eyeData1 = load(fullfile(folderName1,'extractedData','EyeData.mat'));
eyeData2 = load(fullfile(folderName2,'extractedData','EyeData.mat'));

if ~isequal(eyeData1.eyeRangeMS,eyeData2.eyeRangeMS)
    error('Error combining Eye Data: Eye data Time Vals do not match')
end

eyeDataC.eyeRangeMS = eyeData1.eyeRangeMS;
eyeDataC.eyeData = cat(2,eyeData1.eyeData(1:numStims1),eyeData2.eyeData);
save(fullfile(folderNameCombined,'extractedData','EyeData'),'-struct','eyeDataC');


% Segmented Eye Data
disp('Combining Segmented Eye Data')

clear eyeData1 eyeData2 eyeDataC
eyeData1 = load(fullfile(folderName1,'segmentedData','eyeData','eyeDataDeg.mat'));
eyeData2 = load(fullfile(folderName2,'segmentedData','eyeData','eyeDataDeg.mat'));

eyeDataC.eyeDataDegX = cat(1,eyeData1.eyeDataDegX,eyeData2.eyeDataDegX);
eyeDataC.eyeDataDegY = cat(1,eyeData1.eyeDataDegY,eyeData2.eyeDataDegY);

save(fullfile(folderNameCombined,'segmentedData','eyeData','eyeDataDeg'),'-struct','eyeDataC');

clear eyeData1 eyeData2 eyeDataC
eyeData1 = load(fullfile(folderName1,'segmentedData','eyeData','EyeDataStimPos.mat'));
eyeData2 = load(fullfile(folderName2,'segmentedData','eyeData','EyeDataStimPos.mat'));

eyeDataC.durationsMS = eyeData1.durationsMS;
eyeDataC.eyeXAllPos = cellfun(@(x,y) cat(1,x,y),eyeData1.eyeXAllPos(1,:),eyeData2.eyeXAllPos(1,:),'UniformOutput',false);
eyeDataC.eyeYAllPos = cellfun(@(x,y) cat(1,x,y),eyeData1.eyeYAllPos(1,:),eyeData2.eyeYAllPos(1,:),'UniformOutput',false);
eyeDataC.xs = eyeData1.xs;
save(fullfile(folderNameCombined,'segmentedData','eyeData','EyeDataStimPos'),'-struct','eyeDataC');


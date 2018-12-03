function combineProtocols(monkeyName,expDate,protocolNames,gridType,folderSourceString,skipStims,combineSegments)
% Combine two (adjacent) protocols recorded on the same day, having the same parameters.

if ~exist('skipStims','var'); skipStims = 0; end
if ~exist('combineSegments','var'); combineSegments = 0; end

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
makeDirectory(fullfile(folderNameCombined,'segmentedData','LFP'));
makeDirectory(fullfile(folderNameCombined,'segmentedData','Spikes'));
    
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

%% LFP data %%

disp('Combining LFP');
for i = 1:length(lfpElec1.analogChannelsStored)
    elecName = ['elec' num2str(lfpElec1.analogChannelsStored(i))];
    disp(elecName);
    lfpData1 = load(fullfile(folderName1,'segmentedData','LFP',elecName));
    lfpData2 = load(fullfile(folderName2,'segmentedData','LFP',elecName));
    if ~isequal(lfpData1.analogInfo,lfpData2.analogInfo)
        error(['Error combining LFP: ' elecName]);
    end
    clear lfpDataC;
    lfpDataC.analogInfo = lfpData1.analogInfo;
    lfpDataC.analogData = cat(1,lfpData1.analogData(1:numStims1,:),lfpData2.analogData); 
    save(fullfile(folderNameCombined,'segmentedData','LFP',elecName),'-struct','lfpDataC');
end
save(fullfile(folderNameCombined,'segmentedData','LFP','lfpInfo'),'-struct','lfpElec1');

%% Spiking data %%

disp('Combining Spikes');
for i = 1:length(spkElec1.neuralChannelsStored)
    elecName = ['elec' num2str(spkElec1.neuralChannelsStored(i)) '_SID' num2str(spkElec1.SourceUnitID(i))];
    disp(elecName);
    spkData1 = load(fullfile(folderName1,'segmentedData','Spikes',elecName));
    spkData2 = load(fullfile(folderName2,'segmentedData','Spikes',elecName));
    if ~isequal(spkData1.neuralInfo,spkData2.neuralInfo)
        error(['Error combining Spikes: ' elecName]);
    end
    clear spkDataC;
    spkDataC.neuralInfo = spkData1.neuralInfo;
    spkDataC.spikeData = [spkData1.spikeData(1:numStims1) spkData2.spikeData];
    if isfield(spkData1,'spikeID')
        % increment spike IDs in 2nd protocol by spike count of 1st protocol
        spkDataC.spikeID = [spkData1.spikeID cellfun(@(x)x+sgmtElec1.numItems(i),spkData2.spikeID,'UniformOutput',false)];
    end
    save(fullfile(folderNameCombined,'segmentedData','Spikes',elecName),'-struct','spkDataC');
end
save(fullfile(folderNameCombined,'segmentedData','Spikes','spikeInfo'),'-struct','spkElec1');

%% Segment data (for SNR calculations) if requested %%
% NOTE: 'skipStims' argument applies only to extracted data and will not be
% used in the segment data combination below. This means that combined
% segment data may have some segments which are not in spike data above.

if combineSegments
    makeDirectory(fullfile(folderNameCombined,'segmentedData','Segments'));
    
    % find differences in start times for the two protocols to reproduce
    % approximate timestamps for combined data
    fi1 = load(fullfile(folderName1,'extractedData','NEVFileInfo.mat'));
    fi2 = load(fullfile(folderName2,'extractedData','NEVFileInfo.mat'));
    dv1 = [fi1.fileInfo.Time_Year fi1.fileInfo.Time_Month fi1.fileInfo.Time_Day ...
        fi1.fileInfo.Time_Hour fi1.fileInfo.Time_Min fi1.fileInfo.Time_Sec+fi1.fileInfo.Time_MilliSec/1000];
    dv2 = [fi2.fileInfo.Time_Year fi2.fileInfo.Time_Month fi2.fileInfo.Time_Day ...
        fi2.fileInfo.Time_Hour fi2.fileInfo.Time_Min fi2.fileInfo.Time_Sec+fi2.fileInfo.Time_MilliSec/1000];
    timeDiff = etime(dv2,dv1);

    disp('Combining Segments');
    for i = 1:length(sgmtElec1.segmentChannelsStored)
        elecName = ['elec' num2str(sgmtElec1.segmentChannelsStored(i))];
        disp(elecName);
        sgmtData1 = load(fullfile(folderName1,'segmentedData','Segments',elecName));
        sgmtData2 = load(fullfile(folderName2,'segmentedData','Segments',elecName));
        if ~isequal(sgmtData1.segmentInfo,sgmtData2.segmentInfo)
            error(['Error combining Segments: ' elecName]);
        end
        clear sgmtDataC;
        % increment timestamps in 2nd protocol by time elapsed since 1st protocol
        sgmtDataC.timeStamp = [sgmtData1.timeStamp; sgmtData2.timeStamp+timeDiff];
        sgmtDataC.segmentInfo = sgmtData1.segmentInfo;
        if isequal(size(sgmtData1.segmentData),[1 48])
            sgmtDataC.segmentData = cat(2,sgmtData1.segmentData',sgmtData2.segmentData);
        elseif isequal(size(sgmtData2.segmentData),[1 48])
            sgmtDataC.segmentData = cat(2,sgmtData1.segmentData,sgmtData2.segmentData');
        else
            sgmtDataC.segmentData = cat(2,sgmtData1.segmentData,sgmtData2.segmentData);
        end
        sgmtDataC.sampleCount = [sgmtData1.sampleCount; sgmtData2.sampleCount];
        sgmtDataC.unitID = [sgmtData1.unitID; sgmtData2.unitID];
        save(fullfile(folderNameCombined,'segmentedData','Segments',elecName),'-struct','sgmtDataC');
    end
    save(fullfile(folderNameCombined,'segmentedData','Segments','segmentInfo'),'-struct','sgmtElec1');
end

end

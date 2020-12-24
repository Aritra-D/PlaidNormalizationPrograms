function [eyeDataDegX,eyeDataDegY,eyeData,eyeSpeedData,microsaccadeData,timeValsEyeData,trialNums,FsEye] = getEyeDataIndividualMonkey_v2(monkeyName,cleanDataFolder,timeRange,cutOff)

gridType = 'Microelectrode';
FsEye = 200; % Hz % FsEye is constant across protocols.

[expDates, protocolNames,~,~,~] = dataInformationPlaidNorm(monkeyName);

for iProt = 1:length(expDates)
    disp(['Processing Eye Data: Monkey: ' monkeyName, ', expDate: ' expDates{iProt}, ', protocolName: ' protocolNames{iProt}])
    folderExtract = fullfile(cleanDataFolder,'data',monkeyName,gridType,expDates{iProt},protocolNames{iProt},'extractedData');
    folderSegment = fullfile(cleanDataFolder,'data',monkeyName,gridType,expDates{iProt},protocolNames{iProt},'segmentedData');
    
    % Get Combinations
    [parameterCombinations,parameterCombinations2,...
        aValsUnique,eValsUnique,~,~,oValsUnique,cValsUnique,tValsUnique, ...
        aValsUnique2,eValsUnique2,~,~,oValsUnique2,cValsUnique2,tValsUnique2] = loadParameterCombinations(folderExtract); %#ok<*ASGLU>
    
    cListFlipped_Ori2 = flip(1:length(cValsUnique2)); % helps in plotting the responses from low to high contrast
    
    if aValsUnique ~= aValsUnique2 || eValsUnique ~= eValsUnique2
        error('Azimuths and/or elevations do not match!');
    end
    
    if tValsUnique ~= tValsUnique2
        error('TF do not match!');
    else
        tList = 1:length(tValsUnique);
    end
    
    a=1; e=1; s=1; f=1; o=1;
    eyeData_Extracted = load(fullfile(folderExtract,'EyeData'));
    timeRange_EyeData =  eyeData_Extracted.eyeRangeMS;
    timeValsEyeData = (timeRange_EyeData(1):1000/FsEye:timeRange_EyeData(2)-1000/FsEye)/1000;
    
    clear eyeData_TMP
    eyeData_TMP = load(fullfile(folderSegment,'eyeData','eyeDataDeg'));
    
    % Get bad trials
    badTrialFile = fullfile(folderSegment,'badTrials.mat');
    if ~exist(badTrialFile,'file')
        disp('Bad trial file does not exist...');
        badTrials=[];
    else
        badTrials = loadBadTrials(badTrialFile);
        disp([num2str(length(badTrials)) ' bad trials']);
    end
    
    for t = 1:length(tList)
        for c_Ori2 = 1:length(cValsUnique2)
            for c_Ori1 = 1:length(cValsUnique)
                
                clear goodPos
                goodPos = parameterCombinations{a,e,s,f,o,c_Ori1,tList(t)};
                goodPos2 = parameterCombinations2{a,e,s,f,o,cListFlipped_Ori2(c_Ori2),tList(t)};
                goodPos = intersect(goodPos,goodPos2);
                goodPos = setdiff(goodPos,badTrials);
                
                if isempty(goodPos)
                    disp('No entries for this combination..');
                else
                    trialNums{iProt}(c_Ori2,c_Ori1,t) = length(goodPos);
                    if t==1
                    eyeData.eyeDataDegX_static{iProt}{c_Ori2,c_Ori1} = eyeData_TMP.eyeDataDegX(goodPos,:);
                    eyeData.eyeDataDegY_static{iProt}{c_Ori2,c_Ori1} = eyeData_TMP.eyeDataDegY(goodPos,:);
                    elseif t==2
                    eyeData.eyeDataDegX_CP{iProt}{c_Ori2,c_Ori1} = eyeData_TMP.eyeDataDegX(goodPos,:);
                    eyeData.eyeDataDegY_CP{iProt}{c_Ori2,c_Ori1} = eyeData_TMP.eyeDataDegY(goodPos,:);
                    end
                end
            end
        end
    end
end


eyeDataDegX_static = eyeData.eyeDataDegX_static{1,1};
eyeDataDegY_static = eyeData.eyeDataDegY_static{1,1};
eyeDataDegX_CP = eyeData.eyeDataDegX_CP{1,1};
eyeDataDegY_CP = eyeData.eyeDataDegY_CP{1,1};

for c_Ori2 = 1:length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        for iProt = 2:length(expDates)
            eyeDataDegX_static_TMP = eyeData.eyeDataDegX_static{1,iProt};
            eyeDataDegY_static_TMP = eyeData.eyeDataDegY_static{1,iProt};
            eyeDataDegX_CP_TMP = eyeData.eyeDataDegX_CP{1,iProt};
            eyeDataDegY_CP_TMP = eyeData.eyeDataDegY_CP{1,iProt};

            eyeDataDegX_static{c_Ori2,c_Ori1} = cat(1,eyeDataDegX_static{c_Ori2,c_Ori1},eyeDataDegX_static_TMP{c_Ori2,c_Ori1});
            eyeDataDegY_static{c_Ori2,c_Ori1} = cat(1,eyeDataDegY_static{c_Ori2,c_Ori1},eyeDataDegY_static_TMP{c_Ori2,c_Ori1});
            eyeDataDegX_CP{c_Ori2,c_Ori1} = cat(1,eyeDataDegX_CP{c_Ori2,c_Ori1},eyeDataDegX_CP_TMP{c_Ori2,c_Ori1});
            eyeDataDegY_CP{c_Ori2,c_Ori1} = cat(1,eyeDataDegY_CP{c_Ori2,c_Ori1},eyeDataDegY_CP_TMP{c_Ori2,c_Ori1});
        end
    end
end

eyeDataDegX.static = eyeDataDegX_static;
eyeDataDegX.CP = eyeDataDegX_CP;

eyeDataDegY.static = eyeDataDegY_static;
eyeDataDegY.CP = eyeDataDegY_CP;

% eyeSpeed calculated over data combined across sessions
eyeDataDeg{1} = eyeDataDegX.static;
eyeDataDeg{2} = eyeDataDegY.static;
eyeDataDeg{3} = eyeDataDegX.CP;
eyeDataDeg{4} = eyeDataDegY.CP;

eyeSpeed_Data = cell(1,length(eyeDataDeg));
for iData = 1:length(eyeDataDeg)
    for c_Ori2 = 1:length(cValsUnique2)
        for c_Ori1 = 1:length(cValsUnique)
            clear data
            data = eyeDataDeg{iData}{c_Ori2,c_Ori1};
            lengthEyeSignal = size(data,2);
            for j=1:size(data,1)
                eyeSpeed_Data{iData}{c_Ori2,c_Ori1}(j,:) = [data(j,2:lengthEyeSignal)-data(j,1:lengthEyeSignal-1) 0];
            end
        end
    end
end

eyeSpeedData.StaticX = eyeSpeed_Data{1};
eyeSpeedData.StaticY = eyeSpeed_Data{2};
eyeSpeedData.CPX = eyeSpeed_Data{3};
eyeSpeedData.CPY = eyeSpeed_Data{4};

for c_Ori2 = 1:length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        clear eyeSpeedX_static eyeSpeedY_static eyeSpeedX_CP eyeSpeedY_CP
        eyeSpeedX_static = eyeSpeed_Data{1}{c_Ori2,c_Ori1}; eyeSpeedY_static = eyeSpeed_Data{2}{c_Ori2,c_Ori1};
        eyeSpeedX_CP = eyeSpeed_Data{3}{c_Ori2,c_Ori1}; eyeSpeedY_CP = eyeSpeed_Data{4}{c_Ori2,c_Ori1};
        eyeSpeedData.StaticMag{c_Ori2,c_Ori1} = FsEye*sqrt(eyeSpeedX_static.^2 + eyeSpeedY_static.^2);
        eyeSpeedData.CPMag{c_Ori2,c_Ori1} = FsEye*sqrt(eyeSpeedX_CP.^2 + eyeSpeedY_CP.^2);
    end
end

% cutoff = 15;
% timeRange = [-0.5 0.5];
% get microsaccades
for c_Ori2 = 1:length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
            [MS_static{c_Ori2,c_Ori1},numMS_static{c_Ori2,c_Ori1}] = findMicroSaccades(eyeSpeedData.StaticMag{c_Ori2,c_Ori1},cutOff,timeValsEyeData,timeRange); %#ok<*AGROW>
            [MS_CP{c_Ori2,c_Ori1},numMS_CP{c_Ori2,c_Ori1}] = findMicroSaccades(eyeSpeedData.CPMag{c_Ori2,c_Ori1},cutOff,timeValsEyeData,timeRange);
    end
end

microsaccadeData.static = MS_static;
microsaccadeData.CP = MS_CP;
microsaccadeData.static_Count = numMS_static;
microsaccadeData.CP_Count = numMS_CP;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Accessory Functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get parameter combinations
function [parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,...
    cValsUnique,tValsUnique,aValsUnique2,eValsUnique2,sValsUnique2,...
    fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract) %#ok<*STOUT>

load(fullfile(folderExtract,'parameterCombinations.mat')); %#ok<*LOAD>

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;            end

end

% Get Bad Trials
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end

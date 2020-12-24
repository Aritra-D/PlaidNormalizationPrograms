function [eyeDataDegX,eyeDataDegY,eyeRangeMS,trialNums,protNumForTrial,FsEyes] = getEyeDataIndividualMonkey(monkeyName,cleanDataFolder)

FsEye = 200;
OrientationTuningFlag = 0;
gridType = 'microelectrode';
[expDates, protocolNames,~,~,~] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag);

eyeDataDegX = []; eyeDataDegY = []; trialNums = []; protNumForTrial = [];

for iProt = 1:length(expDates)
    disp(['Processing Eye Data: Monkey: ' monkeyName, ', expDate: ' expDates{iProt}, ', protocolName: ' protocolNames{iProt}])
    folderExtract = fullfile(cleanDataFolder,'data',monkeyName,gridType,expDates{iProt},protocolNames{iProt},'extractedData');
    folderSegment = fullfile(cleanDataFolder,'data',monkeyName,gridType,expDates{iProt},protocolNames{iProt},'segmentedData');
    
    clear eyeDataExtracted eyeDataSegmented
    
    % Loading Extracted Eye data % contains eyePosDataX, eyePosDataY,
    % eyeCal values for individual trials and eyeRangeMS 
    eyeDataExtracted = load(fullfile(folderExtract,'EyeData'),'eyeData','eyeRangeMS');
    eyeRangeMS = eyeDataExtracted.eyeRangeMS;
%     timeValsEyeData = (timeRangeEyeData(1):1000/FsEye:timeRangeEyeData(2)-1000/FsEye)/1000;
    
    % Loding Segmented Eye data % contains eyeDataDegX and eyeDa
    eyeDataSegmented = load(fullfile(folderSegment,'eyeData\eyeDataDeg'),'eyeDataDegX','eyeDataDegY');
    eyeData.eyeDataDegX = eyeDataSegmented.eyeDataDegX;
    eyeData.eyeDataDegY = eyeDataSegmented.eyeDataDegY;
    
    % combining eyeData across days/sessions
    eyeDataDegX = cat(1,eyeDataDegX,eyeData.eyeDataDegX);
    eyeDataDegY = cat(1,eyeDataDegY,eyeData.eyeDataDegY);
    trialNums = cat(1,trialNums,(1:size(eyeData.eyeDataDegX,1))');
    protNumForTrial = cat(1,protNumForTrial,repmat(iProt,size(eyeData.eyeDataDegX,1),1));
    FsEyes{iProt} = FsEye; %#ok<AGROW> % single(1/(timeValsEyeData(2)-timeValsEyeData(1))); %#ok<AGROW> % FsEye is constant across protocols.

end
end
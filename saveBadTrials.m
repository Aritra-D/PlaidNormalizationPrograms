function saveBadTrials(monkeyName,gridType,protocolType)

% electrode Parameters set for Plaid Protocols
elecParams.spikeCutoff = 15;
elecParams.snrCutoff = 2;
elecParams.dRange = [0 0.75];
elecParams.unitID = 0;
elecParams.oriSelectiveFlag = 0;
elecParams.getSpikeElectrodesFlag = 1;

timeRangeFRComputation = [0.15 0.4];
versionNum = 2;
contrastIndexList = {[1,1],[5,5]};
OrientationTuningFlag = 0; %  OrientationTuningFlag, 0= Plaid Protocols, 1= OriTuningProtocols prior to plaid Protocol
% find data Information for Ori tuning protocols
[expDates_Plaid, protocolNames_Plaid,~,~,dataFolderSourceString] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag);
% find data Information for Plaid protocols
OrientationTuningFlag = 1;
[expDates_oriTuning, protocolNames_oriTuning,~,~,~] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag);

% Fixed params for finding Bad Trials
processAllElectrodes = 0;
threshold = 6;
minLimit = -2000;
maxLimit = 1000;
saveDataFlag = 1;
checkPeriod = [-0.5 0.5];
rejectTolerance = 1;
showElectrodes = [];


for i = 1:length(expDates_Plaid)
    
    [~,~,~,~,~,~,goodElectrodes] ...
        = getGoodElectrodesSingleSession(monkeyName,expDates_Plaid{i},protocolNames_Plaid{i},...
        gridType,elecParams,timeRangeFRComputation,dataFolderSourceString,versionNum,...
        contrastIndexList);
    
    checkTheseElectrodes = goodElectrodes;
    
    if strcmp(protocolType,'Plaid')
        % find Bad Trials for Plaid Protocols
        findBadTrialsWithLFPv3(monkeyName,expDates_Plaid{i},protocolNames_Plaid{i},dataFolderSourceString,gridType, ...
            checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,saveDataFlag,checkPeriod,rejectTolerance);
            
        
        
    elseif strcmp(protocolType,'oriTuning')
        % find Bad Trials for oriTuning Protocols
        findBadTrialsWithLFPv3(monkeyName,expDates_oriTuning{i},protocolNames_oriTuning{i},dataFolderSourceString,gridType, ...
            checkTheseElectrodes,processAllElectrodes,threshold,maxLimit,minLimit,showElectrodes,saveDataFlag,checkPeriod,rejectTolerance);
           
        
        
    end
end
end
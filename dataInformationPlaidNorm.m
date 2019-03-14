% dataInformation for Plaid Normalization

function [expDates, protocolNames,positionList,oriList,dataFolderSourceString] = dataInformationPlaidNorm(monkeyName,gridType,OrientationTuningFlag)
if ~exist('gridType','var');                            gridType = 'microelectrode'; end
if ~exist('OrientationTuningFlag','var');               OrientationTuningFlag = 0;   end 
% 0 = Plaid Protocols, 
% 1 = OrientationTuningProtocol 
%(grating placed on the grid centre with orientation varying from 0-157.5
% in 8 steps of 22.5 degrees)

if strcmp(getenv('computername'),'RAYLABPC-ARITRA')
    dataFolderSourceString = 'E:\';
else
    dataFolderSourceString = 'M:\Data\PlaidNorm';
end

if strcmp(monkeyName,'alpaH')
    [allExpDates,allProtocolNames,~,~] = getAllProtocols(monkeyName, gridType);
    if OrientationTuningFlag == 1
       protocolList = [306 308 312 314 319 321 324 330 355 359 361 372 374];
        
    elseif OrientationTuningFlag == 0
       protocolList = [307 309 313 315 320 322 325 331 356 360 362 373 375]; 
            
    end
    
    numDays = length(protocolList);
    expDates = cell(1,numDays);
    protocolNames = cell(1,numDays);
    
    for i=1:numDays
    expDates{i} = allExpDates{protocolList(i)};
    protocolNames{i} = allProtocolNames{protocolList(i)};
    end
    
    positionList = [3.3 -2.8];
    
    oriList{1} = [0 90];
    oriList{2} = [45 135];
    oriList{3} = [22.5 112.5];
    oriList{4} = [67.5 157.5]; % add later after combining 315 and 316 into 315
    oriList{5} = [0 90];
    oriList{6} = [45 135];
    oriList{7} = [22.5 112.5];
    oriList{8} = [67.5 157.5];
    oriList{9} = [22.5 112.5];
    oriList{10} = [0 90];
    oriList{11} = [45 135];
    oriList{12} = [22.5 112.5];
    oriList{13} = [67.5 157.5];
    
elseif strcmp(monkeyName,'kesariH')
    
    [allExpDates,allProtocolNames,~,~] = getAllProtocols(monkeyName, gridType);
    if OrientationTuningFlag == 1
       protocolList = [70 94 111 114 116 118 120 122 124]; % 106 not used; OriTuning Protocol for 107
        
    elseif OrientationTuningFlag == 0
       protocolList = [71 95 112 115 117 119 121 123 125]; % 107 not used; driifting instead of counter-phase
            
    end
    
    numDays = length(protocolList);
    expDates = cell(1,numDays);
    protocolNames = cell(1,numDays);
    
    for i=1:numDays
    expDates{i} = allExpDates{protocolList(i)};
    protocolNames{i} = allProtocolNames{protocolList(i)};
    end
    
    positionList = [0.5 -1.75];
    
    oriList{1} = [0 90];
    oriList{2} = [45 135];
%     oriList{3} = [22.5 112.5]; % 107; Drifting grating, no SSVEP measure;
%     excluded
    oriList{3} = [22.5 112.5];
    oriList{4} = [67.5 157.5];
    oriList{5} = [0 90];
    oriList{6} = [45 135];
    oriList{7} = [22.5 112.5];
    oriList{8} = [67.5 157.5];
    oriList{9} = [22.5 112.5];
  
    
end
end
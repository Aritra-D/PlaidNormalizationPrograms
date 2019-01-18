% dataInformation for Plaid Normalization

function [expDates, protocolNames,positionList,oriList] = dataInformationPlaidNorm(monkeyName,gridType, OrientationTuningFlag)
if ~exist('gridType','var');                            gridType = 'microelectrode'; end
if ~exist('OrientationTuningFlag','var');               OrientationTuningFlag = 0;      end 
% 0 = Plaid Protocols, 
% 1 = OrientationTuningProtocol 
%(grating placed on the grid centre with orientation varying from 0-157.5
% in 8 steps of 22.5 degrees)

if strcmp(monkeyName,'alpaH')
    [allExpDates,allProtocolNames,~,~] = getAllProtocols(monkeyName, gridType);
    if OrientationTuningFlag == 1
       protocolList = [306 308 312 319 321 324 330 355 359 361 372 374];
        
    elseif OrientationTuningFlag == 0
       protocolList = [307 309 313 320 322 325 331 356 360 362 373 375]; 
            
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
    oriList{4} = [0 90];
    oriList{5} = [45 135];
    oriList{6} = [22.5 112.5];
    oriList{7} = [67.5 157.5];
    oriList{8} = [22.5 112.5];
    oriList{9} = [0 90];
    oriList{10} = [45 135];
    oriList{11} = [22.5 112.5];
    oriList{12} = [67.5 157.5];
    
elseif strcmp(monkeyName,'kesariH')
    
    [allExpDates,allProtocolNames,~,~] = getAllProtocols(monkeyName, gridType);
    if OrientationTuningFlag == 1
       protocolList = [70 94 106 111 114 116 118 120 122 124];
        
    elseif OrientationTuningFlag == 0
       protocolList = [71 95 107 112 115 117 119 121 123 125]; 
            
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
    oriList{3} = [22.5 112.5]; % Drifting grating, no SSVEP measure
    oriList{4} = [22.5 112.5];
    oriList{5} = [67.5 157.5];
    oriList{6} = [0 90];
    oriList{7} = [45 135];
    oriList{8} = [22.5 112.5];
    oriList{9} = [67.5 157.5];
    oriList{10} = [22.5 112.5];
  
    
end
end
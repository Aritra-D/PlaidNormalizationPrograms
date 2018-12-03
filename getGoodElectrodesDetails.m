function [tmpElectrodeStringList,tmpElectrodeArrayList,allElecs,monkeyNameList] = getGoodElectrodesDetails
% gridType = 'microelectrode';
monkeyNameList{1} = 'alpaH'; monkeyNameList{2} = 'kesariH';
gridType = 'microelectrode';
tmpElectrodeStringList = cell(1,2);
tmpElectrodeArrayList = cell(1,2);
allElecs = zeros(1,2);
for i=1:2
    disp(['MonkeyName: ' monkeyNameList{i}])
    clear expDates protocolNames
    [~,protocolNames,~,~]= dataInformationPlaidNorm(monkeyNameList{i},gridType,0);
    numSessions = size(protocolNames,2);
    tmpElectrodeListArray = cell(1,numSessions);
    tmpElectrodeListStr = cell(1,numSessions);
    numElecs = 0;
    for SessionNum = 1:numSessions
        [tmpElectrodeListArray{SessionNum},tmpElectrodeListStr{SessionNum},goodElectrodes] = getGoodElectrodesSingleSession(monkeyNameList{i},gridType,SessionNum);
        numElecs = numElecs+length(goodElectrodes);
    end
    tmpElectrodeStringList{i} = tmpElectrodeListStr;
    tmpElectrodeArrayList{i} = tmpElectrodeListArray;
    allElecs(i) = numElecs;
end
end
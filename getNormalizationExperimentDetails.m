function [fileNameStringList,monkeyNameList] = getNormalizationExperimentDetails

gridType = 'microelectrode';
monkeyNameList{1} = 'alpaH'; monkeyNameList{2} = 'kesariH';

fileNameStringList = cell(1,2);
for i=1:2
    clear expDates protocolNames
    [expDates,protocolNames,~,~]= dataInformationPlaidNorm(monkeyNameList{i},gridType,0);
    numSessions = size(protocolNames,2);
    tmpFileNameList = cell(1,numSessions);
    for j = 1:numSessions
    tmpFileNameList{j} = [monkeyNameList{i} expDates{j} protocolNames{j}];
    end
    fileNameStringList{i} = tmpFileNameList;
end

end
% Saves Impedance Values for each of the protocols
monkeyName = 'kesariH';
gridType = 'Microelectrode';


[expDates, protocolNames,positionList,folderSourceString] = dataInformationPlaidNorm(monkeyName,gridType);

numDays = length(expDates);

for i=1:numDays
    expDate = expDates{i};
    getImpedanceDataPlaidNorm(monkeyName,expDate,folderSourceString,gridType)
    
end
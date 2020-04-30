% Saves Impedance Values for each of the protocols
monkeyName = 'alpaH';
gridType = 'Microelectrode';


[expDates, protocolNames,oriList,folderSourceString] = dataInformationPlaidProtocolsEEG(monkeyName);

numDays = length(expDates);

for i=1:numDays
    expDate = expDates{i};
    getImpedanceDataPlaidNorm(monkeyName,expDate,folderSourceString,gridType)
    
end
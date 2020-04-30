function getImpedanceDataPlaidNorm(subjectName,expDate,folderSourceString,gridType)

folderSourceStringRawData = 'M:';
folderNameRawData = fullfile(folderSourceStringRawData,'Data','rawData',subjectName,[subjectName expDate]);
folderNameSave = fullfile(folderSourceString,'data',subjectName,gridType,expDate);
makeDirectory(folderNameSave);

fileName = fullfile(folderNameRawData,[subjectName expDate 'impedance']);
if ~exist(fileName,'file')    
    fileName = [fileName '.txt'];             
end

if ~exist(fileName,'file')
    disp([fileName ' does not exist']);
else
    X = textread(fileName,'%s'); %#ok<DTXTRD>
    
    impedanceValues = zeros(1,107);
    for i=1:107
        for j=1:length(X)
            if strcmp(X{j},['elec1-' num2str(i)])||strcmp(X{j},['elec2-' num2str(i)])||strcmp(X{j},['chan' num2str(i)]) %saves in a sorted manner
                impedanceValues(i) = str2double(X{j+1});
                break;
            end
        end
    end
    % Save
    save(fullfile(folderNameSave,'impedanceValues.mat'),'impedanceValues');
end
end
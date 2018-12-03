
% Gets GoodPos for all Stim Combinations (5 cons and 2 TFs)

function [goodPosAll, sizeGoodPos] = getGoodPosPlaid(MonkeyName,expDate,protocolName,folderSourceString,gridType)

if ~exist('folderSourceString','var');  folderSourceString='E:';        end         % Default values to be taken if function arguments (variables) not defined! 
if ~exist('gridType','var');            gridType='microelectrode';      end

folderName = fullfile(folderSourceString,'data',MonkeyName,gridType,expDate,protocolName);
folderExtract = fullfile(folderName,'extractedData');
load(fullfile(folderExtract,'parameterCombinations.mat'));

 
for a = 1: length(aValsUnique)
    for e = 1: length(eValsUnique)
        for s = 1: length(sValsUnique)
            for f = 1: length(fValsUnique)
                for o = 1: length(oValsUnique)
                    for c = 1: length(cValsUnique)
                        for t = 1: length(tValsUnique)
                            goodPos = parameterCombinations{a,e,s,f,o,c,t};
                            goodPos2 = parameterCombinations2{a,e,s,f,o,c,t};
                            LeftGoodPos{c,t} = goodPos; 
                            RightGoodPos{c,t} = goodPos2;
                        end
                    end
                end
            end
        end
    end
end

LeftStimIndex = reshape(LeftGoodPos,[10,1]);
RightStimIndex = reshape(RightGoodPos,[10,1]);
for iCombLeft = 1: length(LeftStimIndex)
    for iCombRight = 1: length(RightStimIndex)
        goodPosAll{iCombLeft,iCombRight} = intersect(LeftStimIndex{iCombLeft,1},RightStimIndex{iCombRight,1}) ;
        sizeGoodPos(iCombLeft,iCombRight) = size(goodPosAll{iCombLeft,iCombRight},2);
    end
end

disp(['max Number of Stim Reps = ', num2str(max(max(sizeGoodPos,[],1),[],2))]);
disp(['min Number of Stim Reps = ', num2str(min(min(sizeGoodPos,[],1),[],2))]);
AllStimReps = reshape(sizeGoodPos,[100,1]);
disp(['avg Number of Stim Reps = ', num2str(mean(AllStimReps))]);
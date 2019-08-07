function [diffOut,predictedResponseMatrix] = getResponseMatrixGrating(params,actualResponseMatrix)

cList = [0 1 2 4 8]/16;

L=params(1);
sigma = min(max(0,params(2)),5);

predictedResponseMatrix = (cList * L) ./ (cList+sigma);

if ~isempty(actualResponseMatrix)         
    diffOut = sum((actualResponseMatrix(:) - predictedResponseMatrix(:)).^2);
else
    diffOut = [];
end
end
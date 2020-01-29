function params = getParametersPlaid(actualResponseMatrix,versionNum)

opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
    'Display','off','LargeScale','off','MaxFunEvals',500);

if versionNum==1
sigmaStart = 0.1;
LStart = (1+2*sigmaStart)*actualResponseMatrix(1,1);
alphaStart = 2*LStart/actualResponseMatrix(1,5) - (1+2*sigmaStart);

startPt = [LStart alphaStart sigmaStart];

elseif versionNum==2
sigmaStart = 0.1;
L1Start = (1+2*sigmaStart)*actualResponseMatrix(5,5);
L2Start = (1+2*sigmaStart)*actualResponseMatrix(1,1);
alphaStart = (L1Start+L2Start)/actualResponseMatrix(1,5) - (1+2*sigmaStart);

startPt = [L1Start L2Start alphaStart sigmaStart];
end

params = fminsearch(@(params) getResponseMatrixPlaid(params,actualResponseMatrix,versionNum),startPt,opts);
end
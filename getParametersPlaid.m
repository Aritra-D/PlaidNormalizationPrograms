function params = getParametersPlaid(actualResponseMatrix)

opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
    'Display','off','LargeScale','off','MaxFunEvals',500);

sigmaStart = 0.1;
LStart = (1+2*sigmaStart)*actualResponseMatrix(1,1);
alphaStart = 2*LStart/actualResponseMatrix(1,5) - (1+2*sigmaStart);

startPt = [LStart alphaStart sigmaStart];

params = fminsearch(@(params) getResponseMatrixPlaid(params,actualResponseMatrix),startPt,opts);
end
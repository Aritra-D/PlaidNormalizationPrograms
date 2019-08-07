function params = getParametersGrating(actualResponseMatrix)

opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
    'Display','off','LargeScale','off','MaxFunEvals',500);

LStart = actualResponseMatrix(end);
sigmaStart = 0.1;

startPt = [LStart sigmaStart];
params = fminsearch(@(params) getResponseMatrixGrating(params,actualResponseMatrix),startPt,opts);
end
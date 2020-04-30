function params = getParametersPlaidV2(actualResponseMatrix,modelNum)

opts = optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
    'Display','off','LargeScale','off','MaxFunEvals',500);

if modelNum==1 % Free parameters L1, L2, alpha, and sigma; Divisive normalization model
    sigmaStart = 0.1;
    L1Start = (1+2*sigmaStart)*actualResponseMatrix(5,5);
    L2Start = (1+2*sigmaStart)*actualResponseMatrix(1,1);
    
    startPt = [L1Start L2Start sigmaStart];
    
elseif modelNum==2 % Free parameters L1, L2, alpha, and sigma; Stimulus tuned normalization model
    sigmaStart = 0.1;
    L1Start = (1+2*sigmaStart)*actualResponseMatrix(5,5);
    L2Start = (1+2*sigmaStart)*actualResponseMatrix(1,1);
    alphaStart = (L1Start+L2Start)/actualResponseMatrix(1,5) - (1+2*sigmaStart);
    
    startPt = [L1Start L2Start alphaStart sigmaStart];
    
elseif modelNum==3 % Free parameters L1, L2, alpha, and sigma; EMS-Stimulus tuned normalization model
    sigmaStart = 0.1;
    L1Start = (1+2*sigmaStart)*actualResponseMatrix(5,5);
    L2Start = (1+2*sigmaStart)*actualResponseMatrix(1,1);
    alphaStart = (L1Start+L2Start)/actualResponseMatrix(1,5) - (1+2*sigmaStart);
    
    startPt = [L1Start L2Start alphaStart sigmaStart];
    
elseif modelNum==4 % Free parameters L, alpha, and sigma; EMS-Stimulus tuned normalization model for symmetrized data
    sigmaStart = 0.1;
    LStart = (1+2*sigmaStart)*actualResponseMatrix(1,1);
    alphaStart = 2*LStart/actualResponseMatrix(1,5) - (1+2*sigmaStart);
    
    startPt = [LStart alphaStart sigmaStart];
    
end

params = fminsearch(@(params) getResponseMatrixPlaidV2(params,actualResponseMatrix,modelNum),startPt,opts);
end
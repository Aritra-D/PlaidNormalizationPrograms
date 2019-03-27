function params = getParameters(actualResponseMatrix)

opts = optimoptions('fmincon');
opts.MaxIterations = 1e4;
% opts.Display = 'off';
opts.MaxFunctionEvaluations = 1e4;
% optimset('TolX',1e-6,'TolFun',1e-6,'MaxIter',5000,...
%     'Display','off','LargeScale','off','MaxFunEvals',500);

startPt = rand(1,5);
% startPt = [(actualResponseMatrix(5,5)+actualResponseMatrix(1,1))/2 (actualResponseMatrix(5,5)+actualResponseMatrix(1,1))/2 1 1 0.2];
params = fmincon(@(params) getResponseMatrix(params,actualResponseMatrix),startPt,[],[],[],[],zeros(1,5),[inf(1,4) 1],[],opts);
end
function [diffOut,predictedResponseMatrix] = getResponseMatrixPlaid(params,actualResponseMatrix)

c1List = [0 1 2 4 8]/16;
c2List = c1List;

L=params(1);
alpha = max(0,params(2));
sigma=min(max(0,params(3)),5);

for i=1:length(c1List)
    c1=c1List(i);
    for j=1:length(c2List)
        c2 = c2List(j);
        predictedResponseMatrix(i,j) = (c1*L)/(c1+alpha*c2+sigma) + (c2*L)/(c2+alpha*c1+sigma); %#ok<AGROW>
    end
end

predictedResponseMatrix = flip(predictedResponseMatrix,1);

if ~isempty(actualResponseMatrix)         
    diffOut = sum((actualResponseMatrix(:) - predictedResponseMatrix(:)).^2);
else
    diffOut = [];
end
end
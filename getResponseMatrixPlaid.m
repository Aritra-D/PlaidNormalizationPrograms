function [diffOut,predictedResponseMatrix] = getResponseMatrixPlaid(params,actualResponseMatrix,versionNum)

c1List = [0 1 2 4 8]/16;
c2List = c1List;

if versionNum==1
    
    L=params(1);
    alpha = params(2);
    sigma=min(max(0,params(3)),5);
    
    for i=1:length(c2List)
        c2=c2List(i);
        for j=1:length(c1List)
            c1 = c1List(j);
            predictedResponseMatrix(i,j) = (c1*L)/(c1+alpha*c2+sigma) + (c2*L)/(c2+alpha*c1+sigma); %#ok<AGROW>
        end
    end
end

if versionNum==2
    
    L1=params(1);
    L2=params(2);
    alpha = params(3);
    sigma=min(max(0,params(4)),5);
    
    for i=1:length(c1List)
        c2=c2List(i);
        for j=1:length(c2List)
            c1 = c1List(j);
            predictedResponseMatrix(i,j) = (c2*L2)/(c2+alpha*c1+sigma)+(c1*L1)/(c1+alpha*c2+sigma); 
        end
    end
end


predictedResponseMatrix = flip(predictedResponseMatrix,1);

if ~isempty(actualResponseMatrix)
    diffOut = sum((actualResponseMatrix(:) - predictedResponseMatrix(:)).^2);
else
    diffOut = [];
end
end
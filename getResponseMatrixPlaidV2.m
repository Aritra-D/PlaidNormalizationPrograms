% DaactualResponseMatrix is organized into Response for Grating 2 (Rows) x Response for
% Grating 1 (Columns); positive x-axis in rows indicate increase in
% contrast of Grating 1, positivey-axis along columns inicate increase in
% contrast of Grating 2

function [diffOut,predictedResponseMatrix] = getResponseMatrixPlaidV2(params,actualResponseMatrix,modelNum)

c1List = [0 1 2 4 8]/16;
c2List = c1List;

if modelNum==1
    L1=params(1);
    L2=params(2);
    sigma=min(max(0,params(3)),5);
    
    for i=1:length(c2List)
        c2=c2List(i);
        for j=1:length(c1List)
            c1 = c1List(j);
            predictedResponseMatrix(i,j) = (c2*L2+c1*L1)/(c1+c2+sigma);
        end
    end
   
elseif modelNum==2
    L1=params(1);
    L2=params(2);
    alpha = params(3);
    sigma=min(max(0,params(4)),5);
    
    for i=1:length(c2List)
        c2=c2List(i);
        for j=1:length(c1List)
            c1 = c1List(j);
            predictedResponseMatrix(i,j) = (c2*L2+c1*L1)/(c1+alpha*c2+sigma); %#ok<*AGROW>
        end
    end
    
elseif modelNum==3
    L1=params(1);
    L2=params(2);
    alpha = params(3);
    sigma=min(max(0,params(4)),5);
    
    for i=1:length(c2List)
        c2=c2List(i);
        for j=1:length(c1List)
            c1 = c1List(j);
            predictedResponseMatrix(i,j) = (c2*L2)/(c2+alpha*c1+sigma)+(c1*L1)/(c1+alpha*c2+sigma);
        end
    end
    
elseif modelNum==4
    L=params(1);
    alpha = params(2);
    sigma=min(max(0,params(3)),5);
    
    for i=1:length(c2List)
        c2=c2List(i);
        for j=1:length(c1List)
            c1 = c1List(j);
            predictedResponseMatrix(i,j) = (c1*L)/(c1+alpha*c2+sigma) + (c2*L)/(c2+alpha*c1+sigma);
        end
    end
    
end


predictedResponseMatrix = flip(predictedResponseMatrix,1); % data is flipped to match organization of actual response matrix (see notes on top)

if ~isempty(actualResponseMatrix)
    diffOut = sum((actualResponseMatrix(:) - predictedResponseMatrix(:)).^2);
else
    diffOut = [];
end
end
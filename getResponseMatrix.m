function [diffOut,predictedResponseMatrix] = getResponseMatrix(params,actualResponseMatrix)

c1List = [0 1 2 4 8]/16;
c2List = c1List;
% c1List = [0:0.1:8]/16;
% c2List = 0;


L1=params(1);
L2=params(2);
alpha1=max(0,params(3));
alpha2=max(0,params(4));
sigma=max(0,params(5));

for i=1:length(c1List)
    c1=c1List(i);
    for j=1:length(c2List)
        c2 = c2List(j);
        predictedResponseMatrix(i,j) = (c1*L1+c2*L2)/(alpha1*(c1)+alpha2*(c2)+sigma); %#ok<AGROW>
%         predictedResponseMatrix(i,j) = ((c1)^2*L1+(c2)^2*L2)/(alpha1*(c1)^2+alpha2*(c2)^2+(sigma)^2); %#ok<AGROW>
    end
end

predictedResponseMatrix = flip(predictedResponseMatrix,1);

% subplot(2,2,1)
% plot(c1List,actualResponseMatrix(end,:),'Marker','o','LineStyle','none','color','r')
% hold on;
% plot(c1List,predictedResponseMatrix(end,:),'LineWidth',2,'color','r')
% plot(c1List,diag(flipud(actualResponseMatrix)),'Marker','o','LineStyle','none','color','b')
% plot(c1List,diag(flipud(predictedResponseMatrix)),'LineWidth',2,'color','b')



if ~isempty(actualResponseMatrix)
    diffOut = sum(sum((actualResponseMatrix - predictedResponseMatrix).^2));
else
    diffOut = [];
end
end
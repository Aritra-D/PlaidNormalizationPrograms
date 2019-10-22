function plotBars(hPlot,fileSave1,fileSave3,colorScheme,Criterion)

if exist(fileSave1,'file')
    load(fileSave1)
    disp(['Loading file ' fileSave1]);
end
spikeRateData = firingRateData;
clear erpData firingRateData fftData energyData energyDataTF oriTuningData NI_Data

if exist(fileSave3,'file')
    load(fileSave3)
    disp(['Loading file ' fileSave3]);
end
clear erpData firingRateData fftData energyDataTF oriTuningData NI_Data

xAll{1} = squeeze(spikeRateData.analysisDataST(:,1,:,:)) - squeeze(spikeRateData.analysisData_cBL(:,1,:,:));
type{1} = 'FR';

for i=2:4
    xAll{i} = squeeze(energyData.analysisDataST{i}) - squeeze(energyData.analysisData_cBL{i});
end

type{2} = 'G'; type{3} = 'HG'; type{4} = 'S';


for i=1:4
    x = xAll{i};
    numElectrodes = length(xAll{i});

    for j=1:numElectrodes
        
        % Grating
        g1 = squeeze(x(j,5,:))';
        g2 = flip(squeeze(x(j,:,1)),2);
        
        g = (g1+g2)/2; % Make symmetric
        parG = getParametersGrating(g);
        [dg,pg] = getResponseMatrixGrating(parG,g);
        eG(j) = 1 - (dg/sum((g-mean(g)).^2));
        sG(j) = parG(2);
        dataG(j,:) = g;
        
        %Plaid
        y = squeeze(x(j,:,:));
        z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
        parP = getParametersPlaid(z);
        [dp,pz] = getResponseMatrixPlaid(parP,z);
        eP(j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
        aP(j) = parP(2);
        sP(j) = parP(3);
        dataP(j,:,:) = z;
    end
    
    eP_All{i} = eP;
    aP_All{i} = aP;
    sP_All{i} = sP;
    dataP_All{i} = dataP;
end

for i=1:4
dataP =  dataP_All{i};   
NI{i} = squeeze(2*dataP(:,1,1)./dataP(:,1,5))';
end

% median_NI = cellfun(@median,NI);
% median_aP
if strcmp(colorScheme,'color')
    colors = [0 0 1; 0 1 0; 1 0 0; 0 0.75 0.75];
elseif strcmp(colorScheme,'grayscale')
    colors = repmat(0.55:-0.1:0.25,[3 1])';
end

% NI_data = cellfun(@median,NI);
% aP_data = cellfun(@median,aP_All);
% sP_data = cellfun(@median,sP_All);
% eP_data = cellfun(@median,eP_All);

% for i=1:4
%     hNI = bar(hPlot.hPlot1(3),1:4,NI_data);
%     hold(hPlot.hPlot1(2),'on')
%     errorbar(hPlot.hPlot1(2),i,median(NI{i}),getSEMedian(NI{i},1000),'color','k','marker','.','lineStyle','-')
% end
%     bar(hPlot.hPlot1(2),i,aP_data(i),'faceColor',colors(i,:,:));
%     hold(hPlot.hPlot1(2),'on')
%     bar(hPlot.hPlot1(3),i,sP_data(i),'faceColor',colors(i,:,:));
%     hold(hPlot.hPlot1(3),'on')
%     bar(hPlot.hPlot1(4),i,eP_data(i),'faceColor',colors(i,:,:));
%     hold(hPlot.hPlot1(4),'on')

% hNI = bar(hPlot.hPlot1(1),1:4,cellfun(@median,NI));
% hNI(1).FaceColor = 'r';
% 
% set(hNI(1),'faceColor','k');
elecNums = size( NI{1},2);
if strcmp(Criterion,'aP')
   elecIDs = (aP_All{1}>0);
elseif strcmp(Criterion,'eP')
   elecIDs = (eP_All{1}>0.5);
elseif strcmp(Criterion,'aP&eP')
   elecIDs = (aP_All{1}>0 & eP_All{1}>0.5);
else
    elecIDs = 1:elecNums;
end

NI{1} = NI{1}(elecIDs);
aP_All{1} = aP_All{1}(elecIDs);
sP_All{1} = sP_All{1}(elecIDs);
eP_All{1} = eP_All{1}(elecIDs);
fraction = size(NI{1},2)/elecNums;





bar(hPlot.hPlot1(1),cellfun(@median,NI),'facecolor',[0.5 0.5 0.5]);
bar(hPlot.hPlot1(2),cellfun(@median,aP_All),'facecolor',[0.5 0.5 0.5]);
bar(hPlot.hPlot1(3),cellfun(@median,sP_All),'facecolor',[0.5 0.5 0.5]);
bar(hPlot.hPlot1(4),cellfun(@median,eP_All),'facecolor',[0.5 0.5 0.5]);
for i=1:4
    hold(hPlot.hPlot1(1),'on');
        
        errorbar(hPlot.hPlot1(1),i,median(NI{i}),getSEMedian(NI{i},1000),'color','k','marker','.','lineStyle','-')

        hold(hPlot.hPlot1(2),'on');
    errorbar(hPlot.hPlot1(2),i,median(aP_All{i}),getSEMedian(aP_All{i},1000),'color','k','marker','.','lineStyle','-')
    hold(hPlot.hPlot1(3),'on');
    errorbar(hPlot.hPlot1(3),i,median(sP_All{i}),getSEMedian(sP_All{i},1000),'color','k','marker','.','lineStyle','-')
    hold(hPlot.hPlot1(4),'on');
    errorbar(hPlot.hPlot1(4),i,median(eP_All{i}),getSEMedian(eP_All{i},1000),'color','k','marker','.','lineStyle','-')
end

% for i=1:4
    for j = 1: size(NI{1},2)
        plot(hPlot.hPlot1(1),1,NI{1}(j),'color','k','marker','o')
        plot(hPlot.hPlot1(2),1,aP_All{1}(j),'color','k','marker','o')
        plot(hPlot.hPlot1(3),1,sP_All{1}(j),'color','k','marker','o')
        plot(hPlot.hPlot1(4),1,eP_All{1}(j),'color','k','marker','o')

    end
% end

barNames = {'spikeRate','gamma','hi-gamma','SSVEP'};

tickLengthPlot = 2*get(hPlot.hPlot1(1),'TickLength');
for i=1:4
    ylabel(hPlot.hPlot1(i),'Population Median')
    set(hPlot.hPlot1(i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlot.hPlot1(i),'xticklabel',barNames)
end

title(hPlot.hPlot1(1),'N.I.','fontSize',20);
title(hPlot.hPlot1(2),'\alpha','fontSize',20);
title(hPlot.hPlot1(3),'\sigma','fontSize',20);
title(hPlot.hPlot1(4),['% Explained Variance,elecsFrac ' num2str(round(fraction,2))] ,'fontSize',20);






end
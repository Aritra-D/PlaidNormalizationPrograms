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

combineUniqueElectrodes=1;
N=15; snr=2; d=0.75;

%%%%%%%%%%%%%%%% Subselect electrodes based on N, snr and d %%%%%%%%%%%%%%%
goodN = (max([elecInfo.N(1,:);elecInfo.N(2,:)])>N);
goodSNR = (elecInfo.SNR>snr);
goodD = (elecInfo.d <= d);
goodElectrodeIDs = (goodN & goodSNR & goodD);

%%%%%%%%%%%%%%%%%%% Get information about the electrodes %%%%%%%%%%%%%%%%%%
allElectrodes = [];
for i=1:length(elecInfo.elecs)
    x = elecInfo.elecs{i}{end};
    if i<=13
        allElectrodes = cat(2,allElectrodes,x);
    else
        allElectrodes = cat(2,allElectrodes,x+100); % We add 100 to each electrode ID in Monkey 2 so that electrode numbers are different from monkey 1
    end
end

%%%%%% For each electrode, find the IDs that need to be averaged %%%%%%%%%%
clear goodIDList
if combineUniqueElectrodes
    
    uniqueGoodElectrodeList = unique(allElectrodes(goodElectrodeIDs));
    numGoodElectrodes = length(uniqueGoodElectrodeList);
    
    goodIDList = cell(1,numGoodElectrodes);
    for i=1:numGoodElectrodes
        goodIDList{i} = find(allElectrodes==uniqueGoodElectrodeList(i));
    end
else
    goodIDs = find(goodElectrodeIDs==1); %#ok<*UNRCH>
    numGoodElectrodes = length(goodIDs);
    
    goodIDList = cell(1,numGoodElectrodes);
    for i=1:numGoodElectrodes
        goodIDList{i} = goodIDs(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%% Combine data if necessary %%%%%%%%%%%%%%%%%%%%%%%%%
% allData = cell(1,4);
clear allDataTMP
allDataTMP = zeros(numGoodElectrodes,5,5);
for i=1:numGoodElectrodes
    clear x
    x = squeeze(spikeRateData.analysisDataST(goodIDList{i},1,:,:) - spikeRateData.analysisData_cBL(goodIDList{i},1,:,:));
    
    if length(goodIDList{i})==1
%         y = x;
%         z = flip(y,1); z = (z+z')/2; z = flip(z,1);
        allDataTMP(i,:,:) = x;
    else
        xs = zeros(size(x));
        for j=1:length(goodIDList{i})
%             y = squeeze(x(j,:,:));
%             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
            xs(j,:,:) = squeeze(x(j,:,:));
        end
        allDataTMP(i,:,:) = squeeze(mean(xs,1));
    end
end
xAll{1} = allDataTMP;
% xAll{1} = squeeze(spikeRateData.analysisDataST(:,1,:,:)) - squeeze(spikeRateData.analysisData_cBL(:,1,:,:));
type{1} = 'FR';

% allDataType{2} = 'G'; allDataType{3} = 'HG'; allDataType{4} = 'S';
for i=2:4
    eDataTMP = squeeze(energyData.analysisDataST{i}) - squeeze(energyData.analysisData_cBL{i});
    
    allDataTMP = zeros(numGoodElectrodes,5,5);
    for j=1:numGoodElectrodes
        clear x
        x = eDataTMP(goodIDList{j},:,:);
        
        xs = zeros(size(x));
        for k=1:length(goodIDList{j})
%             y = squeeze(x(k,:,:));
%             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
            xs(k,:,:) = squeeze(x(k,:,:));
        end

        allDataTMP(j,:,:) = squeeze(mean(xs,1));
    end
    xAll{i} = allDataTMP; 
end

type{2} = 'G'; type{3} = 'HG'; type{4} = 'S';
versionNum = 3;

for i=1:4
    x = xAll{i};
    numElectrodes = length(xAll{i});

    for j=1:numElectrodes
        
%         % Grating
%         g1 = squeeze(x(j,5,:))';
%         g2 = flip(squeeze(x(j,:,1)),2);
%         
%         g = (g1+g2)/2; % Make symmetric
%         parG = getParametersGrating(g);
%         [dg,pg] = getResponseMatrixGrating(parG,g);
%         eG(j) = 1 - (dg/sum((g-mean(g)).^2));
%         sG(j) = parG(2);
%         dataG(j,:) = g;
        
        %Plaid
        z = squeeze(x(j,:,:));
%         z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
        parP = getParametersPlaidV2(z,versionNum);
        [dp,pz] = getResponseMatrixPlaidV2(parP,z,versionNum); %#ok<*ASGLU>
        eP(j) = 1 - (dp/sum((z(:)-mean(z(:))).^2)); %#ok<*AGROW>
        aP(j) = parP(3);
        sP(j) = parP(4);
        dataP(j,:,:) = z;
    end
    
    eP_All{i} = eP;
    aP_All{i} = aP;
    sP_All{i} = sP;
    dataP_All{i} = dataP;
end

for i=1:4
dataP =  dataP_All{i};   
NI{i} = squeeze((dataP(:,1,1)+dataP(:,5,5))./dataP(:,1,5))';
end

% median_NI = cellfun(@median,NI);
% median_aP
% if strcmp(colorScheme,'color')
%     colors = [0 0 1; 0 1 0; 1 0 0; 0 0.75 0.75];
% elseif strcmp(colorScheme,'grayscale')
%     colors = repmat(0.55:-0.1:0.25,[3 1])';
% end

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
% fraction = size(NI{1},2)/elecNums;





bar(hPlot.hPlot1(1),cellfun(@median,NI),'facecolor',[0.5 0.5 0.5]);
bar(hPlot.hPlot1(2),cellfun(@median,aP_All),'facecolor',[0.5 0.5 0.5]);
bar(hPlot.hPlot1(3),cellfun(@median,sP_All),'facecolor',[0.5 0.5 0.5]);
bar(hPlot.hPlot1(4),cellfun(@median,eP_All),'facecolor',[0.5 0.5 0.5]);

for i=1:4
    sem_NI{i} = getSEMedian(NI{i});
    sem_aP{i} = getSEMedian(aP_All{i});
    sem_sP{i} = getSEMedian(sP_All{i});
    sem_eP{i} = getSEMedian(eP_All{i});
end

for i=1:4
    hold(hPlot.hPlot1(1),'on');
    errorbar(hPlot.hPlot1(1),i,median(NI{i}),sem_NI{i},'color','k','marker','.','lineStyle','-')
    hold(hPlot.hPlot1(2),'on');
    errorbar(hPlot.hPlot1(2),i,median(aP_All{i}),sem_aP{i},'color','k','marker','.','lineStyle','-')
    hold(hPlot.hPlot1(3),'on');
    errorbar(hPlot.hPlot1(3),i,median(sP_All{i}),sem_sP{i},'color','k','marker','.','lineStyle','-')
    hold(hPlot.hPlot1(4),'on');
    errorbar(hPlot.hPlot1(4),i,median(eP_All{i}),sem_eP{i},'color','k','marker','.','lineStyle','-')
end


paramType{1}= 'NI';
paramType{2}= 'aP';
paramType{3}= 'sP';
paramType{4}= 'eP';

for i = 1:4
    disp([paramType{i},':'])
    disp('median +/- SEMedian')
    for j = 1:4
        if i==1
            disp([type{j},':' num2str(round(median(NI{j}),2)),' +/- ' num2str(round(sem_NI{j},2))])
        elseif i==2
            disp([type{j},':' num2str(round(median(aP_All{j}),2)),' +/- ' num2str(round(sem_aP{j},2))])
        elseif i==3
            disp([type{j},':' num2str(round(median(sP_All{j}),2)),' +/- ' num2str(round(sem_sP{j},2))])
        elseif i==4
            disp([type{j},':' num2str(round(median(eP_All{j}),2)),' +/- ' num2str(round(sem_eP{j},3))])
        end
    end
end

% for i=1:4
%     for j = 1: size(NI{1},2)
%         plot(hPlot.hPlot1(1),1,NI{1}(j),'color','k','marker','o')
%         plot(hPlot.hPlot1(2),1,aP_All{1}(j),'color','k','marker','o')
%         plot(hPlot.hPlot1(3),1,sP_All{1}(j),'color','k','marker','o')
%         plot(hPlot.hPlot1(4),1,eP_All{1}(j),'color','k','marker','o')
% 
%     end
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


title(hPlot.hPlot1(4),'% Explained Variance' ,'fontSize',20);

xLims = get(hPlot.hPlot1(1),'XLim');
pl = line(xLims,[1 1],'parent',hPlot.hPlot1(1));
set(pl,'color','k','LineStyle','-','LineWidth',2);
p2 = line(xLims,[2 2],'parent',hPlot.hPlot1(1));
set(p2,'color','k','LineStyle','--','LineWidth',2);

% stat tests
pairs = nchoosek(1:length(barNames),2);

disp('Statistical tests for aP')
for i=1:size(pairs,1)
    clear x
    x = [aP_All{1,pairs(i,1)}',aP_All{1,pairs(i,2)}'];
    p_aP(i) = kruskalwallis(x,[],'off');
    disp(['p between ' barNames{pairs(i,1)} ' and ' barNames{pairs(i,2)} ': ' num2str(p_aP(i))])
end


disp('Statistical tests for sP')
for i=1:size(pairs,1)
    clear y
    y = [sP_All{1,pairs(i,1)}',sP_All{1,pairs(i,2)}'];
    p_sP(i) = kruskalwallis(y,[],'off');
    disp(['p between ' barNames{pairs(i,1)} ' and ' barNames{pairs(i,2)} ': ' num2str(p_sP(i))])
end

end
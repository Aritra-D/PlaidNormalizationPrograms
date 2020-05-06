

clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aritra - if you make one file for N=0, S=0 (all electrodes) you can
% always sub-select based on N,snr and d.
load('E:\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\all_N15_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');

N=15; snr=2; d=0.75;
combineUniqueElectrodes=1;

cList = [0 6.25 12.5 25 50]; %[0 1 2 4 8]/16;

% 1- Standard Divisive Normalization (Untuned normalization)
% 2- Tuned Normalization
% 3- EMS-stimulus Normalization on actual dataset
% 4- EMS-stimulus Normalization on symmetrized dataset
modelStringList = {'Untuned','Tuned-alpha_c2','Population'};

for iModel =1:3
    modelNum = iModel;
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
    allData = cell(1,4);
    clear allDataTMP
    allDataTMP = zeros(numGoodElectrodes,5,5);
    for i=1:numGoodElectrodes
        clear x
        x = squeeze(firingRateData.analysisDataST(goodIDList{i},1,:,:) - firingRateData.analysisData_cBL(goodIDList{i},1,:,:));
        
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
    allData{1} = allDataTMP; allDataType{1} = 'FR';
    
    allDataType{2} = 'G'; allDataType{3} = 'HG'; allDataType{4} = 'S';
    for i=2:4
        eDataTMP = 10*(squeeze(energyData.analysisDataST{i}) - squeeze(energyData.analysisData_cBL{i}));
        
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
        allData{i} = allDataTMP;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    eP = zeros(4,numGoodElectrodes);
    aP = zeros(4,numGoodElectrodes);
    sP = zeros(4,numGoodElectrodes);
    dataP = zeros(4,numGoodElectrodes,5,5);
    
    for i=1:4
        for j=1:numGoodElectrodes
            
            %Plaid
            z = squeeze(allData{i}(j,:,:));
            %         z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
            parP = getParametersPlaidV2(z,modelNum);
            [dp,pz] = getResponseMatrixPlaidV2(parP,z,modelNum);
            eP(i,j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
            if modelNum==1
                %             aP(i,j) = 0;
                sP(i,j) = parP(3);
            else
                aP(i,j) = parP(3);
                sP(i,j) = parP(4);
            end
            dataP(i,j,:,:) = z;
            
        end
    end
    
    eP_allModels{iModel} = eP; %#ok<*SAGROW>
    sP_allModels{iModel} = sP;
    
    if modelNum==1
    else
        aP_allModels{iModel-1} = aP;
    end
    
    titleString = {'Spikes','Gamma','Hi-Gamma','SSVEP'};
    paramString = {'N.I.','eP','aP','sP'};
    
    meanP = squeeze(mean(dataP,2));
    hFig=figure;
    set(hFig,'units','normalized','outerposition',[0 0 1 1])
    
    for i=1:4
        %     meanP = mean(squeeze(dataP(i,:,:,:)),1);
        
        subplot(3,4,i);
        imagesc(squeeze(meanP(i,:,:))); colorbar;
        %     imagesc(squeeze(meanP)); colorbar;
        tickLengthPlot = 2*get(subplot(3,4,i),'TickLength');
        
        title(titleString{i});
        
        mp = squeeze(meanP(i,:,:));
        %       mp = squeeze(meanP);
        
        parMP = getParametersPlaidV2(mp,modelNum);
        [dp,pmp] = getResponseMatrixPlaidV2(parMP,mp,modelNum);
        expVar = 1 - (dp/sum((mp(:)-mean(mp(:))).^2));
        
        subplot(3,4,4+i);
        plot(cList,mp(5,:),'color',[0.4 0.4 0.4],'marker','o','linestyle','none','LineWidth',1.5);
        hold on;
        plot(cList,[mp(5,1) mp(4,2) mp(3,3) mp(2,4) mp(1,5)],'color','k','marker','v','linestyle','none','LineWidth',1.5);
        plot(cList,[pmp(5,1) pmp(4,2) pmp(3,3) pmp(2,4) pmp(1,5)],'color','k','LineWidth',1.5);
        plot(cList,pmp(5,:),'color',[0.4 0.4 0.4],'LineWidth',1.5);
        plot(cList,flip(mp(:,1)),'color',[0.6 0.6 0.6],'marker','s','linestyle','none','LineWidth',1.5);
        plot(cList,flip(pmp(:,1)),'color',[0.6 0.6 0.6],'LineWidth',1.5);
        set(subplot(3,4,4+i),'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(subplot(3,4,4+i),'XTick',cList,'XTickLabelRotation',90,'XTickLabel',cList);
        
        if i==1
            xlabel('Contrast(%)','fontSize',14);
            ylabel(subplot(3,4,4+i),{'Change in'  'Spike Rate' '(spike/s)'},'fontSize',14);
        elseif i==2
            ylabel(subplot(3,4,4+i),'Change in Power (dB)','fontSize',14);
            
        end
        
        if modelNum==1
            title({titleString{i}, ['\sigma=' num2str(parMP(3),3) ',ExpVar=' num2str(round(100*expVar)) '%']});
        else
            title({titleString{i}, ['\alpha=' num2str(parMP(3),3) ', \sigma=' num2str(parMP(4),3) ',ExpVar=' num2str(round(100*expVar)) '%']});
        end
        
    end
    
    text(0.65,0.4,'Grating 1','color',[0.4 0.4 0.4],'fontWeight','bold','fontSize',11,'unit','normalized','parent',subplot(3,4,5))
    text(0.65,0.25,'Grating 2','color',[0.6 0.6 0.6],'fontWeight','bold','fontSize',11,'unit','normalized','parent',subplot(3,4,5))
    text(0.65,0.1,'Plaid','color','k','fontWeight','bold','fontSize',11,'unit','normalized','parent',subplot(3,4,5))
    for i=1:4
        xlim(subplot(3,4,4+i),[0 50]);
    end
    subplot(3,4,9)
    NI = squeeze((dataP(:,:,1,1)+dataP(:,:,5,5))./dataP(:,:,1,5));
    bar(median(NI,2));
    hold on;
    errorbar(1:4,median(NI,2),getSEMedian(NI'),'k.');
    
    subplot(3,4,10)
    bar(median(eP,2));
    hold on;
    errorbar(1:4,median(eP,2),getSEMedian(eP'),'k.');
    
    subplot(3,4,11)
    if modelNum==1
    else
        bar(median(aP,2));
        hold on;
        errorbar(1:4,median(aP,2),getSEMedian(aP'),'k.');
    end
    
    subplot(3,4,12)
    bar(median(sP,2));
    hold on;
    errorbar(1:4,median(sP,2),getSEMedian(sP'),'k.');
        params = {'N.I.','aP','sP','eP'};

    % Statistics
    for i = 1:4
        disp([params{i},':'])
        disp('median +/- SEMedian')
        for j = 1:4
            if i==1
                disp([titleString{j},':' num2str(round(median(NI(j,:)),2)),' +/- ' num2str(round(getSEMedian(NI(j,:)'),2))])
            elseif i==2
                if modelNum==2 || modelNum==3
                    disp([titleString{j},':' num2str(round(median(aP(j,:)),2)),' +/- ' num2str(round(getSEMedian(aP(j,:)'),2))])
                end
            elseif i==3
                disp([titleString{j},':' num2str(round(median(sP(j,:)),2)),' +/- ' num2str(round(getSEMedian(sP(j,:)'),2))])
            elseif i==4
                disp([titleString{j},':' num2str(round(median(eP(j,:)),2)),' +/- ' num2str(round(getSEMedian(eP(j,:)'),3))])
            end
        end
    end
    
    % stat tests
    disp(['Model: ' iModel ' ' modelStringList{iModel}])
    pairs = nchoosek(1:length(titleString),2);
    
    if modelNum==2||modelNum==3
        
        disp('Statistical tests for aP- Kruskal-Wallis test')
        for i=1:size(pairs,1)
            clear statData
            statData = [aP(pairs(i,1),:)',aP(pairs(i,2),:)'];
            p_aP(i) = kruskalwallis(statData,[],'off');
            disp(['p between ' titleString{pairs(i,1)} ' and ' titleString{pairs(i,2)} ': ' num2str(p_aP(i))])
        end
    end
    
    disp('Statistical tests for sP- Kruskal-Wallis test')
    for i=1:size(pairs,1)
        clear statData2
        statData2 = [sP(pairs(i,1),:)',sP(pairs(i,2),:)'];
        p_sP(i) = kruskalwallis(statData2,[],'off');
        disp(['p between ' titleString{pairs(i,1)} ' and ' titleString{pairs(i,2)} ': ' num2str(p_sP(i))])
    end
    
    if modelNum==3 %KS Test for NI populations
        disp('Statistical tests for NI- Kruskal-Wallis test')
        for i=1:size(pairs,1)
            clear statData3
            statData3 = [NI(pairs(i,1),:)',NI(pairs(i,2),:)'];
            p_NI(i) = kruskalwallis(statData3,[],'off');
%             statData3 = NI(pairs(i,1),:)';
%             statData4 = NI(pairs(i,2),:)';
%             [h_NI(i),p_NI(i)] = kstest2(statData3,statData4);
            disp(['p between ' titleString{pairs(i,1)} ' and ' titleString{pairs(i,2)} ': ' num2str(p_NI(i))])
        end      
    end
    % PlotLabels
    for i =1:4
        subplot(3,4,8+i);
        title(paramString{i});
    end
    folderSave_Figs = fullfile('E:\Projects\Aritra_PlaidNormalizationProject\Figures\fitting');
    if ~exist(folderSave_Figs,'dir')
        mkdir(folderSave_Figs)
    end
    FigName = fullfile(folderSave_Figs,['N15_S2_d_0.75_T150_400_ModelFitResults_' num2str(modelNum) '_' modelStringList{modelNum}]);
    saveas(hFig,[FigName '.fig'])
    saveas(hFig,[FigName '.tif'])
end

for i=1:2
    aP_median{i} =  median(aP_allModels{i},2)';
    for j = 1:4
        aP_seMedian{i}(j) = getSEMedian(aP_allModels{i}(j,:));
    end
end

for i=1:4
    aPBar_median{i} = [aP_median{1}(i) aP_median{2}(i)];
    aPBar_seMedian{i} = [aP_seMedian{1}(i) aP_seMedian{2}(i)];
end

for i=1:3
    sP_median{i} =  median(sP_allModels{i},2)';
    eP_median{i} =  median(eP_allModels{i},2)';
    for j = 1:4
        sP_seMedian{i}(j) = getSEMedian(sP_allModels{i}(j,:));
        eP_seMedian{i}(j) = getSEMedian(eP_allModels{i}(j,:));
    end
end

for i=1:4
    sPBar_median{i} = [sP_median{1}(i) sP_median{2}(i) sP_median{3}(i)];
    sPBar_seMedian{i} = [sP_seMedian{1}(i) sP_seMedian{2}(i) sP_seMedian{3}(i)];
    ePBar_median{i} = [eP_median{1}(i) eP_median{2}(i) eP_median{3}(i)];
    ePBar_seMedian{i} = [eP_seMedian{1}(i) eP_seMedian{2}(i) eP_seMedian{3}(i)];
end

aP_Bar =[];aP_ErrorBar=[];
for i=1:4
    aP_Bar = cat(1,aP_Bar,aPBar_median{i});
    aP_ErrorBar = cat(1,aP_ErrorBar,aPBar_seMedian{i});
end
sP_Bar =[];sP_ErrorBar=[];
for i=1:4
    sP_Bar = cat(1,sP_Bar,sPBar_median{i});
    sP_ErrorBar = cat(1,sP_ErrorBar,sPBar_seMedian{i});
end
eP_Bar =[];eP_ErrorBar=[];
for i=1:4
    eP_Bar = cat(1,eP_Bar,ePBar_median{i});
    eP_ErrorBar = cat(1,eP_ErrorBar,ePBar_seMedian{i});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Figure 4   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hFig2 = figure;
set(hFig2,'units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
b = bar(median(NI,2));
b.FaceColor = [0 0.4470 0.7410];
hold on;
errorbar(1:4,median(NI,2),getSEMedian(NI'),'k.','lineWidth',1);
set(subplot(2,2,1),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

xLims = get(subplot(2,2,1),'XLim');
pl = line(xLims,[1 1],'parent',subplot(2,2,1));
set(pl,'color','k','LineStyle','-','LineWidth',2);
p2 = line(xLims,[2 2],'parent',subplot(2,2,1));
set(p2,'color','k','LineStyle','--','LineWidth',2);


subplot(2,2,3)
b1 = bar(aP_Bar, 'grouped');
b1(1).FaceColor = [0 0.65 0.65];
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(aP_Bar, 1);
nbars = size(aP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, aP_Bar(:,i), aP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(2,2,3),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
legend(subplot(2,2,3),{'Tuned','Population'},'fontSize',14,'Location','best')


subplot(2,2,2)
b2 = bar(sP_Bar, 'grouped');
b2(2).FaceColor = [0 0.65 0.65];

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(sP_Bar, 1);
nbars = size(sP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, sP_Bar(:,i), sP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(2,2,2),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
legend(subplot(2,2,2),{'Untuned','Tuned','Population'},'fontSize',14,'Location','best')

subplot(2,2,4)
b3 = bar(eP_Bar, 'grouped');
b3(2).FaceColor = [0 0.65 0.65];

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(eP_Bar, 1);
nbars = size(eP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, eP_Bar(:,i), eP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(2,2,4),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

title(subplot(2,2,1),'N.I.','fontSize',14);
title(subplot(2,2,3),'\alpha','fontSize',14);
title(subplot(2,2,2),'\sigma','fontSize',14);
title(subplot(2,2,4),'% Explained Variance' ,'fontSize',14);
for i=1:4
    ylabel(subplot(2,2,i),'Population Median')
    set(subplot(2,2,i),'xticklabel',titleString)
end

% Break axis for subplot (2,2,3)
% breakInfo = breakyaxis(subplot(2,2,3),[4 7]);
% set(breakInfo.lowAxes,'YLim',[0 1.3]);
% set(breakInfo.highAxes,'YLim',[5 10]);
% legend(breakInfo.highAxes,{'Tuned','Population'},'color','w','fontSize',14,'Location','best')
% line('yData',breakInfo.yPoints1(1:10),'xdata',breakInfo.xPoints(1:10),'Parent',breakInfo.breakAxes,'Color','k','LineWidth',1);
% line('yData',breakInfo.yPoints2(1:10),'xdata',breakInfo.xPoints(1:10),'Parent',breakInfo.breakAxes,'Color','k','LineWidth',1);

textH{1} = getPlotHandles(1,1,[0.07 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.07 0.45 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.5 0.95 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.5 0.45 0.01 0.01]);

textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
FigName2 = fullfile(folderSave_Figs,'N15_S2_d_0.75_T150_400_ModelFitResults_allModelsComparison');
saveas(hFig2,[FigName2 '.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Figure 5   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig3 = figure;
set(hFig3,'units','normalized','outerposition',[0 0 1 1])

subplot(1,3,1)
b1 = bar(aP_Bar, 'grouped');
b1(1).FaceColor = [0 0.65 0.65];
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(aP_Bar, 1);
nbars = size(aP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, aP_Bar(:,i), aP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(1,3,1),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
% legend(subplot(1,3,1),{'Tuned','Population'},'fontSize',14,'Location','best')
% Break axis for subplot (2,2,3)

subplot(1,3,2)
b2 = bar(sP_Bar, 'grouped');
b2(2).FaceColor = [0 0.65 0.65];

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(sP_Bar, 1);
nbars = size(sP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, sP_Bar(:,i), sP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(1,3,2),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
legend(subplot(1,3,2),{'Untuned','Tuned','Population'},'fontSize',14,'Location','best')

subplot(1,3,3)
b3 = bar(eP_Bar, 'grouped');
b3(2).FaceColor = [0 0.65 0.65];

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(eP_Bar, 1);
nbars = size(eP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, eP_Bar(:,i), eP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(1,3,3),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

title(subplot(1,3,1),'\alpha','fontSize',14);
title(subplot(1,3,2),'\sigma','fontSize',14);
title(subplot(1,3,3),'% Explained Variance' ,'fontSize',14);

for i=1:3
    ylabel(subplot(1,3,i),'Population Median')
    set(subplot(1,3,i),'xticklabel',titleString,'XTickLabelRotation',45,'fontSize',14)
end
% breakInfo2 = breakyaxis(subplot(1,3,1),[4 7]);
% set(breakInfo2.lowAxes,'YLim',[0 1.3]);
% set(breakInfo2.highAxes,'YLim',[6 10]);
% set(breakInfo2.highAxes,'YTick',7:10);
% legend(breakInfo2.highAxes,{'Tuned','Population'},'color','w','fontSize',14,'Location','best')
% line('yData',breakInfo2.yPoints1(1:16),'xdata',breakInfo2.xPoints(1:16),'Parent',breakInfo2.breakAxes,'Color','k','LineWidth',1);
% line('yData',breakInfo2.yPoints2(1:16),'xdata',breakInfo2.xPoints(1:16),'Parent',breakInfo2.breakAxes,'Color','k','LineWidth',1);
% 

textH{1} = getPlotHandles(1,1,[0.07 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.36 0.95 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.64 0.95 0.01 0.01]);


textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
FigName3 = fullfile(folderSave_Figs,'N15_S2_d_0.75_T150_400_ModelFitResults_allModelsComparison_no_NI');
saveas(hFig3,[FigName3 '.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Figure 6   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig4 = figure;
set(hFig4,'units','normalized','outerposition',[0 0 1 1])

subplot(1,3,1)
b1 = bar(100*eP_Bar, 'grouped');
b1(3).FaceColor = 'y';
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(eP_Bar, 1);
nbars = size(eP_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, 100*eP_Bar(:,i), 100*eP_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
end
set(subplot(1,3,1),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
legend(subplot(1,3,1),{'Untuned','Tuned','Population'},'fontSize',14,'Location','best')


subplot(1,3,2)
b2 = bar(aP_Bar(:,2));
b2.FaceColor = 'y';
hold on
for i=1:4
    errorbar(subplot(1,3,2),i,aP_Bar(i,2),aP_ErrorBar(i,2),'color','k','lineStyle','-')
end

set(subplot(1,3,2),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

subplot(1,3,3)
b3 = bar(sP_Bar(:,3));
b3.FaceColor = 'y';
hold on
for i=1:4
    errorbar(subplot(1,3,3),i,sP_Bar(i,3),sP_ErrorBar(i,3),'color','k','lineStyle','-')
end

set(subplot(1,3,3),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

title(subplot(1,3,1),'% Explained Variance' ,'fontSize',14);
title(subplot(1,3,2),'\alpha','fontSize',14);
title(subplot(1,3,3),'\sigma','fontSize',14);

for i=1:3
    ylabel(subplot(1,3,i),'Population Median')
    set(subplot(1,3,i),'xticklabel',titleString,'XTickLabelRotation',45,'fontSize',14)
end

textH{1} = getPlotHandles(1,1,[0.07 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.36 0.95 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.64 0.95 0.01 0.01]);


textString = {'A','B','C','D'};
for i = 1:4
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
FigName4 = fullfile(folderSave_Figs,'N15_S2_d_0.75_T150_400_ModelFitResults_allModelsComparison_withoutNI');
saveas(hFig4,[FigName4 '.fig'])
saveas(hFig4,[FigName4 '.tif'])
% print(hFig4,'-dtiff','resolution' '-r300') 
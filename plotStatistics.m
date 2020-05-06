function hFigure5 = plotStatistics(fileSave,colorScheme)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aritra - if you make one file for N=0, S=0 (all electrodes) you can
% always sub-select based on N,snr and d.
% load('E:\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\all_N15_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');
load(fileSave)

N=15; snr=2; d=0.75;
combineUniqueElectrodes=1;

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
            [dp,pz] = getResponseMatrixPlaidV2(parP,z,modelNum); %#ok<*ASGLU>
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
    sP_allModels{iModel} = sP; %#ok<*AGROW>
    
    if modelNum==1
    else
        aP_allModels{iModel-1} = aP;
    end
    
    neuralMeasures = {'Spikes','Gamma','Hi-Gamma','SSVEP'};
    params = {'N.I.','aP','sP','eP'};
    
    
    NI = squeeze((dataP(:,:,1,1)+dataP(:,:,5,5))./dataP(:,:,1,5));
    
    
    % Statistics
    for i = 1:4 % params
        disp([params{i},':'])
        disp('median +/- SEMedian')
        for j = 1:4 % Neural measures
            if i==1
                disp([neuralMeasures{j},':' num2str(round(median(NI(j,:)),2)),' +/- ' num2str(round(getSEMedian(NI(j,:)'),2))])
            elseif i==2
                if modelNum==2 || modelNum==3
                    disp([neuralMeasures{j},':' num2str(round(median(aP(j,:)),2)),' +/- ' num2str(round(getSEMedian(aP(j,:)'),2))])
                end
            elseif i==3
                disp([neuralMeasures{j},':' num2str(round(median(sP(j,:)),2)),' +/- ' num2str(round(getSEMedian(sP(j,:)'),2))])
            elseif i==4
                disp([neuralMeasures{j},':' num2str(round(median(eP(j,:)),2)),' +/- ' num2str(round(getSEMedian(eP(j,:)'),3))])
            end
        end
    end
    
    % stat tests
    disp(['Model: ' iModel ' ' modelStringList{iModel}])
    pairs = nchoosek(1:length(neuralMeasures),2);
    
    if modelNum==2||modelNum==3
        
        disp('Statistical tests for aP- Kruskal-Wallis test')
        for i=1:size(pairs,1)
            clear statData
            statData = [aP(pairs(i,1),:)',aP(pairs(i,2),:)'];
            p_aP(i) = kruskalwallis(statData,[],'off');
            disp(['p between ' neuralMeasures{pairs(i,1)} ' and ' neuralMeasures{pairs(i,2)} ': ' num2str(p_aP(i))])
        end
    end
    
    disp('Statistical tests for sP- Kruskal-Wallis test')
    for i=1:size(pairs,1)
        clear statData2
        statData2 = [sP(pairs(i,1),:)',sP(pairs(i,2),:)'];
        p_sP(i) = kruskalwallis(statData2,[],'off');
        disp(['p between ' neuralMeasures{pairs(i,1)} ' and ' neuralMeasures{pairs(i,2)} ': ' num2str(p_sP(i))])
    end
    
    if modelNum==3 %KW Test for NI populations
        disp('Statistical tests for NI- Kruskal-Wallis test')
        for i=1:size(pairs,1)
            clear statData3
            statData3 = [NI(pairs(i,1),:)',NI(pairs(i,2),:)'];
            p_NI(i) = kruskalwallis(statData3,[],'off');
            %             statData3 = NI(pairs(i,1),:)';
            %             statData4 = NI(pairs(i,2),:)';
            %             [h_NI(i),p_NI(i)] = kstest2(statData3,statData4);
            disp(['p between ' neuralMeasures{pairs(i,1)} ' and ' neuralMeasures{pairs(i,2)} ': ' num2str(p_NI(i))])
        end
    end
end

% calculating aP_median and aP_seMedian for two normalization (tuned &
% population models)
for i=1:2 % models
    aP_median{i} =  median(aP_allModels{i},2)';
    for j = 1:4 % Neural measures
        aP_seMedian{i}(j) = getSEMedian(aP_allModels{i}(j,:));
    end
end

% calculating aP_median and aP_seMedian for all normalization (untuned,
% tuned & population models)
for i=1:3 % models
    sP_median{i} =  median(sP_allModels{i},2)';
    eP_median{i} =  median(eP_allModels{i},2)';
    for j = 1:4 % Neural measures
        sP_seMedian{i}(j) = getSEMedian(sP_allModels{i}(j,:));
        eP_seMedian{i}(j) = getSEMedian(eP_allModels{i}(j,:));
    end
end

% grouping parameter values according to neural measures
for i=1:4 % Neural Measures
    aPBar_median{i} = [aP_median{1}(i) aP_median{2}(i)];
    aPBar_seMedian{i} = [aP_seMedian{1}(i) aP_seMedian{2}(i)];
    sPBar_median{i} = [sP_median{1}(i) sP_median{2}(i) sP_median{3}(i)];
    sPBar_seMedian{i} = [sP_seMedian{1}(i) sP_seMedian{2}(i) sP_seMedian{3}(i)];
    ePBar_median{i} = [eP_median{1}(i) eP_median{2}(i) eP_median{3}(i)];
    ePBar_seMedian{i} = [eP_seMedian{1}(i) eP_seMedian{2}(i) eP_seMedian{3}(i)];
end

% converting into arrays from cells that are used to plot barplots
aP_Bar =[];aP_ErrorBar=[];
sP_Bar =[];sP_ErrorBar=[];
eP_Bar =[];eP_ErrorBar=[];
for i=1:4
    aP_Bar = cat(1,aP_Bar,aPBar_median{i});
    aP_ErrorBar = cat(1,aP_ErrorBar,aPBar_seMedian{i});
    
    sP_Bar = cat(1,sP_Bar,sPBar_median{i});
    sP_ErrorBar = cat(1,sP_ErrorBar,sPBar_seMedian{i});
    
    eP_Bar = cat(1,eP_Bar,ePBar_median{i});
    eP_ErrorBar = cat(1,eP_ErrorBar,ePBar_seMedian{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Figure 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFigure5 = figure(5);
set(hFigure5,'units','normalized','outerposition',[0 0 1 1])

subplot(1,3,1)
b1 = bar(100*eP_Bar, 'grouped');

if strcmp(colorScheme,'color')
    b1(3).FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    colors = repmat(0.65:-0.1:0.45,[3 1])';
    for i = 1:size(eP_Bar,2)
        b1(i).FaceColor = colors(i,:,:);
    end
end
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
tickLengthPlot = 2*get(subplot(1,3,1),'TickLength');

set(subplot(1,3,1),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
legend(subplot(1,3,1),{'Untuned','Tuned','Population'},'fontSize',14,'Location','best')


subplot(1,3,2)
b2 = bar(aP_Bar(:,2)); % population model (aP is not present in untuned normalization model)
if strcmp(colorScheme,'color')
    b2.FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    b2.FaceColor = [0.45 0.45 0.45];
end
hold on
for i=1:4% Neural measures
    errorbar(subplot(1,3,2),i,aP_Bar(i,2),aP_ErrorBar(i,2),'color','k','lineStyle','-')
end

set(subplot(1,3,2),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

subplot(1,3,3)
b3 = bar(sP_Bar(:,3)); % population model
if strcmp(colorScheme,'color')
    b3.FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    b3.FaceColor = [0.45 0.45 0.45];
end
hold on
for i=1:4 % Neural measures
    errorbar(subplot(1,3,3),i,sP_Bar(i,3),sP_ErrorBar(i,3),'color','k','lineStyle','-')
end

set(subplot(1,3,3),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

title(subplot(1,3,1),'% Explained Variance' ,'fontSize',14);
title(subplot(1,3,2),'\alpha','fontSize',14);
title(subplot(1,3,3),'\sigma','fontSize',14);

textH{1} = getPlotHandles(1,1,[0.07 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.36 0.95 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.64 0.95 0.01 0.01]);
textString = {'A','B','C'};

for i=1:3 % Plot numbers
    ylabel(subplot(1,3,i),'Population Median')
    set(subplot(1,3,i),'xticklabel',neuralMeasures,'XTickLabelRotation',45,'fontSize',14)
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end

end

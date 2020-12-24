function hFigure = plotStatistics(monkeyName,folderSourceString,stimType,elecParams,timeRangeForComputation,tapers_MT,combineUniqueElectrodeData,colorScheme)

close all; % uncomment if this function is run separately! You will also
% need elecParams from plotFigures
if ~exist(folderSourceString,'var'); folderSourceString = 'E:\'; end

if strcmp(stimType,'Static')
    folderData = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures');
elseif strcmp(stimType,'Counterphase')
    folderData = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures\Counterphase');
end

gridType = 'microelectrode';
removeERPFlag = 0;
% data is saved for all non-unique electrodes; cne is set as 0;
dataParams = ['_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(0) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) '.mat'];

N_Spike=elecParams.spikeCutoff; snr=elecParams.snrCutoff; d=elecParams.dRange(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aritra - if you make one file for N=0, S=0 (all electrodes) you can
% always sub-select based on N,snr and d.
fileSave = fullfile(folderData,'',[monkeyName dataParams]);
disp(['Loading File: ' fileSave])
load(fileSave) %#ok<*LOAD>

%%%%%%%%%%%%%%%% Subselect electrodes based on N, snr and d %%%%%%%%%%%%%%%
goodN = (max([elecInfo.N(1,:);elecInfo.N(2,:)])>N_Spike);
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
if combineUniqueElectrodeData
    
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
        %         z = flip(y,1); z = (z+z')/2; z = flip(z,1); % provisions to make data symmetric
        allDataTMP(i,:,:) = x;
    else
        xs = zeros(size(x));
        for j=1:length(goodIDList{i})
            %             y = squeeze(x(j,:,:));
            %             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % provisions to make data symmetric
            xs(j,:,:) = squeeze(x(j,:,:));
        end
        allDataTMP(i,:,:) = squeeze(mean(xs,1));
    end
end

allData{1} = allDataTMP;
allDataType{1} = 'FR'; allDataType{2} = 'G'; allDataType{3} = 'HG'; allDataType{4} = 'S';

for i=2:4
    eDataTMP = 10*(squeeze(energyData.analysisDataST{i}) - squeeze(energyData.analysisData_cBL{i})); % in dB
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1- Standard Divisive Normalization (Untuned normalization)
% 2- Tuned Normalization
% 3- EMS-stimulus Normalization on actual dataset
% 4- EMS-stimulus Normalization on symmetrized dataset (not used in this analysis)

models = {'Untuned','Tuned-alpha_c2','Population'};
neuralMeasures = {'Spikes','Gamma','Hi-Gamma','SSVEP'};
params = {'N.I.','aP','sP','eP','AIC'};
modelPairs = nchoosek(1:length(models),2);
neuralMeasurePairs = nchoosek(1:length(neuralMeasures),2);


for iModel =1:3
    modelNum = iModel;
    clear p_aP p_sP p_AIC p_NI
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    eP = zeros(4,numGoodElectrodes);
    aP = zeros(4,numGoodElectrodes);
    sP = zeros(4,numGoodElectrodes);
    dataP = zeros(4,numGoodElectrodes,5,5);
    rSS = zeros(4,numGoodElectrodes);
    aicVals = zeros(4,numGoodElectrodes);
    
    for i=1:4 % Neural Measures
        for j=1:numGoodElectrodes
            clear z
            %Plaid
            z = squeeze(allData{i}(j,:,:));
            %         z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
            parP = getParametersPlaidV2(z,modelNum);
            [dp,pz] = getResponseMatrixPlaidV2(parP,z,modelNum); %#ok<*ASGLU>
            eP(i,j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
            rSS(i,j) = dp; % storing dps as residual squared sum for estimation of maximum likelihood for each model
            aicVals(i,j) = calculateAIC(dp,numel(z),size(parP,2));
            
            if modelNum==1
                sP(i,j) = min(max(0,parP(3)),5);
            else
                aP(i,j) = parP(3);
                sP(i,j) = min(max(0,parP(4)),5);
            end
            dataP(i,j,:,:) = z;
        end
    end
    
    if modelNum==2 || modelNum==3
        aP_Models{modelNum-1} = aP; %#ok<*NASGU>
    end
    sP_Models{iModel} = sP;
    eP_Models{iModel} = eP;
    rSS_Models{iModel} = rSS;
    nObs_Models{iModel} = numel(z);
    params_Models{iModel} = size(parP,2);
    aic_Models{iModel} = aicVals;
    
    % Statistics
    disp(['Monkey: ' monkeyName])
    disp(['Model: ' iModel ' ' models{iModel}])
    
    for i = 1:5 % params
        disp([params{i},': median +/- SEMedian'])
        for j = 1:4 % Neural measures
            clear aP_median_temp aP_seMedian_temp 
            clear sP_median_temp sP_seMedian_temp
            clear eP_median_temp eP_seMedian_temp
            clear aic_median_temp aic_seMedian_temp
            
            if i==1
                if modelNum==3
                    NI = squeeze((dataP(:,:,1,1)+dataP(:,:,5,5))./dataP(:,:,1,5));
                    disp([neuralMeasures{j},':' num2str(round(median(NI(j,:)),3)),' +/- ' num2str(round(getSEMedian(NI(j,:)',1000),2))])
                end
                
            elseif i==2
                if modelNum==2 || modelNum==3
                    aP_median_temp = median(aP(j,:));
                    aP_seMedian_temp = getSEMedian(aP(j,:)',1000);
                    disp([neuralMeasures{j},':' num2str(round(aP_median_temp,2)),' +/- ' num2str(round(aP_seMedian_temp,2))])
                    aP_median{modelNum-1}(j) = aP_median_temp;
                    aP_seMedian{modelNum-1}(j)= aP_seMedian_temp;
                end
                
            elseif i==3
                sP_median_temp = median(sP(j,:));
                sP_seMedian_temp = getSEMedian(sP(j,:)',1000);
                disp([neuralMeasures{j},':' num2str(round(sP_median_temp,2)),' +/- ' num2str(round(sP_seMedian_temp,2))])
                sP_median{iModel}(j) = sP_median_temp;
                sP_seMedian{iModel}(j)= sP_seMedian_temp;
            
            elseif i==4
                eP_median_temp = median(eP(j,:));
                eP_seMedian_temp = getSEMedian(eP(j,:)',1000);
                disp([neuralMeasures{j},':' num2str(round(eP_median_temp,2)),' +/- ' num2str(round(eP_seMedian_temp,3))])
                eP_median{iModel}(j) = eP_median_temp;
                eP_seMedian{iModel}(j)= eP_seMedian_temp;
            
            elseif i==5
                aic_median_temp = median(aicVals(j,:));
                aic_seMedian_temp = getSEMedian(aicVals(j,:)',1000);
                disp([neuralMeasures{j},':' num2str(round(aic_median_temp,2)),' +/- ' num2str(round(aic_seMedian_temp,3))])
                aic_median{iModel}(j) = aic_median_temp;
                aic_seMedian{iModel}(j)= aic_seMedian_temp;
            end
        end
    end
    
    % stat tests
    disp(['Model: ' iModel ' ' models{iModel}])
    if modelNum==2||modelNum==3
        disp('Statistical tests for aP- Kruskal-Wallis test')
        for i=1:size(neuralMeasurePairs,1)
            clear statData
            statData = [aP(neuralMeasurePairs(i,1),:)',aP(neuralMeasurePairs(i,2),:)'];
            p_aP(i) = kruskalwallis(statData,[],'off');
            disp(['p between ' neuralMeasures{neuralMeasurePairs(i,1)} ' and ' neuralMeasures{neuralMeasurePairs(i,2)} ': ' num2str(p_aP(i))])
        end
    end
    
    disp('Statistical tests for sP- Kruskal-Wallis test')
    for i=1:size(neuralMeasurePairs,1)
        clear statData2
        statData2 = [sP(neuralMeasurePairs(i,1),:)',sP(neuralMeasurePairs(i,2),:)'];
        p_sP(i) = kruskalwallis(statData2,[],'off');
        disp(['p between ' neuralMeasures{neuralMeasurePairs(i,1)} ' and ' neuralMeasures{neuralMeasurePairs(i,2)} ': ' num2str(p_sP(i))])
    end
    
    if modelNum==3 %KW Test for NI populations
        disp('Statistical tests for NI- Kruskal-Wallis test')
        for i=1:size(neuralMeasurePairs,1)
            clear statData3
            statData3 = [NI(neuralMeasurePairs(i,1),:)',NI(neuralMeasurePairs(i,2),:)'];
            p_NI(i) = kruskalwallis(statData3,[],'off');
            disp(['p between ' neuralMeasures{neuralMeasurePairs(i,1)} ' and ' neuralMeasures{neuralMeasurePairs(i,2)} ': ' num2str(p_NI(i))])
        end
        
        % Spearman's Correlation analysis between N.I. values of different Neural
        % Measures
        for i=1:size(neuralMeasurePairs,1)
            clear rho pVals_Corr
            [rho,pVals_Corr] = corr(NI(neuralMeasurePairs(i,1),:)',NI(neuralMeasurePairs(i,2),:)','Type','Spearman');
            disp([neuralMeasures{neuralMeasurePairs(i,1)} ' and ' neuralMeasures{neuralMeasurePairs(i,2)} ': rho: ' num2str(rho,2) ', p_Corr: ' num2str(pVals_Corr,2)])
        end
    end
end

% Significance Test across normalization models for each neural measure
disp('Statistical tests for AIC- Kruskal-Wallis test')
for i=1:length(neuralMeasures)
    for j= 1:size(modelPairs,1)
        clear statData4
        statData4 = [aic_Models{modelPairs(j,1)}(i,:)',aic_Models{modelPairs(j,2)}(i,:)'];
        p_AIC(i,j) = kruskalwallis(statData4,[],'off');
        disp([neuralMeasures{i} ' : p between ' models{modelPairs(j,1)} ' and ' models{modelPairs(j,2)} ' : ' num2str(p_AIC(i,j))])
    end
end


disp('Statistical tests for Wilcoxon Rank Sum test')
for i=1:length(neuralMeasures)
    for j= 1:size(modelPairs,1)
        [p_AIC_Wil(i,j),h_AIC_Wil(i,j),stats_AIC_wil{i,j}] = ranksum(aic_Models{modelPairs(j,1)}(i,:),aic_Models{modelPairs(j,2)}(i,:),'method','approximate');
        disp([neuralMeasures{i} ' : p between ' models{modelPairs(j,1)} ' and ' models{modelPairs(j,2)} ' : ' num2str(p_AIC_Wil(i,j))])
    end
end

% % calculating aP_median and aP_seMedian for two normalization (tuned &
% % population models)
% for i=1:2 % models
%     aP_median{i} =  median(aP_allModels{i},2)';
%     for j = 1:4 % Neural measures
%         aP_seMedian{i}(j) = getSEMedian(aP_allModels{i}(j,:));
%     end
% end
%
% % calculating eP/sP_median and eP/sP_seMedian for all normalization (untuned,
% % tuned & population models)
% for i=1:3 % models
%     sP_median{i} =  median(sP_allModels{i},2)';
%     eP_median{i} =  median(eP_allModels{i},2)';
%     for j = 1:4 % Neural measures
%         sP_seMedian{i}(j) = getSEMedian(sP_allModels{i}(j,:));
%         eP_seMedian{i}(j) = getSEMedian(eP_allModels{i}(j,:));
%     end
% end

% % statistical distribution and Krushkal-Wallis test of AIC values of normalization models
% disp('AIC for Normalization models:')
% disp('all monkeys:' )
% for i = 1:4 % Neural measures
%     for iModel =1:3 % Models
%         aicVals_median{iModel} = median(aic_Models{iModel}(i,:));
%         aicVals_seMedian{iModel} = getSEMedian(aic_Models{iModel}(i,:)',1000);
%         disp([neuralMeasures{i}, ' -- Model: ' iModel ': ' models{iModel} ': median +/- SEMedian: ' num2str(round(aicVals_median,2)),' +/- ' num2str(round(aicVals_seMedian,2))])
%
%
%         if iModel == 3
%             models = 1:1:length(models);
%             aic_grouped{i}= [aic_Models{models(1)}(i,:)', aic_Models{models(2)}(i,:)', aic_Models{models(3)}(i,:)'];
%             p_aic(i) = kruskalwallis(aic_grouped{i},[],'off');
%             disp(['Krushkal Wallis test, data compared to all models: p =  ' num2str(p_aic(i))])
%
%             modelPairs = nchoosek(1:length(models),2);
%             for j=1:size(modelPairs,1)
%                 aic_groupedModelPairs{i}{j} = [aic_Models{modelPairs(j,1)}(i,:)', aic_Models{modelPairs(j,2)}(i,:)'];
%                 p_aic_ModelPairs(i,j) = kruskalwallis(aic_groupedModelPairs{i}{j},[],'off');
%                 disp(['Krushkal Wallis test, data compared between Model ' num2str(modelPairs(j,1)) ' and Model ' num2str(modelPairs(j,2)) ': p =  ' num2str(p_aic_ModelPairs(i,j))])
%
%             end
%         end
%     end
% end



% grouping parameter values according to neural measures
for i=1:4 % Neural Measures
    aPBar_median{i} = [aP_median{1}(i) aP_median{2}(i)]; %#ok<*AGROW>
    aPBar_seMedian{i} = [aP_seMedian{1}(i) aP_seMedian{2}(i)];
    sPBar_median{i} = [sP_median{1}(i) sP_median{2}(i) sP_median{3}(i)];
    sPBar_seMedian{i} = [sP_seMedian{1}(i) sP_seMedian{2}(i) sP_seMedian{3}(i)];
    ePBar_median{i} = [eP_median{1}(i) eP_median{2}(i) eP_median{3}(i)];
    ePBar_seMedian{i} = [eP_seMedian{1}(i) eP_seMedian{2}(i) eP_seMedian{3}(i)];
    aicBar_median{i} = [aic_median{1}(i) aic_median{2}(i) aic_median{3}(i)];
    aicBar_seMedian{i} = [aic_seMedian{1}(i) aic_seMedian{2}(i) aic_seMedian{3}(i)];
end

% converting into arrays from cells that are used to plot barplots
aP_Bar =[];aP_ErrorBar=[];
sP_Bar =[];sP_ErrorBar=[];
eP_Bar =[];eP_ErrorBar=[];
aic_Bar = []; aic_ErrorBar = [];
for i=1:4
    aP_Bar = cat(1,aP_Bar,aPBar_median{i});
    aP_ErrorBar = cat(1,aP_ErrorBar,aPBar_seMedian{i});
    
    sP_Bar = cat(1,sP_Bar,sPBar_median{i});
    sP_ErrorBar = cat(1,sP_ErrorBar,sPBar_seMedian{i});
    
    eP_Bar = cat(1,eP_Bar,ePBar_median{i});
    eP_ErrorBar = cat(1,eP_ErrorBar,ePBar_seMedian{i});
    
    aic_Bar = cat(1,aic_Bar,aicBar_median{i});
    aic_ErrorBar = cat(1,aic_ErrorBar,aicBar_seMedian{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Figure 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFigure = figure(5);
set(hFigure,'units','normalized','outerposition',[0 0 1 1])

jitterAmount = 0.03;
jitterAmount2 = 0.3;
%%%%%%%%%%%%%%%%%%% (A) BarPlots of Explained Variance %%%%%%%%%%%%%%%%%%%%

subplot(2,2,1)
b1 = bar(100*eP_Bar, 'grouped');

if strcmp(colorScheme,'color')
    b1(3).FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    colorVals = flip(gray(4)); colorVals(4,:) = [];
    colors = colorVals;
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
    %     for j = 1:ngroups
    %         scatter(x(j)*ones(1,size(eP_Models{i},2)),100*eP_Models{i}(j,:),'k','filled','jitter','on','jitterAmount',jitterAmount)
    %     end
end

tickLengthPlot = 3*get(subplot(2,2,1),'TickLength');
yTicks = 0:25:100;
set(subplot(2,2,1),'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
ylim(subplot(2,2,1),[0 100]);
legend(subplot(2,2,1),{'Standard','Tuned','Population'},'fontSize',14,'Location','best')


%%%%%%%%%%%%%%%%%%% (B) BarPlots of AIC values %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
b2 = bar(aic_Bar, 'grouped');

if strcmp(colorScheme,'color')
    b2.FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    colorVals = flip(gray(4)); colorVals(4,:) = [];
    colors = colorVals;
    for i = 1:size(aic_Bar,2)
        b2(i).FaceColor = colors(i,:,:);
    end
end
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(aic_Bar, 1);
nbars = size(aic_Bar, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, aic_Bar(:,i), aic_ErrorBar(:,i), 'k', 'linestyle', 'none','lineWidth',1);
    %     for j = 1:ngroups
    %         scatter(x(j)*ones(1,size(aic_Models{i},2)),100*aic_Models{i}(j,:),'k','filled','jitter','on','jitterAmount',jitterAmount)
    %     end
end

set(subplot(2,2,2),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
% ylim(subplot(2,2,2),[0 100]);
legend(subplot(2,2,2),{'Standard','Tuned','Population'},'fontSize',14,'Location','best')

%%%%%%%%%%%%%%%%%%% (C) BarPlots of Tuned Normalization parameters for all Neural Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
b3 = bar(aP_Bar(:,2)); % population model (aP is not present in untuned normalization model)
if strcmp(colorScheme,'color')
    b3.FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    b3.FaceColor = [0.33 0.33 0.33];
end
hold on
for i=1:4% Neural measures
    errorbar(subplot(2,2,3),i,aP_Bar(i,2),aP_ErrorBar(i,2),'color','k','lineStyle','-')
    %     scatter(subplot(2,2,3),i*ones(1,size(aP_Models{1,2},2)),aP_Models{2}(i,:),'k','filled','jitter', 'on', 'jitterAmount', jitterAmount2);
end
set(subplot(2,2,3),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

%%%%%%%%%%%%%%%%%%% (C) BarPlots of Semi-saturation constants for all Neural Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4)
b3 = bar(sP_Bar(:,3)); % population model for sP
if strcmp(colorScheme,'color')
    b3.FaceColor = 'y';
elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
    b3.FaceColor = [0.33 0.33 0.33];
end
hold on
for i=1:4 % Neural measures
    errorbar(subplot(2,2,4),i,sP_Bar(i,3),sP_ErrorBar(i,3),'color','k','lineStyle','-')
    %     scatter(subplot(2,2,4),i*ones(1,size(sP_Models{1,3},2)),sP_Models{3}(i,:),'k','filled','jitter', 'on', 'jitterAmount', jitterAmount2);
end

set(subplot(2,2,4),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

title(subplot(2,2,1),'% Explained Variance' ,'fontSize',14);
title(subplot(2,2,2),'AIC Values' ,'fontSize',14);
title(subplot(2,2,3),'Tuned Normalization Parameter (\alpha)','fontSize',14);
title(subplot(2,2,4),'Semi-saturation Constant (\sigma)','fontSize',14);

textH{1} = getPlotHandles(1,1,[0.07 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.51 0.95 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.07 0.5 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.51 0.5 0.01 0.01]);

textString = {'A','B','C','D'};

for i=1:4 % Plot numbers
    ylabel(subplot(2,2,i),'Population Median')
    set(subplot(2,2,i),'xticklabel',neuralMeasures,'XTickLabelRotation',45,'fontSize',14)
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end



end

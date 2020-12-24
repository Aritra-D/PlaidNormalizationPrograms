function hFigure = plotStatisticsForSeparateMonkeys(folderSourceString,stimType,elecParams,timeRangeForComputation,tapers_MT,combineUniqueElectrodeData,colorScheme)

% close all; % uncomment if this function is run separately!You will also
% need elecParams from plotFigures
if ~exist(folderSourceString,'var'); folderSourceString = 'E:\'; end
monkeyName{1} = 'alpaH'; monkeyName{2} = 'kesariH';


if strcmp(stimType,'Static')
    folderData = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures');
elseif strcmp(stimType,'Counterphase')
    folderData = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures\Counterphase');
end

gridType = 'microelectrode';
removeERPFlag = 0;

dataParams = ['_N' num2str(elecParams.spikeCutoff) '_S' num2str(elecParams.snrCutoff) '_allElecs'...
    '_T' num2str(round(1000*timeRangeForComputation(1))) '_' num2str(round(1000*timeRangeForComputation(2))) ...
    '_d' num2str(elecParams.dRange(1)) '_' num2str(elecParams.dRange(2))...
    '_tapers' num2str(tapers_MT(2)) '_removeERP' num2str(removeERPFlag) '_cne' num2str(0) ...
    '_gse' num2str(elecParams.getSpikeElectrodesFlag) '_' gridType '_UnitID' num2str(elecParams.unitID) '.mat'];

N_Spike=elecParams.spikeCutoff; snr=elecParams.snrCutoff; d=elecParams.dRange(2);



models = {'Untuned','Tuned-alpha_c2','Population'};
neuralMeasures = {'Spikes','Gamma','Hi-Gamma','SSVEP'};
params = {'N.I.','aP','sP','eP','AIC'};
modelPairs = nchoosek(1:length(models),2);
neuralMeasurePairs = nchoosek(1:length(neuralMeasures),2);

for iMonkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Aritra - if you make one file for N=0, S=0 (all electrodes) you can
    % always sub-select based on N,snr and d.
    
    fileSave = fullfile(folderData,[monkeyName{iMonkey} dataParams]);
    load(fileSave) %#ok<*LOAD>
    
    
    %%%%%%%%%%%%%%%% Subselect electrodes based on N, snr and d %%%%%%%%%%%%%%%
    goodN = (max([elecInfo.N(1,:);elecInfo.N(2,:)])>N_Spike);
    goodSNR = (elecInfo.SNR>snr);
    goodD = (elecInfo.d <= d);
    goodElectrodeIDs = (goodN & goodSNR & goodD);
    
    %%%%%%%%%%%%%%%%%%% Get information about the electrodes %%%%%%%%%%%%%%%%%%
    allElectrodes =[];
    for i=1:length(elecInfo.elecs)
        x = elecInfo.elecs{i}{end};
        allElectrodes = cat(2,allElectrodes,x);
    end
    
    %%%%%% For each electrode, find the IDs that need to be averaged %%%%%%%%%%
    clear goodIDList
    if combineUniqueElectrodeData
        
        uniqueGoodElectrodeList = unique(allElectrodes(goodElectrodeIDs));
        numGoodElectrodes = length(uniqueGoodElectrodeList);
        
        goodIDList{iMonkey} = cell(1,numGoodElectrodes);
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
            %         y = x; % use for symmetric model
            %         z = flip(y,1); z = (z+z')/2; z = flip(z,1);
            allDataTMP(i,:,:) = x;
        else
            xs = zeros(size(x));
            for j=1:length(goodIDList{i})
                %             y = squeeze(x(j,:,:)); % use for symmetric model
                %             z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
                xs(j,:,:) = squeeze(x(j,:,:));
            end
            allDataTMP(i,:,:) = squeeze(mean(xs,1));
        end
    end
    allData{1} = allDataTMP;
    allDataType{1} = 'FR'; allDataType{2} = 'G'; allDataType{3} = 'HG'; allDataType{4} = 'S';
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1- Standard Divisive Normalization (Untuned normalization)
    % 2- Tuned Normalization
    % 3- EMS-stimulus Normalization on actual dataset
    % 4- EMS-stimulus Normalization on symmetrized dataset (not used in this analysis)
    
    for iModel =1:3
        modelNum = iModel;
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Get Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        eP = zeros(4,numGoodElectrodes);
        aP = zeros(4,numGoodElectrodes);
        sP = zeros(4,numGoodElectrodes);
        dataP = zeros(4,numGoodElectrodes,5,5);
        rSS = zeros(4,numGoodElectrodes);
        aicVals = zeros(4,numGoodElectrodes);
        
        for i=1:4
            for j=1:numGoodElectrodes
                
                %Plaid
                z = squeeze(allData{i}(j,:,:));
                %         z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
                parP = getParametersPlaidV2(z,modelNum);
                [dp,pz] = getResponseMatrixPlaidV2(parP,z,modelNum); %#ok<*ASGLU>
                eP(i,j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
                rSS(i,j) = dp; % storing dps as residual squared sum for estimation of maximum likelihood for each model
                aicVals(i,j) = calculateAIC(dp,numel(z),size(parP,2));
                
                if modelNum==1
                    %             aP(i,j) = 0;
                    sP(i,j) = min(max(0,parP(3)),5);
                else
                    aP(i,j) = parP(3);
                    sP(i,j) = min(max(0,parP(4)),5);
                end
                dataP(i,j,:,:) = z;
                
            end
        end
        
        if modelNum==2 || modelNum==3
            aP_Models{iMonkey}{modelNum-1} = aP;
        end
        sP_Models{iMonkey}{iModel} = sP;
        eP_Models{iMonkey}{iModel} = eP;
        rSS_Models{iMonkey}{iModel} = rSS; %#ok<*NASGU>
        nObs_Models{iMonkey}{iModel} = numel(z);
        params_Models{iMonkey}{iModel} = size(parP,2);
        aic_Models{iMonkey}{iModel} = aicVals;
        
        
        
        
        NI = squeeze((dataP(:,:,1,1)+dataP(:,:,5,5))./dataP(:,:,1,5));
        NI_monkeys{iMonkey} = NI;
        % Statistics
        disp(['Monkey: ' monkeyName{iMonkey}])
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
                        disp([neuralMeasures{j},':' num2str(round(median(NI(j,:)),3)),' +/- ' num2str(round(getSEMedian(NI(j,:)',1000),3))])
                    end
                    
                elseif i==2
                    if modelNum==2 || modelNum==3
                        aP_median_temp = median(aP(j,:));
                        aP_seMedian_temp = getSEMedian(aP(j,:)',1000);
                        disp([neuralMeasures{j},':' num2str(round(aP_median_temp,2)),' +/- ' num2str(round(aP_seMedian_temp,2))])
                        aP_median{iMonkey}{modelNum-1}(j) = aP_median_temp;
                        aP_seMedian{iMonkey}{modelNum-1}(j)= aP_seMedian_temp;
                    end
                    
                elseif i==3
                    sP_median_temp = median(sP(j,:));
                    sP_seMedian_temp = getSEMedian(sP(j,:)',1000);
                    disp([neuralMeasures{j},':' num2str(round(sP_median_temp,2)),' +/- ' num2str(round(sP_seMedian_temp,2))])
                    sP_median{iMonkey}{iModel}(j) = sP_median_temp;
                    sP_seMedian{iMonkey}{iModel}(j)= sP_seMedian_temp;
                    
                elseif i==4
                    eP_median_temp = median(eP(j,:));
                    eP_seMedian_temp = getSEMedian(eP(j,:)',1000);
                    disp([neuralMeasures{j},':' num2str(round(eP_median_temp,2)),' +/- ' num2str(round(eP_seMedian_temp,3))])
                    eP_median{iMonkey}{iModel}(j) = eP_median_temp;
                    eP_seMedian{iMonkey}{iModel}(j)= eP_seMedian_temp;
                    
                elseif i==5
                    aic_median_temp = median(aicVals(j,:));
                    aic_seMedian_temp = getSEMedian(aicVals(j,:)',1000);
                    disp([neuralMeasures{j},':' num2str(round(aic_median_temp,2)),' +/- ' num2str(round(aic_seMedian_temp,3))])
                    aic_median{iMonkey}{iModel}(j) = aic_median_temp;
                    aic_seMedian{iMonkey}{iModel}(j)= aic_seMedian_temp;
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
end

% Significance Test across normalization models for each neural measure

for iMonkey =1:2
    disp(['Monkey: ' monkeyName{iMonkey}])
    disp('Statistical tests for AIC- Kruskal-Wallis test')
    for i=1:length(neuralMeasures)
        for j= 1:size(modelPairs,1)
            clear statData4
            statData4 = [aic_Models{iMonkey}{modelPairs(j,1)}(i,:)',aic_Models{iMonkey}{modelPairs(j,2)}(i,:)'];
            p_AIC{iMonkey}(i,j) = kruskalwallis(statData4,[],'off');
            disp([neuralMeasures{i} ' : p between ' models{modelPairs(j,1)} ' and ' models{modelPairs(j,2)} ' : ' num2str(p_AIC{iMonkey}(i,j))])
        end
    end
    
    
    disp('Statistical tests for Wilcoxon Signed Rank test')
    for i=1:length(neuralMeasures)
        for j= 1:size(modelPairs,1)
            [p_AIC_Wil{iMonkey}(i,j),h_AIC_Wil{iMonkey}(i,j),stats_AIC_wil{iMonkey}{i,j}] = signrank(aic_Models{iMonkey}{modelPairs(j,2)}(i,:),aic_Models{iMonkey}{modelPairs(j,1)}(i,:));
            disp([neuralMeasures{i} ' : p between ' models{modelPairs(j,2)} ' and ' models{modelPairs(j,1)} ' : ' num2str(p_AIC_Wil{iMonkey}(i,j))])
        end
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

% % calculate aic for different models and Neural measures
% for iMonkey = 1:2
%     for iModel =1:3
%         aic_Models{iMonkey}{iModel} = nObs_Models{iMonkey}{iModel} .* log(rSS_Models{iMonkey}{iModel}./nObs_Models{iMonkey}{iModel})+2.*(params_Models{iMonkey}{iModel}+1);
%     end
% end
%
% % statistical distribution and Krushkal-Wallis test of AIC values of normalization models
% disp('AIC for Normalization models:')
% for iMonkey = 1:2
%     disp([monkeyName{iMonkey} ':' ])
%     for i = 1:4 % Neural measures
%         for iModel =1:3
%             aic_median = median(aic_Models{iMonkey}{iModel}(i,:));
%             aic_seMedian = getSEMedian(aic_Models{iMonkey}{iModel}(i,:)',1000);
%             disp([neuralMeasures{i}, ' -- Model: ' iModel ': ' modelStringList{iModel} ': median +/- SEMedian: ' num2str(round(aic_median,2)),' +/- ' num2str(round(aic_seMedian,2))])
%
%
%             if iModel == 3
%                 models = 1:1:length(modelStringList);
%                 aic_grouped{iMonkey,i}= [aic_Models{iMonkey}{models(1)}(i,:)', aic_Models{iMonkey}{models(2)}(i,:)', aic_Models{iMonkey}{models(3)}(i,:)'];
%                 p_aic(iMonkey,i) = kruskalwallis(aic_grouped{iMonkey,i},[],'off');
%                 disp(['Krushkal Wallis test, data compared to all models: p =  ' num2str(p_aic(iMonkey,i))])
%
%                 modelPairs = nchoosek(1:length(modelStringList),2);
%                 for j=1:size(modelPairs,1)
%                     aic_groupedModelPairs{iMonkey}{i}{j} = [aic_Models{iMonkey}{modelPairs(j,1)}(i,:)', aic_Models{iMonkey}{modelPairs(j,2)}(i,:)'];
%                     p_aic_ModelPairs(iMonkey,i,j) = kruskalwallis(aic_groupedModelPairs{iMonkey}{i}{j},[],'off');
%                     disp(['Krushkal Wallis test, data compared between Model ' num2str(modelPairs(j,1)) ' and Model ' num2str(modelPairs(j,2)) ': p =  ' num2str(p_aic_ModelPairs(iMonkey,i,j))])
%
%                 end
%             end
%         end
%     end
% end



% grouping parameter values according to neural measures
for iMonkey = 1:2
    for i=1:4 % Neural Measures
        aPBar_median{iMonkey}{i} = [aP_median{iMonkey}{1}(i) aP_median{iMonkey}{2}(i)]; %#ok<*AGROW>
        aPBar_seMedian{iMonkey}{i} = [aP_seMedian{iMonkey}{1}(i) aP_seMedian{iMonkey}{2}(i)];
        sPBar_median{iMonkey}{i} = [sP_median{iMonkey}{1}(i) sP_median{iMonkey}{2}(i) sP_median{iMonkey}{3}(i)];
        sPBar_seMedian{iMonkey}{i} = [sP_seMedian{iMonkey}{1}(i) sP_seMedian{iMonkey}{2}(i) sP_seMedian{iMonkey}{3}(i)];
        ePBar_median{iMonkey}{i} = [eP_median{iMonkey}{1}(i) eP_median{iMonkey}{2}(i) eP_median{iMonkey}{3}(i)];
        ePBar_seMedian{iMonkey}{i} = [eP_seMedian{iMonkey}{1}(i) eP_seMedian{iMonkey}{2}(i) eP_seMedian{iMonkey}{3}(i)];
        aicBar_median{iMonkey}{i} = [aic_median{iMonkey}{1}(i) aic_median{iMonkey}{2}(i) aic_median{iMonkey}{3}(i)];
        aicBar_seMedian{iMonkey}{i} = [aic_seMedian{iMonkey}{1}(i) aic_seMedian{iMonkey}{2}(i) aic_seMedian{iMonkey}{3}(i)];
        
    end
end

% converting into arrays from cells that are used to plot barplots



for iMonkey = 1:2 % Monkeys
    aP_Bar{iMonkey} =[];aP_ErrorBar{iMonkey}=[];
    sP_Bar{iMonkey} =[];sP_ErrorBar{iMonkey}=[];
    eP_Bar{iMonkey} =[];eP_ErrorBar{iMonkey}=[];
    aic_Bar{iMonkey} = []; aic_ErrorBar{iMonkey} = [];
    
    for i=1:4 % Neural Measures
        aP_Bar{iMonkey} = cat(1,aP_Bar{iMonkey},aPBar_median{iMonkey}{i});
        aP_ErrorBar{iMonkey} = cat(1,aP_ErrorBar{iMonkey},aPBar_seMedian{iMonkey}{i});
        
        sP_Bar{iMonkey} = cat(1,sP_Bar{iMonkey},sPBar_median{iMonkey}{i});
        sP_ErrorBar{iMonkey} = cat(1,sP_ErrorBar{iMonkey},sPBar_seMedian{iMonkey}{i});
        
        eP_Bar{iMonkey} = cat(1,eP_Bar{iMonkey},ePBar_median{iMonkey}{i});
        eP_ErrorBar{iMonkey} = cat(1,eP_ErrorBar{iMonkey},ePBar_seMedian{iMonkey}{i});
        
        aic_Bar{iMonkey} = cat(1,aic_Bar{iMonkey},aicBar_median{iMonkey}{i});
        aic_ErrorBar{iMonkey} = cat(1,aic_ErrorBar{iMonkey},aicBar_seMedian{iMonkey}{i});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Figure 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFigure = figure(6);
set(hFigure,'units','normalized','outerposition',[0 0 1 1])

jitterAmount = 0.03;

for iMonkey =1:2
    
    %%%%%%%%%%%%%%%%%%% (A) BarPlots of Explained Variance, Monkey 1 %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%  BarPlots of Explained Variance, Monkey 2 %%%%%%%%%%%%%%%%%%%%
    subplot(3,2,iMonkey)
    b1 = bar(100*eP_Bar{iMonkey}, 'grouped');
    
    if strcmp(colorScheme,'color')
        b1(3).FaceColor = 'y';
    elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
        colorVals = flip(gray(4)); colorVals(4,:) = [];
        colors = colorVals;
        for i = 1:size(eP_Bar{iMonkey},2)
            b1(i).FaceColor = colors(i,:,:);
        end
    end
    
    hold on
    
    % Find the number of groups and the number of bars in each group
    ngroups = size(eP_Bar{iMonkey}, 1);
    nbars = size(eP_Bar{iMonkey}, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars % number of bars = number of models
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        
        if iMonkey ==2
            if i==3
                err_Bar = errorbar(x(1), 100*eP_Bar{iMonkey}(1,i), 100*eP_ErrorBar{iMonkey}(1,i),'k', 'linestyle', 'none','lineWidth',1);
                err_Bar.YPositiveDelta=[];
                errorbar(x(2:4), 100*eP_Bar{iMonkey}(2:4,i), 100*eP_ErrorBar{iMonkey}(2:4,i),'k', 'linestyle', 'none','lineWidth',1);
            else
                errorbar(x, 100*eP_Bar{iMonkey}(:,i), 100*eP_ErrorBar{iMonkey}(:,i),'k', 'linestyle', 'none','lineWidth',1);
            end
        else
            errorbar(x, 100*eP_Bar{iMonkey}(:,i), 100*eP_ErrorBar{iMonkey}(:,i),'k', 'linestyle', 'none','lineWidth',1);
        end
        
        
        for j = 1:ngroups
            scatter(x(j)*ones(1,size(eP_Models{iMonkey}{i},2)),100*eP_Models{iMonkey}{i}(j,:),'k','filled','jitter', 'on', 'jitterAmount', jitterAmount);
        end
    end
    err_Bar.YPositiveDelta = [];
    tickLengthPlot = 3*get(subplot(3,2,iMonkey),'TickLength');
    set(subplot(3,2,iMonkey),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
    ylim(subplot(3,2,iMonkey),[0 100]);
    legend(subplot(3,2,1),{'Standard','Tuned','Population'},'fontSize',12,'Location','best')
    
    %%%%%%%%%%%%%%%%%%% (B) BarPlots of AIC values, Monkey 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%  BarPlots of AIC values, Monkey 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,2,2+iMonkey)
    b2 = bar(aic_Bar{iMonkey}, 'grouped');
    
    if strcmp(colorScheme,'color')
        b2.FaceColor = 'y';
    elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
        colorVals = flip(gray(4)); colorVals(4,:) = [];
        colors = colorVals;
        for i = 1:size(aic_Bar{iMonkey},2)
            b2(i).FaceColor = colors(i,:,:);
        end
    end
    hold on
    % Find the number of groups and the number of bars in each group
    ngroups = size(aic_Bar{iMonkey}, 1);
    nbars = size(aic_Bar{iMonkey}, 2);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, aic_Bar{iMonkey}(:,i), aic_ErrorBar{iMonkey}(:,i), 'k', 'linestyle', 'none','lineWidth',1);
        for j = 1:ngroups
            scatter(x(j)*ones(1,size(aic_Models{iMonkey}{i},2)),aic_Models{iMonkey}{i}(j,:),'k','filled','jitter','on','jitterAmount',jitterAmount)
        end
    end
    
    set(subplot(3,2,2+iMonkey),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)
    % ylim(subplot(2,2,2),[0 100]);
    if iMonkey == 1
        legend(subplot(3,2,2+iMonkey),{'Standard','Tuned','Population'},'fontSize',12,'Location','northeast')
    end
end

%%%%%%%%%%%%%%%%%%% (C) BarPlots of Tuned Normalization parameters for all Neural Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,2,5)
modelNum = 2; % population model for aP % first model does not have aP
aP_Bar_BothMonkeys = [aP_Bar{1}(:,modelNum), aP_Bar{2}(:,modelNum)];
aP_ErrorBar_BothMonkeys = [aP_ErrorBar{1}(:,modelNum), aP_ErrorBar{2}(:,modelNum)];

b3 = bar(aP_Bar_BothMonkeys,'grouped'); % population model (aP is not present in untuned normalization model)
% if strcmp(colorScheme,'color')
% %     b2(3).FaceColor = 'y';
%
% elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
% %     colors = repmat(0.65:-0.1:0.45,[3 1])';
%     for i = 1:size(aP_Bar_BothMonkeys,2)
%         b2(i).FaceColor = {'k','w'};
%     end
% end
colors = [0.45 0.45 0.45; 0.75 0.75 0.75];
for i = 1:size(aP_Bar_BothMonkeys,2)
    b3(i).FaceColor = colors(i,:,:);
end
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(aP_Bar_BothMonkeys, 1);
nbars = size(aP_Bar_BothMonkeys, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars % nbars = nMonkeys
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,aP_Bar_BothMonkeys(:,i),aP_ErrorBar_BothMonkeys(:,i),'color','k','lineStyle','none')
    for j = 1:ngroups
        scatter(x(j)*ones(1,size(aP_Models{i}{modelNum},2)),aP_Models{i}{modelNum}(j,:),'k','filled','jitter', 'on', 'jitterAmount', jitterAmount);
    end
end
legend(subplot(3,2,5),{'Monkey 1','Monkey 2'},'fontSize',12,'Location','northeast')
set(subplot(3,2,5),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

%%%%%%%%%%%%%%%%%%% (D) BarPlots of Semi-saturation constants for all Neural Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,2,6)
modelNum =3; % population model for sP
sP_Bar_BothMonkeys = [sP_Bar{1}(:,modelNum), sP_Bar{2}(:,modelNum)];
sP_ErrorBar_BothMonkeys = [sP_ErrorBar{1}(:,modelNum), sP_ErrorBar{2}(:,modelNum)];

b4 = bar(sP_Bar_BothMonkeys,'grouped'); % population model (aP is not present in untuned normalization model)
% if strcmp(colorScheme,'color')
% %     b2(3).FaceColor = 'y';
%
% elseif strcmp(colorScheme,'grayscale')||strcmp(colorScheme,'greyscale')
% %     colors = repmat(0.65:-0.1:0.45,[3 1])';
%     for i = 1:size(aP_Bar_BothMonkeys,2)
%         b2(i).FaceColor = {'k','w'};
%     end
% end
colors = [0.45 0.45 0.45; 0.75 0.75 0.75];
for i = 1:size(aP_Bar_BothMonkeys,2)
    b4(i).FaceColor = colors(i,:,:);
end
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(sP_Bar_BothMonkeys, 1);
nbars = size(sP_Bar_BothMonkeys, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars % nbars = nMonkeys
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,sP_Bar_BothMonkeys(:,i),sP_ErrorBar_BothMonkeys(:,i),'color','k','lineStyle','none')
    for j = 1:ngroups
        scatter(x(j)*ones(1,size(sP_Models{i}{modelNum},2)),sP_Models{i}{modelNum}(j,:),'k','filled','jitter', 'on', 'jitterAmount', jitterAmount);
    end
end
% legend(subplot(1,4,4),{'Monkey 1','Monkey 2'},'fontSize',14,'Location','best')

set(subplot(3,2,6),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',14)

% Set Axes Properties
title(subplot(3,2,1),'% Explained Variance, Monkey 1','fontSize',14);
title(subplot(3,2,2),'% Explained Variance, Monkey 2','fontSize',14);
title(subplot(3,2,3),'AIC Values, Monkey 1' ,'fontSize',14);
title(subplot(3,2,4),'AIC Values, Monkey 2' ,'fontSize',14);
title(subplot(3,2,5),'Tuned Normalization Parameter (\alpha)','fontSize',14);
title(subplot(3,2,6),'Semi-saturation Constant (\sigma)','fontSize',14);


textH{1} = getPlotHandles(1,1,[0.07 0.95 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.07 0.64 0.01 0.01]);
textH{3} = getPlotHandles(1,1,[0.07 0.34 0.01 0.01]);
textH{4} = getPlotHandles(1,1,[0.51 0.35 0.01 0.01]);

textString = {'A','B','C','D'};

for i=1:4 % Plot numbers
    ylabel(subplot(3,2,i),'Population Median')
    set(subplot(3,2,i),'xticklabel',[],'fontSize',14)
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
ylabel(subplot(3,2,5),'Population Median')
set(subplot(3,2,5),'xticklabel',neuralMeasures,'XTickLabelRotation',45,'fontSize',14)
ylabel(subplot(3,2,6),'Population Median')
set(subplot(3,2,6),'xticklabel',neuralMeasures,'XTickLabelRotation',45,'fontSize',14)
if combineUniqueElectrodeData == 0
    ylim(subplot(3,2,5),[-10 50]);    ylim(subplot(3,2,6),[0 10]);
end
end

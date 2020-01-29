clear;
combineUniqueElectrodes=1;
N=15; snr=2; d=0.75;
versionNum = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aritra - if you make one file for N=0, S=0 (all electrodes) you can
% always sub-select based on N,snr and d.
load('E:\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\all_N15_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat'); 
N=15; snr=2; d=0.75;
% firingRateData = firingRateData2;
% load('E:\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\savedData_Figures\elecInfo_all_N15_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');

% load('savedData_Figures\all_N15_S2_allElecs_T50_250_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');
% load('savedData_Figures\elecInfo_all_N15_S2_allElecs_T50_250_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');

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
        parP = getParametersPlaid(z,versionNum);
        [dp,pz] = getResponseMatrixPlaid(parP,z,versionNum);
        eP(i,j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
        aP(i,j) = parP(2);
        sP(i,j) = parP(3);
        dataP(i,j,:,:) = z;
    end
end

cList = [0 1 2 4 8]/16;


meanP = squeeze(mean(dataP,2));
figure(6);
for i=1:4
    subplot(3,4,i);
    imagesc(squeeze(meanP(i,:,:))); colorbar;
    
    mp = squeeze(meanP(i,:,:));
    parMP = getParametersPlaid(mp,versionNum);
    [dp,pmp] = getResponseMatrixPlaid(parMP,mp,versionNum);
    expVar = 1 - (dp/sum((mp(:)-mean(mp(:))).^2));
    
    subplot(3,4,4+i);
    plot(cList,mp(5,:),'color',[0.5 0.5 0.5],'marker','o','linestyle','none');
    hold on;
    plot(cList,[mp(5,1) mp(4,2) mp(3,3) mp(2,4) mp(1,5)],'color','k','marker','v','linestyle','none');
    plot(cList,pmp(5,:),'color',[0.5 0.5 0.5]);
    plot(cList,[pmp(5,1) pmp(4,2) pmp(3,3) pmp(2,4) pmp(1,5)],'color','k');
    title(['\alpha=' num2str(parMP(2),3) ', \sigma=' num2str(parMP(3),3) ',ExpVar=' num2str(round(100*expVar)) '%']);
end

subplot(3,4,9)
NI = squeeze(2*dataP(:,:,1,1)./dataP(:,:,1,5));
bar(median(NI,2));
hold on;
errorbar(1:4,median(NI,2),getSEMedian(NI'),'k.');

subplot(3,4,10)
bar(median(eP,2));
hold on;
errorbar(1:4,median(eP,2),getSEMedian(eP'),'k.');

subplot(3,4,11)
bar(median(aP,2));
hold on;
errorbar(1:4,median(aP,2),getSEMedian(aP'),'k.');

subplot(3,4,12)
bar(median(sP,2));
hold on;
errorbar(1:4,median(sP,2),getSEMedian(sP'),'k.');
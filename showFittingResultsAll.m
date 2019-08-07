%load('C:\Users\Supratim Ray\OneDrive - Indian Institute of Science\Supratim\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\all_N20_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_gridType_Microelectrode_UnitID0.mat');

if strcmp(getenv('username'),'Aritra') || strcmp(getenv('username'),'Lab Computer-Aritra')
    folderSourceString_Project = 'E:\';
elseif strcmp(getenv('username'),'Supratim Ray')
    folderSourceString_Project = 'M:\';
end

folderName = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures');
fileName = fullfile(folderName,'all_N20_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');

if exist(fileName,'file')
    load(fullfile(fileName));
    disp(['Loading file ' fileName]);
end

figure;
xAll{1} = squeeze(firingRateData.analysisDataST(:,1,:,:)) - squeeze(firingRateData.analysisData_cBL(:,1,:,:));
type{1} = 'FR';

for i=2:4
    xAll{i} = squeeze(energyData.analysisDataST{i}) - squeeze(energyData.analysisData_cBL{i});
end

type{2} = 'G'; type{3} = 'HG'; type{4} = 'S';

numElectrodes = length(xAll{1});

for i=1:4
    for j=1:numElectrodes
        
        % Grating
        g1 = squeeze(xAll{i}(j,5,:))';
        g2 = flip(squeeze(xAll{i}(j,:,1)),2);
        
        g = (g1+g2)/2; % Make symmetric
        parG = getParametersGrating(g);
        [dg,pg] = getResponseMatrixGrating(parG,g);
        eG(i,j) = 1 - (dg/sum((g-mean(g)).^2));
        sG(i,j) = parG(2);
        dataG(i,j,:) = g;
        
        %Plaid
        y = squeeze(xAll{i}(j,:,:));
        z = flip(y,1); z = (z+z')/2; z = flip(z,1); % Make symmetric
        parP = getParametersPlaid(z);
        [dp,pz] = getResponseMatrixPlaid(parP,z);
        eP(i,j) = 1 - (dp/sum((z(:)-mean(z(:))).^2));
        aP(i,j) = parP(2);
        sP(i,j) = parP(3);
        dataP(i,j,:,:) = z;
    end
end

cList = [0 1 2 4 8]/16;

% subplot(131);
% m = squeeze(median(dataG,2));
% colorNames = jet(4);
% for i=1:4
%     nm = m(i,:)/m(i,5);
%     plot(cList,nm,'marker','o','color',colorNames(i,:),'linestyle','none');
%     hold on;
%     
%     pnm(i,:) = getParametersGrating(nm);
% end
% 
% for i=1:4
%     [~,pg] = getResponseMatrixGrating(pnm(i,:),[]);
%     plot(cList,pg,'color',colorNames(i,:));
%     hold on;
%     legendStr{i} = [type{i} ',s=' num2str(pnm(i,2),3)];
% end
% 
% legend(legendStr,'location','best');
% 
% subplot(132);
% barm = median(eG,2);
% bar(barm);
% 
% subplot(133);
% barm = median(sG,2);
% bar(barm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

meanP = squeeze(mean(dataP,2));
for i=1:4
    subplot(3,4,i);
    imagesc(squeeze(meanP(i,:,:))); colorbar;
    
    mp = squeeze(meanP(i,:,:));
    parMP = getParametersPlaid(mp);
    [dp,pmp] = getResponseMatrixPlaid(parMP,mp);
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

subplot(3,4,10)
bar(median(eP,2));

subplot(3,4,11)
bar(median(aP,2));

subplot(3,4,12)
bar(median(sP,2));
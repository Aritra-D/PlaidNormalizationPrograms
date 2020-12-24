function displayNIVersusPrefOri

close all;
load('E:\Projects\Aritra_PlaidNormalizationProject\savedData_Figures\all_N15_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat'); %#ok<*LOAD>

electrodeList{1} = 1:143; % Non-Unique Elecs- Monkey 1
electrodeList{2} = 144:191; % Non-Unique Elecs- Monkey 2

hFig = figure(1);
set(hFig,'units','normalized','outerposition',[0 0 1 1])
xTicks_PO = [0 90 180]; yTicks_PO = xTicks_PO; 
xLims_PO = [-20 200]; yLims_PO = xLims_PO;
fontSize = 14; 

xTicks_OS = [0 0.25]; yTicks_OS = 0:0.25:0.75; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
xLims_OS = [-0.05 0.3]; yLims_OS = [-0.1 0.85];

x = oriTuningData.PO{1};
y = oriTuningData.oValsUnique_Plaid;
% z = oriTuningData.oValsUnique2_Plaid;

% for i =1:4
% oriDiffTMP(i,:) = abs(x-(y+(i-1)*90));
% 
% end
oriDiffTMP = [abs(x-y); abs(x-(y+90)); abs(y-x+180); abs(x-y+90)]; 
oriDiff = min(oriDiffTMP,[],1);


for iMonkey = 1:2
clear idx rho pVals rho2 pVals2 rho3 pVals3
idx = 3*iMonkey-2;
subplot(2,3,idx); 
tickLengthPlot = 4*get(subplot(2,3,idx),'TickLength');
xlim(subplot(2,3,idx),xLims_PO); ylim(subplot(2,3,idx),yLims_PO);
% title (subplot(2,3,iMonkey),['Monkey ' num2str(iMonkey)]);
scatter(oriTuningData.PO{2}(electrodeList{iMonkey}),oriTuningData.PO{1}(electrodeList{iMonkey}),'k','filled'); axis square
[rho,pVals] = corr(oriTuningData.PO{2}(electrodeList{iMonkey})',oriTuningData.PO{1}(electrodeList{iMonkey})','Type','Spearman');
title(subplot(2,3,idx),['r = ' num2str(round(rho,2)) ', p = ' num2str(pVals)])
xlabel('Pref. Ori. (LFP Gamma)'); ylabel('Pref. Ori. (spikes)');
set(subplot(2,3,idx),'xTick',xTicks_PO,'yTick',yTicks_PO,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);

subplot(2,3,idx+1);
xlim(subplot(2,3,idx+1),xLims_OS); ylim(subplot(2,3,idx+1),yLims_OS);
scatter(oriTuningData.OS{2}(electrodeList{iMonkey}),oriTuningData.OS{1}(electrodeList{iMonkey}),'k','filled'); axis square
[rho2,pVals2] = corr(oriTuningData.OS{2}(electrodeList{iMonkey})',oriTuningData.OS{1}(electrodeList{iMonkey})','Type','Spearman');
title(subplot(2,3,idx+1),['r = ' num2str(round(rho2,2)) ', p = ' num2str(pVals2)])
xlabel('OS (LFP Gamma)'); ylabel('OS (spikes)');% title ('Monkey 1');
set(subplot(2,3,idx+1),'xTick',xTicks_OS,'yTick',yTicks_OS,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);

subplot(2,3,idx+2);hold on;
xlim(subplot(2,3,idx+2),[-5 50]); ylim(subplot(2,3,idx+2),[-2 10]);
scatter(oriDiff(electrodeList{iMonkey}),NI_Data.dfiringRate_Cohen(electrodeList{iMonkey}),'k','filled'); axis square
[rho3,pVals3] = corr(oriDiff(electrodeList{iMonkey})',NI_Data.dfiringRate_Cohen(electrodeList{iMonkey}),'Type','Spearman');
title(subplot(2,3,idx+2),['r = ' num2str(round(rho3,2)) ', p = ' num2str(pVals3)])
xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Spikes)');% title ('Monkey 1');

set(subplot(2,3,idx+2),'xTick',[0 22.5 45],'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);

textString = ['Monkey ' num2str(iMonkey)];
textH{iMonkey} = getPlotHandles(1,1,[0.05 0.7-0.48*(iMonkey-1) 0.01 0.01]); %#ok<*AGROW>
set(textH{iMonkey},'Visible','Off');
text(0.1,0.5,textString,'unit','normalized','fontsize',fontSize+4,'fontweight','bold','rotation',90,'parent',textH{iMonkey});


end
%%
% hFig2 = figure(2);

% Spikes
% subplot(3,2,1)
% 
% y = oriTuningData.PO{1};
% y1 = oriTuningData.oValsUnique_Plaid;
% y2 = oriTuningData.oValsUnique2_Plaid;
% 
% yminusy1 = y-y1;
% % yminusy1comp = y-(180-y1);
% yminusy2 = y-y2;
% % yminusy2comp = y-(180-y2);
% alldiff = [yminusy1; yminusy2]; % yminusy1comp; yminusy2comp];
% 
% for i=1:191
%     clear idx
%     idx = find(abs(alldiff(:,i))==min(abs(alldiff(:,i))),1,'first');
%     OriDiff(i) = alldiff(idx,i); %#ok<*SAGROW>
%     ids(i) = idx;
%     if idx == 1
%         PrefOriCloserToOri1orOri2 (i)= 1; 
%     elseif idx == 2
%         PrefOriCloserToOri1orOri2 (i)= 2; 
%     end
% end
% 
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==1),NI_Data.dfiringRate_Cohen(PrefOriCloserToOri1orOri2(1:143)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==2),NI_Data.dfiringRate_Cohen(PrefOriCloserToOri1orOri2(1:143)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Spikes)');% title ('Monkey 1');
% % xTicks = [-45 0 45]; yTicks = [1 3]; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
% xLims = [-180 180]; 
% yLims = [-1 5];
% xlim(xLims); 
% ylim(yLims);
% title(subplot(3,2,1),'Monkey 1')
% set(subplot(3,2,1),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% subplot(3,2,2)
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==1),NI_Data.dfiringRate_Cohen(PrefOriCloserToOri1orOri2(144:end)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==2),NI_Data.dfiringRate_Cohen(PrefOriCloserToOri1orOri2(144:end)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Spikes)');% title ('Monkey 1');
% title(subplot(3,2,2),'Monkey 2')
% xlim(xLims); 
% ylim(yLims);
% set(subplot(3,2,2),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% 
% % Gamma
% subplot(3,2,3)
% 
% y = oriTuningData.PO{2};
% y1 = oriTuningData.oValsUnique_Plaid;
% y2 = oriTuningData.oValsUnique2_Plaid;
% 
% yminusy1 = y-y1;
% % yminusy1comp = y-(180-y1);
% yminusy2 = y-y2;
% % yminusy2comp = y-(180-y2);
% alldiff = [yminusy1; yminusy2];% yminusy1comp; yminusy2comp];
% 
% for i=1:191
%     clear idx
%     idx = find(abs(alldiff(:,i))==min(abs(alldiff(:,i))),1,'first');
%     OriDiff(i) = alldiff(idx,i);
%     ids(i) = idx;
%     if idx == 1
%         PrefOriCloserToOri1orOri2 (i)= 1; 
%     elseif idx == 2
%         PrefOriCloserToOri1orOri2 (i)= 2; 
%     end
% end
% 
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==1),NI_Data.denergy_Cohen{2}(PrefOriCloserToOri1orOri2(1:143)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==2),NI_Data.denergy_Cohen{2}(PrefOriCloserToOri1orOri2(1:143)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Gamma)');% title ('Monkey 1');
% % xTicks = [-45 0 45]; yTicks = [1 3]; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
% % xLims = [-50 50]; yLims = [-1 10];
% xlim(xLims); 
% % ylim(yLims);
% set(subplot(3,2,3),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% subplot(3,2,4)
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==1),NI_Data.denergy_Cohen{2}(PrefOriCloserToOri1orOri2(144:end)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==2),NI_Data.denergy_Cohen{2}(PrefOriCloserToOri1orOri2(144:end)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Gamma)');% title ('Monkey 1');
% xlim(xLims); 
% % ylim(yLims);
% set(subplot(3,2,4),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% 
% % Hi-Gamma
% subplot(3,2,5)
% 
% y = oriTuningData.PO{3};
% y1 = oriTuningData.oValsUnique_Plaid;
% y2 = oriTuningData.oValsUnique2_Plaid;
% 
% yminusy1 = y-y1;
% % yminusy1comp = y-(180-y1);
% yminusy2 = y-y2;
% % yminusy2comp = y-(180-y2);
% alldiff = [yminusy1; yminusy2];% yminusy1comp; yminusy2comp];
% 
% for i=1:191
%     clear idx
%     idx = find(abs(alldiff(:,i))==min(abs(alldiff(:,i))),1,'first');
%     OriDiff(i) = alldiff(idx,i);
%     ids(i) = idx;
%     if idx == 1
%         PrefOriCloserToOri1orOri2 (i)= 1; 
%     elseif idx == 2
%         PrefOriCloserToOri1orOri2 (i)= 2; 
%     end
% end
% 
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==1),NI_Data.denergy_Cohen{3}(PrefOriCloserToOri1orOri2(1:143)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==2),NI_Data.denergy_Cohen{3}(PrefOriCloserToOri1orOri2(1:143)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Hi-Gamma)');% title ('Monkey 1');
% % xTicks = [-45 0 45]; yTicks = [1 3]; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
% % xLims = [-50 50]; yLims = [-1 10];
% xlim(xLims); ylim(yLims);
% legend(subplot(3,2,5),{'Ori 1','Ori 2'},'fontSize',12,'Location','best')
% set(subplot(3,2,5),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% subplot(3,2,6)
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==1),NI_Data.denergy_Cohen{3}(PrefOriCloserToOri1orOri2(144:end)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==2),NI_Data.denergy_Cohen{3}(PrefOriCloserToOri1orOri2(144:end)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Hi-Gamma)');% title ('Monkey 1');
% xlim(xLims); 
% % ylim(yLims);
% set(subplot(3,2,6),'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);


%%
% figure(3)
% 
% % Spikes
% subplot(3,2,1)
% 
% y = oriTuningData.PO{1};
% y1 = oriTuningData.oValsUnique_Plaid;
% y2 = oriTuningData.oValsUnique2_Plaid;
% 
% yminusy1 = y-y1;
% % yminusy1comp = y-(180-y1);
% yminusy2 = y-y2;
% % yminusy2comp = y-(180-y2);
% alldiff = [yminusy1; yminusy2];%; yminusy1comp; yminusy2comp];
% 
% for i=1:191
%     clear idx
%     idx = find(abs(alldiff(:,i))==min(abs(alldiff(:,i))),1,'first');
%     OriDiff(i) = alldiff(idx,i);
%     ids(i) = idx;
%     if idx == 1||idx == 3
%         PrefOriCloserToOri1orOri2 (i)= 1; 
%     elseif idx == 2||idx == 4
%         PrefOriCloserToOri1orOri2 (i)= 2; 
%     end
% end
% 
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==1),NI_Data.firingRate_ST_Cohen(PrefOriCloserToOri1orOri2(1:143)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==2),NI_Data.firingRate_ST_Cohen(PrefOriCloserToOri1orOri2(1:143)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Spikes)');% title ('Monkey 1');
% xTicks = [-45 0 45]; yTicks = [1 3]; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
% xLims = [-50 50]; yLims = [-1 5];
% xlim(xLims); ylim(yLims);
% % set(subplot(3,2,1),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% subplot(3,2,2)
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==1),NI_Data.firingRate_ST_Cohen(PrefOriCloserToOri1orOri2(144:end)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==2),NI_Data.firingRate_ST_Cohen(PrefOriCloserToOri1orOri2(144:end)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Spikes)');% title ('Monkey 1');
% xlim(xLims); ylim(yLims);
% % set(subplot(3,2,2),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% 
% % Gamma
% subplot(3,2,3)
% 
% y = oriTuningData.PO{2};
% y1 = oriTuningData.oValsUnique_Plaid;
% y2 = oriTuningData.oValsUnique2_Plaid;
% 
% yminusy1 = y-y1;
% % yminusy1comp = y-(180-y1);
% yminusy2 = y-y2;
% % yminusy2comp = y-(180-y2);
% alldiff = [yminusy1; yminusy2];%; yminusy1comp; yminusy2comp];
% 
% for i=1:191
%     clear idx
%     idx = find(abs(alldiff(:,i))==min(abs(alldiff(:,i))),1,'first');
%     OriDiff(i) = alldiff(idx,i);
%     ids(i) = idx;
%     if idx == 1||idx == 3
%         PrefOriCloserToOri1orOri2 (i)= 1; 
%     elseif idx == 2||idx == 4
%         PrefOriCloserToOri1orOri2 (i)= 2; 
%     end
% end
% 
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==1),NI_Data.energy_ST_Cohen{2}(PrefOriCloserToOri1orOri2(1:143)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==2),NI_Data.energy_ST_Cohen{2}(PrefOriCloserToOri1orOri2(1:143)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Gamma)');% title ('Monkey 1');
% xTicks = [-45 0 45]; yTicks = [1 3]; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
% xLims = [-50 50]; yLims = [-1 10];
% % xlim(xLims); ylim(yLims);
% % set(subplot(3,2,3),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% subplot(3,2,4)
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==1),NI_Data.energy_ST_Cohen{2}(PrefOriCloserToOri1orOri2(144:end)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==2),NI_Data.energy_ST_Cohen{2}(PrefOriCloserToOri1orOri2(144:end)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Gamma)');% title ('Monkey 1');
% % xlim(xLims); ylim(yLims);
% % set(subplot(3,2,4),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% 
% % Hi-Gamma
% subplot(3,2,5)
% 
% y = oriTuningData.PO{3};
% y1 = oriTuningData.oValsUnique_Plaid;
% y2 = oriTuningData.oValsUnique2_Plaid;
% 
% yminusy1 = y-y1;
% % yminusy1comp = y-(180-y1);
% yminusy2 = y-y2;
% % yminusy2comp = y-(180-y2);
% alldiff = [yminusy1; yminusy2];%; yminusy1comp; yminusy2comp];
% 
% for i=1:191
%     clear idx
%     idx = find(abs(alldiff(:,i))==min(abs(alldiff(:,i))),1,'first');
%     OriDiff(i) = alldiff(idx,i);
%     ids(i) = idx;
%     if idx == 1||idx == 3
%         PrefOriCloserToOri1orOri2 (i)= 1; 
%     elseif idx == 2||idx == 4
%         PrefOriCloserToOri1orOri2 (i)= 2; 
%     end
% end
% 
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==1),NI_Data.energy_ST_Cohen{3}(PrefOriCloserToOri1orOri2(1:143)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(1:143)==2),NI_Data.energy_ST_Cohen{3}(PrefOriCloserToOri1orOri2(1:143)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Hi-Gamma)');% title ('Monkey 1');
% xTicks = [-45 0 45]; yTicks = [1 3]; %tickLengthPlot = 2*get(subplot(2,2,1),'TickLength');
% xLims = [-50 50]; yLims = [-1 10];
% % xlim(xLims); ylim(yLims);
% % set(subplot(3,2,5),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% subplot(3,2,6)
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==1),NI_Data.energy_ST_Cohen{3}(PrefOriCloserToOri1orOri2(144:end)==1),'filled'); axis square; hold on;
% scatter(OriDiff(PrefOriCloserToOri1orOri2(144:end)==2),NI_Data.energy_ST_Cohen{3}(PrefOriCloserToOri1orOri2(144:end)==2),'filled'); 
% xlabel('Pref. Ori. - Nearest Ori1/Ori2'); ylabel('N.I. (Hi-Gamma)');% title ('Monkey 1');
% % xlim(xLims); ylim(yLims);
% % set(subplot(3,2,6),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% 
% % scatter(oriTuningData.PO{2}(1:143),),'k','filled'); axis square
% % xlabel('OS (LFP Gamma)'); ylabel('OS (spikes)');
% % xlim(xLims); ylim(yLims);
% % set(subplot(2,2,4),'xTick',xTicks,'yTick',yTicks,'TickDir','out','Ticklength',tickLengthPlot,'box','off','fontSize',fontSize);
% 
% figure(2)
% histogram(oriTuningData.PO{1})
% hold on;
% histogram((oriTuningData.PO{2}))
end
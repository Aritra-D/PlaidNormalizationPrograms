function displayEyeDataIndividualMonkey(monkeyName,folderSourceString,stimType,timeRange,cutOff,colorScheme)

close all; % closes any open figure to avoid any overlaying issues

hFigure1 = figure(1); % Eye data DegX
set(hFigure1,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig1.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.6],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot1);
hPlotsFig1.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.1125],0.01,0.01,1); linkaxes(hPlotsFig1.hPlot2);

hFigure2 = figure(2); %Eye data DegY
set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig2.hPlot1 = getPlotHandles(5,5,[0.2 0.3 0.6 0.6],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot1);
hPlotsFig2.hPlot2 = getPlotHandles(1,5,[0.2 0.08 0.6 0.1125],0.01,0.01,1); linkaxes(hPlotsFig2.hPlot2);

folderSourceString_Project = strtok(folderSourceString,'\');
folderSave = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures\EyeData\');
if ~exist(folderSave,'dir')
    mkdir(folderSave)
end

fileSave = fullfile(folderSave,['eyeData_' monkeyName '_T_' num2str(1000*timeRange(1)) '_ms_' num2str(1000*timeRange(2)) '_ms_CutOff_', num2str(cutOff)]);
if exist(fileSave,'file')
    disp(['Loading file ' fileSave]);
    load(fileSave); %#ok<LOAD>
else
    [eyeDataDegX,eyeDataDegY,eyeData,eyeSpeedData,microsaccadeData,timeValsEyeData,trialNums,FsEye] = getEyeDataIndividualMonkey(monkeyName,folderSourceString,timeRange,cutOff);
    save(fileSave,'eyeDataDegX','eyeDataDegY','eyeData','eyeSpeedData','microsaccadeData','timeValsEyeData','trialNums','FsEye')
end

cValsUnique = [0 6.25 12.5 25 50];
cValsUnique2 = cValsUnique;
if strcmp(colorScheme,'color')
    colors = jet(length(cValsUnique));
elseif strcmp(colorScheme,'grayscale')
    colors = repmat(0.85:-0.1:0.45,[3 1])';
end
cFlipped_Indices = flip(1:length(cValsUnique2));

tickLengthPlot = 2*get(hPlotsFig1.hPlot1(1),'TickLength');

% Plot Eye-Position Horizontal (deg) and Eye-Position vertical (deg)
for c_Ori2 = 1:5
    for c_Ori1 = 1:5
        if strcmp(stimType,'static')
            plot(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),timeValsEyeData,mean(eyeDataDegX.static{c_Ori2,c_Ori1},1),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
            plot(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),timeValsEyeData,mean(eyeDataDegY.static{c_Ori2,c_Ori1},1),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
            
        elseif strcmp(stimType,'CP')
            plot(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),timeValsEyeData,mean(eyeDataDegX.CP{c_Ori2,c_Ori1},1),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
            plot(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),timeValsEyeData,mean(eyeDataDegY.CP{c_Ori2,c_Ori1},1),'color',colors(cFlipped_Indices(c_Ori2),:,:),'LineWidth',2)
        end
        set(hPlotsFig1.hPlot1(c_Ori2,c_Ori1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
        set(hPlotsFig2.hPlot1(c_Ori2,c_Ori1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    end
end

for cOri1 = 1:5
    title(hPlotsFig1.hPlot1(1,cOri1),[num2str(cValsUnique(cOri1)) ' %'],'fontsize',16)
    title(hPlotsFig2.hPlot1(1,cOri1),[num2str(cValsUnique(cOri1)) ' %'],'fontsize',16)
end


textString = {'50 %','25 %','12.5 %','6.25 %','0 %'};
for i=1:2
    figure(i)
    
    for cOri2 = 1:5
        textH{cOri2} = getPlotHandles(1,1,[0.81 0.82-(cOri2-1)*0.115 0.01 0.01]); %#ok<AGROW>
        set(textH{cOri2},'Visible','Off')
        text(0.35,1.15,textString{cOri2},'unit','normalized','fontsize',16,'fontweight','bold','parent',textH{cOri2});
    end
    
    textString1 ='Contrast of Grating 1';  textString2 = 'Contrast of Grating 2';
    textH1 = getPlotHandles(1,1,[0.36 0.97 0.01 0.01]);
    textH2 = getPlotHandles(1,1,[0.88 0.85 0.01 0.01]);
    set(textH1,'Visible','Off');set(textH2,'Visible','Off');
    text(0.35,1.15,textString1,'unit','normalized','fontsize',16,'fontweight','bold','parent',textH1);
    text(0.35,1.15,textString2,'unit','normalized','fontsize',16,'fontweight','bold','rotation',270,'parent',textH2);
end
rescaleData(hPlotsFig1.hPlot1,-0.5,0.5,[-0.5 0.5]+ getYLims(hPlotsFig1.hPlot1),14);
rescaleData(hPlotsFig2.hPlot1,-0.5,0.5,[-0.5 0.5]+ getYLims(hPlotsFig1.hPlot1),14);

for c_Ori2 = 1: length(cValsUnique2)
    for c_Ori1 = 1:length(cValsUnique)
        if strcmp(stimType,'static')
            plot(hPlotsFig1.hPlot2(1,c_Ori1),timeValsEyeData,mean(eyeDataDegX.static{cFlipped_Indices(c_Ori2),c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
            plot(hPlotsFig2.hPlot2(1,c_Ori1),timeValsEyeData,mean(eyeDataDegY.static{cFlipped_Indices(c_Ori2),c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
            
        elseif strcmp(stimType,'CP')
            plot(hPlotsFig1.hPlot2(1,c_Ori1),timeValsEyeData,mean(eyeDataDegX.CP{cFlipped_Indices(c_Ori2),c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
            plot(hPlotsFig2.hPlot2(1,c_Ori1),timeValsEyeData,mean(eyeDataDegY.CP{cFlipped_Indices(c_Ori2),c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
        end
        hold(hPlotsFig1.hPlot2(1,c_Ori1),'on');
        hold(hPlotsFig2.hPlot2(1,c_Ori1),'on');
        %         plot(hPlotsFig1.hPlot2(2,c_Ori1),firingRateData.timeVals,squeeze(combinedData{1}(cFlipped_Indices(c_Ori2),c_Ori1,:)),'color',colors(c_Ori2,:,:),'LineWidth',2);
        %         hold(hPlotsFig1.hPlot2(2,c_Ori1),'on');
    end
end

for i = 1:length(cValsUnique)
    set(hPlotsFig1.hPlot2(1,i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    %     set(hPlotsFig1.hPlot2(2,i),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end

xlabel(hPlotsFig1.hPlot1(5,1),'Time (s)'); ylabel(hPlotsFig1.hPlot1(5,1),[{'Eye-Position'} {'Horizontal (deg)'}]);
xlabel(hPlotsFig2.hPlot1(5,1),'Time (s)'); ylabel(hPlotsFig2.hPlot1(5,1),[{'Eye-Position'} {'Vertical (deg)'}]);

xlabel(hPlotsFig1.hPlot2(1,1),'Time (s)'); ylabel(hPlotsFig1.hPlot2(1,1),[{'Eye-Position'} {'Horizontal (deg)'}]);
xlabel(hPlotsFig2.hPlot2(1,1),'Time (s)'); ylabel(hPlotsFig2.hPlot2(1,1),[{'Eye-Position'} {'Vertical (deg)'}]);

rescaleData(hPlotsFig1.hPlot2,-0.5,0.5,[-0.5 0.5],14);
rescaleData(hPlotsFig2.hPlot2,-0.5,0.5,[-0.5 0.5],14);

rescaleData(hPlotsFig1.hPlot1,-0.5,0.5,[-0.5 0.5],14);
rescaleData(hPlotsFig2.hPlot1,-0.5,0.5,[-0.5 0.5],14);

if strcmp(colorScheme,'color')
    folderSaveFigs = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\Figures\Color');
elseif strcmp(colorScheme,'greyscale')||strcmp(colorScheme,'grayscale')
    folderSaveFigs = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\Figures\GrayScale');
end

FigName1 = fullfile(folderSaveFigs,['Figure_' monkeyName '_' stimType 'EyePosition_Horizontal']);
FigName2 = fullfile(folderSaveFigs,['Figure_' monkeyName '_' stimType 'EyePosition_Vertical']);
saveas(hFigure1,[FigName1 '.fig'])
saveas(hFigure1,[FigName1,'.tif'])
saveas(hFigure2,[FigName2 '.fig'])
saveas(hFigure2,[FigName2,'.tif'])
print(hFigure1,[FigName1,'HQ.tif'],'-dtiff','-r300')
print(hFigure2,[FigName2,'HQ.tif'],'-dtiff','-r300')


end

% Draw lines for timeTange or FreqRange
function displayRange(plotHandles,range,yLims,colorName) %#ok<DEFNU>
[nX,nY] = size(plotHandles);
%yLims = getYLims(plotHandles);

yVals = yLims(1):(yLims(2)-yLims(1))/100:yLims(2);
xVals1 = range(1) + zeros(1,length(yVals));
xVals2 = range(2) + zeros(1,length(yVals));

for i=1:nX
    for j=1:nY
        hold(plotHandles(i,j),'on');
        plot(plotHandles(i,j),xVals1,yVals,'color',colorName,'LineWidth',1.5);
        plot(plotHandles(i,j),xVals2,yVals,'color',colorName,'LineWidth',1.5);
    end
end
end

% get Y limits for an axis
function yLims = getYLims(plotHandles)

[numRows,numCols] = size(plotHandles);
% Initialize
yMin = inf;
yMax = -inf;

for row=1:numRows
    for column=1:numCols
        % get positions
        axis(plotHandles(row,column),'tight');
        tmpAxisVals = axis(plotHandles(row,column));
        if tmpAxisVals(3) < yMin
            yMin = tmpAxisVals(3);
        end
        if tmpAxisVals(4) > yMax
            yMax = tmpAxisVals(4);
        end
    end
end

yLims=[yMin yMax];
end

% Rescale data
function rescaleData(plotHandles,xMin,xMax,yLims,labelSize)

[numRows,numCols] = size(plotHandles);
% labelSize=14;
for i=1:numRows
    for j=1:numCols
        hold(plotHandles(i,j),'on');
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j==1
                set(plotHandles(i,j),'fontSize',labelSize);
            elseif j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        elseif (rem(i,2)==0 && j~=1)
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        elseif (rem(i,2)==1 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
        elseif (rem(i,2)==1 && j~=1)
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        else
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
% set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
% set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
% set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
% set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
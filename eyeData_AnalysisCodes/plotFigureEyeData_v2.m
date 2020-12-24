function plotFigureEyeData_v2(monkeyName,folderSourceString,timeRange,cutOff,colorScheme)
close all; % closes any open figure to avoid any overlaying issues

hFigure1 = figure(1); % Eye data DegX
set(hFigure1,'units','normalized','outerposition',[0 0 1 1])
hPlotsFig1.hPlot1 = getPlotHandles(3,5,[0.15 0.55 0.6 0.4],0.01,0.01,1); 
hPlotsFig1.hPlot2 = getPlotHandles(3,5,[0.15 0.08 0.6 0.4],0.01,0.01,1); 
hPlotsFig1.hPlot3 = getPlotHandles(3,1,[0.8 0.55 0.11 0.4],0.01,0.01,1); 
hPlotsFig1.hPlot4 = getPlotHandles(3,1,[0.8 0.08 0.11 0.4],0.01,0.01,1); 

for i=3
    linkaxes(hPlotsFig1.hPlot1(i,(1:5)))
    linkaxes(hPlotsFig1.hPlot2(i,(1:5)))
end

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
    colors = flip(jet(length(cValsUnique)));
elseif strcmp(colorScheme,'grayscale')
    colors = flip(repmat(0.85:-0.1:0.45,[3 1])');
end
% cFlipped_Indices = flip(1:length(cValsUnique2));

count = 1; %cLen = length(cValsUnique2);
% Plot Eye-Position Horizontal (deg) and Eye-Position vertical (deg)
for c_Ori2 = length(cValsUnique2):-1:1
    for c_Ori1 = 1:length(cValsUnique)
        plot(hPlotsFig1.hPlot1(1,c_Ori1),timeValsEyeData,mean(eyeDataDegX.static{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig1.hPlot1(1,c_Ori1),'on');
        plot(hPlotsFig1.hPlot1(2,c_Ori1),timeValsEyeData,mean(eyeDataDegY.static{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig1.hPlot1(2,c_Ori1),'on');
        clear H
        [H,timeValsMS] = getPSTH(microsaccadeData.static{c_Ori2, c_Ori1},20,[timeValsEyeData(1) timeValsEyeData(end)]);
        plot(hPlotsFig1.hPlot1(3,c_Ori1),timeValsMS,H,'color',colors(c_Ori2,:,:),'LineWidth',2);   
        hold(hPlotsFig1.hPlot1(3,c_Ori1),'on');
        
        clear mXDataAcrossTimePoints mYDataAcrossTimePoints mMSCountAcrossTimePoints
        mXDataAcrossTimePoints = mean(eyeDataDegX.static{c_Ori2,c_Ori1},2);
        mXDataAcrossStims = mean(mXDataAcrossTimePoints);
        semXDataAcrossStims = std(mXDataAcrossTimePoints,[],1)/sqrt(size(mXDataAcrossTimePoints,1));
        XData{1}{count} = mXDataAcrossTimePoints;
        errorbar(hPlotsFig1.hPlot3(1,1),count,mXDataAcrossStims,semXDataAcrossStims,'-o','color',colors(c_Ori2,:,:));
        hold(hPlotsFig1.hPlot3(1,1),'on');

        mYDataAcrossTimePoints = mean(eyeDataDegY.static{c_Ori2,c_Ori1},2);
        mYDataAcrossStims = mean(mYDataAcrossTimePoints);
        semYDataAcrossStims = std(mYDataAcrossTimePoints,[],1)/sqrt(size(mYDataAcrossTimePoints,1));
        YData{1}{count} = mYDataAcrossTimePoints;
        errorbar(hPlotsFig1.hPlot3(2,1),count,mYDataAcrossStims,semYDataAcrossStims,'-o','color',colors(c_Ori2,:,:));
        hold(hPlotsFig1.hPlot3(2,1),'on');

        
        mMSCountAcrossTimePoints = microsaccadeData.static_Count{c_Ori2,c_Ori1};
        mMSCountAcrossStims = mean(mMSCountAcrossTimePoints,2);
        semMSDataAcrossStims = std(mMSCountAcrossTimePoints,[],2)/sqrt(size(mMSCountAcrossTimePoints,2));
        MSData{1}{count} = mMSCountAcrossTimePoints;
        errorbar(hPlotsFig1.hPlot3(3,1),count,mMSCountAcrossStims,semMSDataAcrossStims,'-o','color',colors(c_Ori2,:,:));
        hold(hPlotsFig1.hPlot3(3,1),'on');

        
        plot(hPlotsFig1.hPlot2(1,c_Ori1),timeValsEyeData,mean(eyeDataDegX.CP{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig1.hPlot2(1,c_Ori1),'on');
        plot(hPlotsFig1.hPlot2(2,c_Ori1),timeValsEyeData,mean(eyeDataDegY.CP{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
        hold(hPlotsFig1.hPlot2(2,c_Ori1),'on');
        clear H
        [H,timeValsMS] = getPSTH(microsaccadeData.CP{c_Ori2, c_Ori1},20,[timeValsEyeData(1) timeValsEyeData(end)]);
        plot(hPlotsFig1.hPlot2(3,c_Ori1),timeValsMS,H,'color',colors(c_Ori2,:,:),'LineWidth',2);   
        hold(hPlotsFig1.hPlot2(3,c_Ori1),'on');
        
        clear mXDataAcrossTimePoints mYDataAcrossTimePoints mMSCountAcrossTimePoints
        mXDataAcrossTimePoints = mean(eyeDataDegX.CP{c_Ori2,c_Ori1},2);
        mXDataAcrossStims = mean(mXDataAcrossTimePoints);
        semXDataAcrossStims = std(mXDataAcrossTimePoints,[],1)/sqrt(size(mXDataAcrossTimePoints,1));
        XData{2}{count} = mXDataAcrossTimePoints;
        errorbar(hPlotsFig1.hPlot4(1,1),count,mXDataAcrossStims,semXDataAcrossStims,'-o','color',colors(c_Ori2,:,:));
        hold(hPlotsFig1.hPlot4(1,1),'on');


        mYDataAcrossTimePoints = mean(eyeDataDegY.CP{c_Ori2,c_Ori1},2);
        mYDataAcrossStims = mean(mYDataAcrossTimePoints);
        semYDataAcrossStims = std(mYDataAcrossTimePoints,[],1)/sqrt(size(mYDataAcrossTimePoints,1));
        YData{2}{count} = mYDataAcrossTimePoints;
        errorbar(hPlotsFig1.hPlot4(2,1),count,mYDataAcrossStims,semYDataAcrossStims,'-o','color',colors(c_Ori2,:,:));
        hold(hPlotsFig1.hPlot4(2,1),'on');

        
        mMSCountAcrossTimePoints = microsaccadeData.CP_Count{c_Ori2,c_Ori1};
        mMSCountAcrossStims = mean(mMSCountAcrossTimePoints,2);
        semMSDataAcrossStims = std(mMSCountAcrossTimePoints,[],2)/sqrt(size(mMSCountAcrossTimePoints,2));
        MSData{2}{count} = mMSCountAcrossTimePoints;
        errorbar(hPlotsFig1.hPlot4(3,1),count,mMSCountAcrossStims,semMSDataAcrossStims,'-o','color',colors(c_Ori2,:,:));
        hold(hPlotsFig1.hPlot4(3,1),'on');

        count = count+1;
    end
end

% One-Way Anova for 25 contrast conditions

% data grouping
for i=1:2
    XData_CombinedAcrossContrasts{i}=[]; %#ok<*AGROW>
    YData_CombinedAcrossContrasts{i}=[];
    MSData_CombinedAcrossContrasts{i}=[];
    contrastFactor{i} =[];
    for j=1:length(XData{1,1})
        XData_CombinedAcrossContrasts{i} = cat(1,XData_CombinedAcrossContrasts{i},XData{i}{j});
        YData_CombinedAcrossContrasts{i} = cat(1,YData_CombinedAcrossContrasts{i},YData{i}{j});
        MSData_CombinedAcrossContrasts{i} = cat(2,MSData_CombinedAcrossContrasts{i},MSData{i}{j});
        contrastConditionTMP = j*ones(1,size(XData{i}{j},1));
        contrastFactor{i} = cat(2,contrastFactor{i},contrastConditionTMP);
    end
end

% Statistical test on Factor of contrast across 25 conditions
[p_XData_static,tbl_XData_static,stats_XData_static] = anovan(XData_CombinedAcrossContrasts{1,1},{contrastFactor{1,1}},'display','off');
[p_XData_CP,tbl_XData_CP,stats_XData_CP] = anovan(XData_CombinedAcrossContrasts{1,2},{contrastFactor{1,2}},'display','off');
% 
[p_YData_static,tbl_YData_static,stats_YData_static] = anovan(YData_CombinedAcrossContrasts{1,1},{contrastFactor{1,1}},'display','off');
[p_YData_CP,tbl_YData_CP,stats_YData_CP] = anovan(YData_CombinedAcrossContrasts{1,2},{contrastFactor{1,2}},'display','off');

[p_MSData_static,tbl_MSData_static,stats_MSData_static] = anovan(MSData_CombinedAcrossContrasts{1,1},{contrastFactor{1,1}},'display','off');
[p_MSData_CP,tbl_MSData_CP,stats_MSData_CP] = anovan(MSData_CombinedAcrossContrasts{1,2},{contrastFactor{1,2}},'display','off');

tickLengthPlot = 2*get(hPlotsFig1.hPlot1(1),'TickLength');
fontSize = 14;
for i = 1:length(cValsUnique)
    set(hPlotsFig1.hPlot1(:,i),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlotsFig1.hPlot2(:,i),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end

for i=1:3
    set(hPlotsFig1.hPlot3(i,1),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlotsFig1.hPlot4(i,1),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end


rescaleData(hPlotsFig1.hPlot1(1,1:5),timeRange(1),timeRange(2),[-0.5 0.5],14);
yTicks = [-0.3 0 0.3];set(hPlotsFig1.hPlot1(1,1),'yTick',yTicks,'yTicklabel',yTicks);
rescaleData(hPlotsFig1.hPlot1(2,1:5),timeRange(1),timeRange(2),[-0.5 0.5],14);
set(hPlotsFig1.hPlot1(2,1),'yTick',[-0.3 0 0.3],'yTicklabel',[-0.3 0 0.3]);
rescaleData(hPlotsFig1.hPlot1(3,1:5),timeRange(1),timeRange(2),[0,4],14);
set(hPlotsFig1.hPlot1(3,1:2:5),'xTicklabel',[-0.5 0 0.5]);
set(hPlotsFig1.hPlot2(3,1),'yTick',[1 3],'yTicklabel',[1 3]);

rescaleData(hPlotsFig1.hPlot2(1,1:5),-0.5,0.5,[-0.5 0.5],14);
set(hPlotsFig1.hPlot2(1,1),'yTick',[-0.3 0 0.3],'yTicklabel',[-0.3 0 0.3]);
rescaleData(hPlotsFig1.hPlot2(2,1:5),-0.5,0.5,[-0.5 0.5],14);
set(hPlotsFig1.hPlot2(2,1),'yTick',[-0.3 0 0.3],'yTicklabel',[-0.3 0 0.3]);
rescaleData(hPlotsFig1.hPlot2(3,1:5),-0.5,0.5,[0,4],14);
set(hPlotsFig1.hPlot2(3,1:2:5),'xTicklabel',[-0.5 0 0.5]);
set(hPlotsFig1.hPlot2(3,1),'yTick',[1 3],'yTicklabel',[1 3]);


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
function yLims = getYLims(plotHandles) %#ok<DEFNU>

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
                set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
            elseif j~=1
                set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
            end
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
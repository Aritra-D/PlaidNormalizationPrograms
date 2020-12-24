function plotFigureEyeData(monkeyName,colorScheme)
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

folderSourceString= 'E:\';
dataFolder = fullfile(folderSourceString,'data\PlaidNorm\');
folderEyeData = fullfile(folderSourceString,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures\eyeData\',monkeyName);

eyeDataFile = fullfile(folderEyeData,'eyeData.mat');
microsaccadesFile = fullfile(folderEyeData,'microsaccades.mat');
microsaccadesStatsFile = fullfile(folderEyeData,'microsaccadesStats.mat');
protNums = fullfile(folderEyeData,'protNumFortrial.mat');

if exist(eyeDataFile,'file')&& exist(microsaccadesFile,'file')&& exist(microsaccadesStatsFile,'file')
    disp(['Loading file ' eyeDataFile]);disp(['Loading file ' microsaccadesFile]);disp(['Loading file ' microsaccadesStatsFile])
    eyeData = load(eyeDataFile); microsaccadeData = load(microsaccadesFile); microsaccadeStats = load(microsaccadesStatsFile);  %#ok<NASGU>
    load(protNums);
else
    error('eyeData files not found!')
end

cValsUnique = [0 6.25 12.5 25 50];
cValsUnique2 = cValsUnique;
if strcmp(colorScheme,'color')
    colors = flip(jet(length(cValsUnique)));
elseif strcmp(colorScheme,'grayscale')
    colors = flip(repmat(0.85:-0.1:0.45,[3 1])');
end



gridType = 'microelectrode';
[expDates, protocolNames,~,~,~] = dataInformationPlaidNorm(monkeyName,gridType,0);

badTrialsFolder = 'E:\Projects\Aritra_PlaidNormalizationProject\badTrials';

%% Segregating eyeDataDegX and eyeDataDegY for different contrast conditions for different sessions

for iProt = 1:length(expDates)
    folderExtract = fullfile(dataFolder,'data',monkeyName,gridType,expDates{iProt},protocolNames{iProt},'extractedData');
    %     folderSegment = fullfile(folderSourceString,'data',monkeyName,gridType,expDates{iProt},protocolNames{iProt},'segmentedData');
    % Get Combinations
    [parameterCombinations,parameterCombinations2,...
        aValsUnique,eValsUnique,~,~,oValsUnique,cValsUnique,tValsUnique, ...
        aValsUnique2,eValsUnique2,~,~,oValsUnique2,cValsUnique2,tValsUnique2] = loadParameterCombinations(folderExtract); %#ok<*ASGLU>
    
    
    cListFlipped_Ori2 = flip(1:length(cValsUnique2)); % helps in plotting the responses from low to high contrast
    
    if aValsUnique ~= aValsUnique2 || eValsUnique ~= eValsUnique2
        error('Azimuths and/or elevations do not match!');
    end
    
    if tValsUnique ~= tValsUnique2
        error('TF do not match!');
    else
        tList = 1:length(tValsUnique);
    end
    
    a=1; e=1; s=1; f=1; o=1;
    % Get bad trials
    badTrialFile = fullfile(badTrialsFolder,[monkeyName, expDates{iProt} protocolNames{iProt} '_badTrials.mat']);
    if ~exist(badTrialFile,'file')
        disp('Bad trial file does not exist...');
        badTrials=[];
    else
        disp(['Loading' badTrialFile])
        badTrials = loadBadTrials(badTrialFile);
        disp([num2str(length(badTrials)) ' bad trials']);
    end
    
    
    for t = 1:length(tList)
        for c_Ori2 = 1:length(cValsUnique2)
            for c_Ori1 = 1:length(cValsUnique)
                clear goodPos  goodPosIDs eyeDataDegX_TMP eyeDataDegY_TMP numMS_TMP
                
                goodPos = parameterCombinations{a,e,s,f,o,c_Ori1,tList(t)};
                goodPos2 = parameterCombinations2{a,e,s,f,o,cListFlipped_Ori2(c_Ori2),tList(t)};
                goodPos = intersect(goodPos,goodPos2);
                goodPos = setdiff(goodPos,badTrials); % discarding bad trials of LFPs
                
                if isempty(goodPos)
                    disp('No entries for this combination..');
                else
                    disp(['iProt: ' num2str(iProt) ' t: ' num2str(t) ' cOri2: ' num2str(c_Ori2) ' cOri1: ' num2str(c_Ori1)])
                    goodPos = intersect(goodPos,microsaccadeData.trialNumsForMSAnalysisPerProtocol{iProt}); %finding intersect of goodEyeTrials and goodStims
                    trialNums{iProt}{t}(c_Ori2,c_Ori1) = length(goodPos);
                    
                    for i = 1:length(goodPos)
                        goodPosIDs(i) = find(goodPos(i)==microsaccadeData.trialNumsForMSAnalysisPerProtocol{iProt});
                    end
                    eyeDataDegX_TMP = eyeData.eyeDataDegX(protNumForTrial==iProt,:);
                    eyeDataDegX{iProt}{t}{c_Ori2,c_Ori1} = eyeDataDegX_TMP(goodPosIDs,:);
                    
                    eyeDataDegY_TMP = eyeData.eyeDataDegY(protNumForTrial==iProt,:);
                    eyeDataDegY{iProt}{t}{c_Ori2,c_Ori1} = eyeDataDegY_TMP(goodPosIDs,:);
                    
                    MS_TMP = microsaccadeData.MS(protNumForTrial==iProt);
                    MS{iProt}{t}{c_Ori2,c_Ori1} = MS_TMP(goodPosIDs);
                    
                    numMS_TMP = microsaccadeData.numMSInRange(protNumForTrial==iProt);
                    numMS{iProt}{t}{c_Ori2,c_Ori1} = numMS_TMP(goodPosIDs);
                end
            end
        end
    end
end

%% Combining eyeDataDegX and eyeDataDegY for different contrast conditions across different sessions
eyeDataDegX_Combined = eyeDataDegX{1}; eyeDataDegY_Combined = eyeDataDegY{1};
trialNums_Combined = trialNums{1}; numMS_Combined = numMS{1}; MS_Combined = MS{1};
for iProt = 2:length(expDates)
    for t = 1:length(tList)
        for c_Ori2 = 1:length(cValsUnique2)
            for c_Ori1 = 1:length(cValsUnique)
                eyeDataDegX_Combined{t}{c_Ori2,c_Ori1} = cat(1,eyeDataDegX_Combined{t}{c_Ori2,c_Ori1},eyeDataDegX{iProt}{t}{c_Ori2,c_Ori1});
                eyeDataDegY_Combined{t}{c_Ori2,c_Ori1} = cat(1,eyeDataDegY_Combined{t}{c_Ori2,c_Ori1},eyeDataDegY{iProt}{t}{c_Ori2,c_Ori1});
                MS_Combined{t}{c_Ori2,c_Ori1} = cat(2,MS_Combined{t}{c_Ori2,c_Ori1},MS{iProt}{t}{c_Ori2,c_Ori1});
                
                numMS_Combined{t}{c_Ori2,c_Ori1} = cat(2,numMS_Combined{t}{c_Ori2,c_Ori1},numMS{iProt}{t}{c_Ori2,c_Ori1});
            end
        end
        trialNums_Combined{t}= trialNums_Combined{t}+ trialNums{iProt}{t};
    end
end


%%         % Grand Mean of eyeDegX,eyeDegY and numMS
grandmean_eyeDegX = zeros(length(t),length(cValsUnique2),length(cValsUnique));
sem_eyeDegX= zeros(length(t),length(cValsUnique2),length(cValsUnique));
grandmean_eyeDegY = zeros(length(t),length(cValsUnique2),length(cValsUnique));
sem_eyeDegY= zeros(length(t),length(cValsUnique2),length(cValsUnique));
grandmean_numMS = zeros(length(t),length(cValsUnique2),length(cValsUnique));
sem_numMS = zeros(length(t),length(cValsUnique2),length(cValsUnique));


for t=1:length(tList)
    count =1;
    for c_Ori2 = 1:length(cValsUnique2)
        for c_Ori1 = 1:length(cValsUnique)
            clear m_eyeDegX_acrossTimePoints m_eyeDegY_acrossTimePoints m_numMS_acrossTimePoints
            m_eyeDegX_acrossTimePoints = mean(eyeDataDegX_Combined{t}{c_Ori2,c_Ori1},2);
            grandmean_eyeDegX(t,c_Ori2,c_Ori1) = mean(m_eyeDegX_acrossTimePoints,1);
            sem_eyeDegX(t,c_Ori2,c_Ori1) = std(m_eyeDegX_acrossTimePoints,[],1)/sqrt(size(m_eyeDegX_acrossTimePoints,1));
            XData{t}{count} = m_eyeDegX_acrossTimePoints;
            
            m_eyeDegY_acrossTimePoints = mean(eyeDataDegY_Combined{t}{c_Ori2,c_Ori1},2);
            grandmean_eyeDegY(t,c_Ori2,c_Ori1) = mean(m_eyeDegY_acrossTimePoints,1);
            sem_eyeDegY(t,c_Ori2,c_Ori1) = std(m_eyeDegY_acrossTimePoints,[],1)/sqrt(size(m_eyeDegY_acrossTimePoints,1));
            YData{t}{count} = m_eyeDegY_acrossTimePoints;

            
            m_numMS_acrossTimePoints = numMS_Combined{t}{c_Ori2,c_Ori1};
            grandmean_numMS(t,c_Ori2,c_Ori1) = mean(m_numMS_acrossTimePoints,2);
            sem_numMS(t,c_Ori2,c_Ori1) = std(m_numMS_acrossTimePoints,[],2)/sqrt(size(m_numMS_acrossTimePoints,2));
            MSData{t}{count} = m_numMS_acrossTimePoints;
            count= count+1;
        end
    end
end

%% Plots
eyeRangeMS = [-495 495]; FsEye = 200;
timeValsEyePos = (eyeRangeMS(1):1000/FsEye:eyeRangeMS(2)-1000/FsEye)/1000;

for t = 1:length(tList)
    count = 1;
    for c_Ori2 = length(cValsUnique2):-1:1
        for c_Ori1 = 1:length(cValsUnique)
            
            
            if t == 1 % static
                % Plot Eye-Position Horizontal (deg)
                plot(hPlotsFig1.hPlot1(1,c_Ori1),timeValsEyePos,mean(eyeDataDegX_Combined{t}{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
                hold(hPlotsFig1.hPlot1(1,c_Ori1),'on');
                % Plot Eye-Position vertical (deg)
                plot(hPlotsFig1.hPlot1(2,c_Ori1),timeValsEyePos,mean(eyeDataDegX_Combined{t}{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
                hold(hPlotsFig1.hPlot1(2,c_Ori1),'on');
                % Plot MS psth
                clear H timeValsMS
                [H,timeValsMS] = getPSTH(MS_Combined{t}{c_Ori2, c_Ori1},20,[timeValsEyePos(t) timeValsEyePos(end)]);
                plot(hPlotsFig1.hPlot1(3,c_Ori1),timeValsMS,H,'color',colors(c_Ori2,:,:),'LineWidth',2);
                hold(hPlotsFig1.hPlot1(3,c_Ori1),'on');
                
                                
                errorbar(hPlotsFig1.hPlot3(1,1),count,grandmean_eyeDegX(t,c_Ori2,c_Ori1),sem_eyeDegX(t,c_Ori2,c_Ori1),'-o','color', colors(c_Ori2,:,:)); 
                hold(hPlotsFig1.hPlot3(1,1),'on');
                
                errorbar(hPlotsFig1.hPlot3(2,1),count,grandmean_eyeDegY(t,c_Ori2,c_Ori1),sem_eyeDegY(t,c_Ori2,c_Ori1),'-o','color', colors(c_Ori2,:,:));  hold on;
                hold(hPlotsFig1.hPlot3(2,1),'on');
                
                errorbar(hPlotsFig1.hPlot3(3,1),count,grandmean_numMS(t,c_Ori2,c_Ori1),sem_numMS(t,c_Ori2,c_Ori1),'-o','color', colors(c_Ori2,:,:));  hold on;
                hold(hPlotsFig1.hPlot3(3,1),'on');

            end
            
            
            if t == 2 % Counter-phase flickers
                % Plot Eye-Position Horizontal (deg)
                plot(hPlotsFig1.hPlot2(1,c_Ori1),timeValsEyePos,mean(eyeDataDegX_Combined{t}{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
                hold(hPlotsFig1.hPlot2(1,c_Ori1),'on');
                % Plot Eye-Position vertical (deg)
                plot(hPlotsFig1.hPlot2(2,c_Ori1),timeValsEyePos,mean(eyeDataDegX_Combined{t}{c_Ori2,c_Ori1},1),'color',colors(c_Ori2,:,:),'LineWidth',2);
                hold(hPlotsFig1.hPlot2(2,c_Ori1),'on');
                % Plot MS psth
                clear H
                [H,timeValsMS] = getPSTH(MS_Combined{t}{c_Ori2, c_Ori1},20,[timeValsEyePos(1) timeValsEyePos(end)]);
                plot(hPlotsFig1.hPlot2(3,c_Ori1),timeValsMS,H,'color',colors(c_Ori2,:,:),'LineWidth',2);
                hold(hPlotsFig1.hPlot2(3,c_Ori1),'on');

                errorbar(hPlotsFig1.hPlot4(1,1),count,grandmean_eyeDegX(t,c_Ori2,c_Ori1),sem_eyeDegX(t,c_Ori2,c_Ori1),'-o','color', colors(c_Ori2,:,:));
                hold(hPlotsFig1.hPlot4(1,1),'on');

                errorbar(hPlotsFig1.hPlot4(2,1),count,grandmean_eyeDegY(t,c_Ori2,c_Ori1),sem_eyeDegY(t,c_Ori2,c_Ori1),'-o','color', colors(c_Ori2,:,:)); 
                hold(hPlotsFig1.hPlot4(2,1),'on');

                errorbar(hPlotsFig1.hPlot4(3,1),count,grandmean_numMS(t,c_Ori2,c_Ori1),sem_numMS(t,c_Ori2,c_Ori1),'-o','color', colors(c_Ori2,:,:)); 
                hold(hPlotsFig1.hPlot4(3,1),'on');

            end
            count = count+1;
        end
    end
end


%% One-Way Anova for 25 contrast conditions

% data grouping
for t=1:2 % Two TFs
    XData_CombinedAcrossContrasts{t}=[]; %#ok<*AGROW>
    YData_CombinedAcrossContrasts{t}=[];
    MSData_CombinedAcrossContrasts{t}=[];
    contrastFactor{t} =[];
    for j=1:length(XData{1,1})
        XData_CombinedAcrossContrasts{t} = cat(1,XData_CombinedAcrossContrasts{t},XData{t}{j});
        YData_CombinedAcrossContrasts{t} = cat(1,YData_CombinedAcrossContrasts{t},YData{t}{j});
        MSData_CombinedAcrossContrasts{t} = cat(2,MSData_CombinedAcrossContrasts{t},MSData{t}{j});
        contrastConditionTMP = j*ones(1,size(XData{t}{j},1));
        contrastFactor{t} = cat(2,contrastFactor{t},contrastConditionTMP);
    end
end

% Statistical test on Factor of contrast across 25 conditions
[p_XData_static,tbl_XData_static,stats_XData_static] = anovan(XData_CombinedAcrossContrasts{1,1},{contrastFactor{1,1}},'display','off');
disp(['Anova Test: EyeDataDegX.Static, p= ' num2str(p_XData_static)])
[p_XData_CP,tbl_XData_CP,stats_XData_CP] = anovan(XData_CombinedAcrossContrasts{1,2},{contrastFactor{1,2}},'display','off');
disp(['Anova Test: EyeDataDegX.CP, p= ' num2str(p_XData_CP)])

%
[p_YData_static,tbl_YData_static,stats_YData_static] = anovan(YData_CombinedAcrossContrasts{1,1},{contrastFactor{1,1}},'display','off');
[p_YData_CP,tbl_YData_CP,stats_YData_CP] = anovan(YData_CombinedAcrossContrasts{1,2},{contrastFactor{1,2}},'display','off');
disp(['Anova Test: EyeDataDegY.Static, p= ' num2str(p_YData_static)]);
disp(['Anova Test: EyeDataDegY.CP, p= ' num2str(p_YData_CP)])

[p_MSData_static,tbl_MSData_static,stats_MSData_static] = anovan(MSData_CombinedAcrossContrasts{1,1},{contrastFactor{1,1}},'display','off');
[p_MSData_CP,tbl_MSData_CP,stats_MSData_CP] = anovan(MSData_CombinedAcrossContrasts{1,2},{contrastFactor{1,2}},'display','off');
disp(['Anova Test: EyeMS.Static, p= ' num2str(p_MSData_static)]);
disp(['Anova Test: EyeMS.CP, p= ' num2str(p_MSData_CP)])

tickLengthPlot = 4*get(hPlotsFig1.hPlot1(1),'TickLength');
fontSize = 12;
for i = 1:length(cValsUnique)
    set(hPlotsFig1.hPlot1(:,i),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlotsFig1.hPlot2(:,i),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end

for i=1:3
    set(hPlotsFig1.hPlot3(i,1),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
    set(hPlotsFig1.hPlot4(i,1),'fontSize',fontSize,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end



%% Rescaling plots
rescaleData(hPlotsFig1.hPlot1(1,1:5),microsaccadeData.timeRange(1),microsaccadeData.timeRange(2),[-0.5 0.5],12);
yTicks = [-0.3 0 0.3];set(hPlotsFig1.hPlot1(1,1),'yTick',yTicks,'yTicklabel',yTicks);
rescaleData(hPlotsFig1.hPlot1(2,1:5),microsaccadeData.timeRange(1),microsaccadeData.timeRange(2),[-0.5 0.5],12);
set(hPlotsFig1.hPlot1(2,1),'yTick',[-0.3 0 0.3],'yTicklabel',[-0.3 0 0.3]);
rescaleData(hPlotsFig1.hPlot1(3,1:5),microsaccadeData.timeRange(1),microsaccadeData.timeRange(2),[0,3],12);
set(hPlotsFig1.hPlot1(3,1:2:5),'xTicklabel',[-0.5 0 0.5]);
set(hPlotsFig1.hPlot1(3,1),'yTick',[1 2],'yTicklabel',[1 2]);

rescaleData(hPlotsFig1.hPlot2(1,1:5),-0.5,0.5,[-0.5 0.5],12);
set(hPlotsFig1.hPlot2(1,1),'yTick',[-0.3 0 0.3],'yTicklabel',[-0.3 0 0.3]);
rescaleData(hPlotsFig1.hPlot2(2,1:5),-0.5,0.5,[-0.5 0.5],12);
set(hPlotsFig1.hPlot2(2,1),'yTick',[-0.3 0 0.3],'yTicklabel',[-0.3 0 0.3]);
rescaleData(hPlotsFig1.hPlot2(3,1:5),-0.5,0.5,[0,3],12);
set(hPlotsFig1.hPlot2(3,1:2:5),'xTicklabel',[-0.5 0 0.5]);
set(hPlotsFig1.hPlot2(3,1),'yTick',[1 2],'yTicklabel',[1 2]);




%% Labels
ylabel(hPlotsFig1.hPlot1(1,1),{'Eye-Position','Horizontal','(deg)'})
ylabel(hPlotsFig1.hPlot1(2,1),{'Eye-Position','Vertical','(deg)'})
ylabel(hPlotsFig1.hPlot1(3,1),{'Microsaccade' 'Rate', '(/s)'})

ylabel(hPlotsFig1.hPlot2(1,1),{'Eye-Position','Horizontal','(deg)'})
ylabel(hPlotsFig1.hPlot2(2,1),{'Eye-Position','Vertical','(deg)'})
ylabel(hPlotsFig1.hPlot2(3,1),{'Microsaccade' 'Rate' '(/s)'})

xlabel(hPlotsFig1.hPlot2(3,1),'Time(s)')
xlabel(hPlotsFig1.hPlot4(3,1),'Contrast Condition')

set(hPlotsFig1.hPlot3(1,1),'xTickLabel',[]);
set(hPlotsFig1.hPlot3(2,1),'xTickLabel',[]);
set(hPlotsFig1.hPlot4(1,1),'xTickLabel',[]);
set(hPlotsFig1.hPlot4(2,1),'xTickLabel',[]);



for i=1:3
    xlim(hPlotsFig1.hPlot3(i),[0 26]);
    xlim(hPlotsFig1.hPlot4(i),[0 26]);
    ylim(hPlotsFig1.hPlot3(i),[-0.1 0.1]);
        ylim(hPlotsFig1.hPlot4(i),[-0.1 0.1]);

    set(hPlotsFig1.hPlot3(i),'YTick',[-0.05 0 0.05]);
    set(hPlotsFig1.hPlot4(i),'YTick',[-0.05 0 0.05]);

    if i==3
    ylim(hPlotsFig1.hPlot3(i),[0 3]);
    set(hPlotsFig1.hPlot3(i),'YTick',[1 2]);
        ylim(hPlotsFig1.hPlot4(i),[0 3]);
    set(hPlotsFig1.hPlot4(i),'YTick',[1 2]);
    end
end
textH{1} = getPlotHandles(1,1,[0.05 0.97 0.01 0.01]);
textH{2} = getPlotHandles(1,1,[0.05 0.5 0.01 0.01]);
textString = {'A','B'};
for i = 1:2
    set(textH{i},'Visible','Off')
    text(0.35,1.15,textString{i},'unit','normalized','fontsize',20,'fontweight','bold','parent',textH{i});
end
folderSaveFigs = 'E:\Projects\Aritra_PlaidNormalizationProject\Figures\Color\';
FigName = fullfile(folderSaveFigs,[monkeyName 'eyeData']);
saveas(hFigure1,[FigName '.fig'])
print(hFigure1,[FigName,'.tif'],'-dtiff','-cmyk','-r300')

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Accessory Functions  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get parameter combinations
function [parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,...
    cValsUnique,tValsUnique,aValsUnique2,eValsUnique2,sValsUnique2,...
    fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract) %#ok<*STOUT>

load(fullfile(folderExtract,'parameterCombinations.mat')); %#ok<*LOAD>

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;            end

end

% Get Bad Trials
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
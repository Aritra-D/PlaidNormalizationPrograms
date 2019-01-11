function displayNormalizationData(folderSourceString)

if ~exist('folderSourceString','var');      folderSourceString = 'E:\'; end

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.72; backgroundColor = 'w'; % Panels

hFigure = figure(1);
set(hFigure,'units','normalized','outerposition',[0 0 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Session(s) panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% FileNameString
[fileNameStringAll,fileNameStringListArray] = getFileNameStringList;

hSessionPanel= uipanel('Title','Session(s)','titleposition','centertop',...
                'fontSize',fontSizeLarge,'Unit','Normalized','Position',...
                [0.05 0.925 0.9 0.08]);
            
hSession = uicontrol('Parent',hSessionPanel,'Unit','Normalized', ...
                    'BackgroundColor', backgroundColor, 'Position', ...
                    [0.4 0.83 0.2 0.16], 'Style','popup','String',...
                    fileNameStringAll,'FontSize',fontSizeLarge);

uicontrol('Parent',hSessionPanel,'Unit','Normalized',...
    'Position',[0.7 0.1 0.2 0.8],'Style','pushbutton','String','Select',...
    'FontSize',fontSizeMedium,'Callback',{@selectSession_Callback});

    function selectSession_Callback(~,~)
        session = get(hSession,'val'); 
        fileNameStringTMP = fileNameStringListArray{session};
        if strcmp(fileNameStringTMP{1}(1:5),'alpaH')
            monkeyName = 'alpaH'; 
        elseif strcmp(fileNameStringTMP{1}(1:7),'kesariH')
            monkeyName = 'kesariH';
        end
        gridType = 'Microelectrode';

        % Show electrodes on Grid
        electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
        hElectrodesonGrid = showElectrodeLocations(electrodeGridPos,[], ...
        [],[],1,0,gridType,monkeyName); %#ok<NASGU>

        %%%%%%%%%%%%%%%%%%%%%%%%%% Find Good Electrodes %%%%%%%%%%%%%%%%%%%
        [ElectrodeStringListAll,ElectrodeArrayListAll]= ...
            getElectrodesList(folderSourceString);
        ElectrodeListSession = ElectrodeArrayListAll{session};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters panel %%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        hParameterPanel = uipanel('Title','Parameters','fontSize', ...
            fontSizeLarge,'Unit','Normalized','Position',...
            [0.25 panelStartHeight 0.25 panelHeight]);
        paramsHeight=1/6;
        
        % ElectrodeString
        uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
            'Position',[0 1-paramsHeight 0.5 paramsHeight],...
            'Style','text','String','Electrode','FontSize',fontSizeMedium);
        hElectrode = uicontrol('Parent',hParameterPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor, ...
            'Position',[0.5 1-paramsHeight 0.5 paramsHeight],...
            'Style','popup','String',ElectrodeStringListAll{session},...
            'FontSize',fontSizeMedium);
        
        % Analysis Method
        analysisMethodString ='FFT Amplitude|Multi-Taper Power';

        uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
            'Position',[0 1-2*paramsHeight 0.5 paramsHeight],...
            'Style','text','String','AnalysisMethod',...
            'FontSize',fontSizeMedium);
        hAnalysisMethod = uicontrol('Parent',hParameterPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[0.5 1-2*paramsHeight 0.5 paramsHeight], ...
            'Style','popup','String',analysisMethodString,...
            'FontSize',fontSizeMedium);          
        
        % Analysis Type
        analysisTypeString = ...
        ['ERP|Firing Rate|Raster|Alpha [8-12 Hz]|Gamma Power [30-60 Hz]|'...
         'SSVEP (16 Hz)|STA'];

        uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
            'Position',[0 1-3*paramsHeight 0.5 paramsHeight],...
            'Style','text','String','AnalysisMeasure',...
            'FontSize',fontSizeMedium);
        hAnalysisType = uicontrol('Parent',hParameterPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[0.5 1-3*paramsHeight 0.5 paramsHeight], ...
            'Style','popup','String',analysisTypeString,...
            'FontSize',fontSizeMedium);  
        
        hAbsoluteMeasures = uicontrol('Parent',hParameterPanel,...
            'Unit','Normalized', ...
            'Position',[0 1-5*paramsHeight 1 paramsHeight], ...
            'Style','togglebutton',...
            'String','Show Absolute Measures',...
            'FontSize',fontSizeMedium);
        
        hNormalizeData = uicontrol('Parent',hParameterPanel,...
            'Unit','Normalized', ...
            'Position',[0 1-6*paramsHeight 1 paramsHeight], ...
            'Style','togglebutton',...
            'String','Normalize Data across sessions',...
            'FontSize',fontSizeMedium);              
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        timingTextWidth = 0.5; timingBoxWidth = 0.25;
        hTimingPanel = uipanel('Title','X and Y Limits',...
            'fontSize', fontSizeLarge,'Unit','Normalized',...
            'Position',[0.5 panelStartHeight 0.25 panelHeight]);
        timingHeight = 1/6; 

        signalRange = [-0.1 0.5]; 
        erpPeriod = [0.05 0.2]; 
        fftRange = [0 100];

        % Signal Range
        uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'Position',[0 1-timingHeight timingTextWidth timingHeight],...
            'Style','text','String','Parameter','FontSize',fontSizeMedium);

        uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'Position',[timingTextWidth 1-timingHeight timingBoxWidth...
            timingHeight], 'Style','text','String','Min',...
            'FontSize',fontSizeMedium);

        uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'Position',[timingTextWidth+timingBoxWidth 1-timingHeight ...
            timingBoxWidth timingHeight], ...
            'Style','text','String','Max','FontSize',fontSizeMedium);

        % Stim Range
        uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'Position',[0 1-2*timingHeight timingTextWidth timingHeight]...
            ,'Style','text','String','Stim Range (s)',...
            'FontSize',fontSizeSmall);
        hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'BackgroundColor', backgroundColor, ...
            'Position',[timingTextWidth 1-2*timingHeight timingBoxWidth...
            timingHeight], ...
            'Style','edit','String',num2str(signalRange(1)),...
            'FontSize',fontSizeSmall);
        hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
            'BackgroundColor', backgroundColor, ...
            'Position',[timingTextWidth+timingBoxWidth 1-2*timingHeight...
            timingBoxWidth timingHeight],'Style','edit',...
            'String',num2str(signalRange(2)),'FontSize',fontSizeSmall);
        
        % ERP Range
        uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'Position',[0 1-3*timingHeight timingTextWidth timingHeight]...
            ,'Style','text','String','ERP Range (s)',...
            'FontSize',fontSizeSmall);
        hERPMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'BackgroundColor', backgroundColor, ...
            'Position',[timingTextWidth 1-3*timingHeight ...
            timingBoxWidth timingHeight], ...
            'Style','edit','String',num2str(erpPeriod(1)),...
            'FontSize',fontSizeSmall);
        hERPMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'BackgroundColor', backgroundColor, ...
            'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight...
            timingBoxWidth timingHeight], ...
            'Style','edit','String',num2str(erpPeriod(2)),...
            'FontSize',fontSizeSmall);

        % FFT Range
        uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
            'Position',[0 1-4*timingHeight timingTextWidth timingHeight]...
            ,'Style','text','String','FFT Range (Hz)',...
            'FontSize',fontSizeSmall);
        hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
            'BackgroundColor', backgroundColor, ...
            'Position',[timingTextWidth 1-4*timingHeight...
            timingBoxWidth timingHeight], ...
            'Style','edit','String',num2str(fftRange(1)),...
            'FontSize',fontSizeSmall);
        hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
            'BackgroundColor', backgroundColor, ...
            'Position',[timingTextWidth+timingBoxWidth 1-4*timingHeight...
            timingBoxWidth timingHeight],'Style','edit',...
            'String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

        % Baseline Range
        baseline = [-0.25 0];
        uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
            'Position',[0 1-5*timingHeight timingTextWidth timingHeight]...
            ,'Style','text','String','Baseline (s)',...
            'FontSize',fontSizeSmall);
        hBaselineMin = uicontrol('Parent',hTimingPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[timingTextWidth 1-5*timingHeight...
            timingBoxWidth timingHeight],'Style','edit',...
            'String',num2str(baseline(1)),'FontSize',fontSizeSmall);
        hBaselineMax = uicontrol('Parent',hTimingPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight...
            timingBoxWidth timingHeight], ...
            'Style','edit','String',num2str(baseline(2)),...
            'FontSize',fontSizeSmall);

        % Stimulus Range
        stimPeriod = [0.15 0.4];
        uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
            'Position',[0 1-6*timingHeight timingTextWidth timingHeight]...
            ,'Style','text','String','Stim Period (s)',...
            'FontSize',fontSizeSmall);
        hStimPeriodMin = uicontrol('Parent',hTimingPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[timingTextWidth 1-6*timingHeight...
            timingBoxWidth timingHeight],'Style','edit',...
            'String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
        hStimPeriodMax = uicontrol('Parent',hTimingPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight...
            timingBoxWidth timingHeight],'Style','edit',...
            'String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hPlotOptionsPanel = uipanel('Title','Plotting Options',...
            'fontSize', fontSizeLarge, 'Unit','Normalized',...
            'Position',[0.75 panelStartHeight 0.2 panelHeight]);
        plotOptionsHeight = 1/6;

        % Button for Plotting
        [colorString, colorNames] = getColorString;
        uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
            'Position',[0 5*plotOptionsHeight 0.6 plotOptionsHeight],...
            'Style','text','String','Color','FontSize',fontSizeSmall);
        hChooseColor = uicontrol('Parent',hPlotOptionsPanel,...
            'Unit','Normalized','BackgroundColor', backgroundColor,...
            'Position',[0.6 5*plotOptionsHeight 0.4 plotOptionsHeight],...
            'Style','popup','String',colorString,'FontSize',fontSizeSmall);

        uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
            'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight],...
            'Style','pushbutton','String','cla',...
            'FontSize',fontSizeMedium,'Callback',{@cla_Callback});

        hHoldOn = uicontrol('Parent',hPlotOptionsPanel,...
            'Unit','Normalized','Position',[0 2*plotOptionsHeight ...
            1 plotOptionsHeight],'Style','togglebutton',...
            'String','hold on','FontSize',fontSizeMedium, ...
            'Callback',{@holdOn_Callback});

        uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
            'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
            'Style','pushbutton','String','Rescale',...
            'FontSize',fontSizeMedium,'Callback',{@rescaleXY_Callback});

        uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
            'Position',[0 0 1 plotOptionsHeight],...
            'Style','pushbutton','String','plot',...
            'FontSize',fontSizeMedium,'Callback',{@plotData_Callback});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots
        % Main Plot showing Neural mesaures in a 5 x 5 contrast conditions
        staticStartPos = 0.05;
        startXPos = staticStartPos; startYPos = 0.05; mainTabHeight = 0.55; 
        mainTabWidth = 0.4;     
        
        cValsUnique = [0 12.5 25 50 100];
        cValsUnique2 = [0 12.5 25 50 100];
        numRows = length(cValsUnique); numCols = length(cValsUnique2);
        gridPos=[0.02+startXPos startYPos mainTabWidth mainTabHeight];
        gap = 0.002;
        plotHandles = getPlotHandles(numRows,numCols,gridPos,gap);
        
        % RF Positions with Stimulus Centre
        hRFCentresandStim = ...
            getPlotHandles(1,1,[0.52 0.35 0.14 0.25],0.05,0.05,1);
        
        % Color Matrix of Neural measures in a 5x5 contrast conditions
        hNeuralMeasureColorMatrix = ...
            getPlotHandles(1,1,[0.52 0.05 0.14 0.25],0.001,0.001,1);
        
        plotHandles2= getPlotHandles(1,5,[0.7 0.5 0.25 0.1],0.001,0.001,1);
        plotHandles3= getPlotHandles(1,5,[0.7 0.2 0.25 0.1],0.001,0.001,1);
        
        hRowCRF = getPlotHandles(1,1,[0.7 0.35 0.1 0.1],0.001,0.001,1);
        hColumnCRF = getPlotHandles(1,1,[0.7 0.05 0.1 0.1],0.001,0.001,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        freqRanges{1} = [8 12]; % alpha
        freqRanges{2} = [30 60]; % gamma
        freqRanges{3} = [16 16];         % SSVEP
        
%         freqRangeStr = {'alpha','gamma','SSVEP'};
%         numFreqRanges = length(freqRanges);        
        
        % Plotting Functions
        function plotData_Callback(~,~)
            
            %%%%%%%%%%%%%%%%%%%%%%%% Read values %%%%%%%%%%%%%%%%%%%%%%%%%%
            electrodeString = get(hElectrode,'val'); 
            if length(fileNameStringTMP) ==1
                ElectrodeListTMP = ElectrodeListSession(electrodeString);
            elseif length(fileNameStringTMP)>1
                ElectrodeListTMP = ElectrodeListSession;
            end
            if isempty(ElectrodeListTMP{1})
                error(['No electrode found for analysis!'...
                        'Please try another session!'])
            end
            analysisMeasure = get(hAnalysisType,'val');
            NormalizeDataFlag = get(hNormalizeData,'val');
            AbsoluteMeasuresFlag = get(hAbsoluteMeasures,'val');
            
            erpRange = [str2double(get(hERPMin,'String'))...
                        str2double(get(hERPMax,'String'))];
            blRange = [str2double(get(hBaselineMin,'String'))...
                       str2double(get(hBaselineMax,'String'))];
            stRange = [str2double(get(hStimPeriodMin,'String'))...
                       str2double(get(hStimPeriodMax,'String'))];

            plotColor = colorNames(get(hChooseColor,'val'));
            holdOnState = get(hHoldOn,'val'); %#ok<NASGU>
            
            [erpData,firingRateData,fftData,energyData,~] = getData(folderSourceString,...
             fileNameStringTMP,ElectrodeListTMP,erpRange,blRange,...
             stRange,freqRanges); 

%             PlotConstant = [4 2 0 -2 -4];
            if analysisMeasure == 1 % computing ERP
                plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,erpData.timeVals,erpData,plotColor,analysisMeasure,AbsoluteMeasuresFlag)
            elseif analysisMeasure == 2 % computing Firing rate
                plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,firingRateData.timeVals,firingRateData,plotColor,analysisMeasure,AbsoluteMeasuresFlag)
            elseif analysisMeasure == 3 % computing Raster Plot from spike data
               error('Still working on raster data!')
            elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6 % computing alpha
                plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,fftData.freqVals,fftData,plotColor,analysisMeasure,AbsoluteMeasuresFlag)
            elseif analysisType == 7 % need to work on STA!
                error('STA computation method not found') 

            
                
%                 if analysisType == 1 % compute ERP
%                     DataSize = size(erpData);
%                     if DataSize(1)==1
%                         erpData = squeeze(erpData(:,1,:,:,:));
%                     elseif DataSize(1)>1
%                         erpData = squeeze(mean(squeeze(erpData(:,1,:,:,:)),1));
%                     end
%                     for c1 = 1:5
%                         for c2=1:5
%                             plot(plotHandles(c1+PlotConstant(c1),c2),timeVals,...
%                                 squeeze(erpData(c1,c2,:)),'color',plotColor);
%                         end
%                     end
%                     %color Coded Matrix of firing Rate
%                     if DataSize(1)==1
%                         RMSvalsERP = squeeze(RMSvalsERP(:,1,:,:));
%                     elseif DataSize(1)>1
%                         RMSvalsERP = squeeze(mean(squeeze(RMSvalsERP(:,1,:,:)),1));
%                     end
%                     RMSvalsFlipped = ...
%                         [RMSvalsERP(5,:);RMSvalsERP(4,:);...
%                         RMSvalsERP(3,:);RMSvalsERP(2,:);RMSvalsERP(1,:);];
%                     imagesc(RMSvalsFlipped,'parent',hNeuralMeasureColorMatrix);
%                     colorbar(hNeuralMeasureColorMatrix);
%                     
%                     % CRF Row-wise & Column-wise
%                      CRFColors = jet(length(cValsUnique));
%                         for iCon = 1:5
%                             plot(plotHandles2(1,iCon),cValsUnique,RMSvalsERP(iCon,:),...
%                                 'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                             plot(plotHandles3(1,iCon),cValsUnique2,RMSvalsERP(:,iCon),...
%                                 'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                             plot(hRowCRF(1,1),cValsUnique,RMSvalsERP(iCon,:),...
%                                 'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                             hold(hRowCRF(1,1),'on')
%                             plot(hColumnCRF(1,1),cValsUnique2,RMSvalsERP(:,iCon),...
%                                 'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                             hold(hColumnCRF(1,1),'on')
%                         end
%                      hold(hRowCRF(1,1),'off'); hold(hColumnCRF(1,1),'off')
                     
%                 elseif analysisType == 4 || analysisType == 5
%                         if size(fftDataBL)==size(fftDataST)
%                             DataSize = size(fftDataST);
%                         else
%                             error('Size of fftDataBL and fftDataST do not match!')
%                         end
%                         
%                         if DataSize(1)==1
%                         fftDataBL = squeeze(fftDataBL(:,1,:,:,:));
%                         fftDataST = squeeze(fftDataST(:,1,:,:,:));
%                         elseif DataSize(1)>1
%                         fftDataBL = squeeze(mean(squeeze(fftDataBL(:,1,:,:,:)),1));
%                         fftDataST = squeeze(mean(squeeze(fftDataST(:,1,:,:,:)),1));
%                         end
%                         
%                         for c1 = 1:5
%                             for c2=1:5
%                                 plot(plotHandles(c1+PlotConstant(c1),c2),freqVals,...
%                                     squeeze(fftDataBL(c1,c2,:)),'g');
%                                 hold(plotHandles(c1+PlotConstant(c1),c2),'on')
%                                 plot(plotHandles(c1+PlotConstant(c1),c2),freqVals,...
%                                     squeeze(fftDataST(c1,c2,:)),'k');
%                                 hold(plotHandles(c1+PlotConstant(c1),c2),'off')
%                             end
%                         end
%                         
%                         if analysisType == 4 
%                         %color Coded Matrix of firing Rate
%                         if DataSize(1)==1
%                             alphaData = squeeze(mean(alphaData(:,1,:,:,:),5));
%                         elseif DataSize(1)>1
%                             alphaData =...
%                                 squeeze(mean(squeeze(mean(alphaData(:,1,:,:,:),5)),1));
%                         end
%                         alphaDataFlipped = ...
%                             [alphaData(5,:);alphaData(4,:);...
%                             alphaData(3,:);alphaData(2,:);alphaData(1,:);];
%                         imagesc(alphaDataFlipped,'parent',hNeuralMeasureColorMatrix);
%                         colorbar(hNeuralMeasureColorMatrix);
%                     
%                         % CRF Row-wise & Column-wise
%                          CRFColors = jet(length(cValsUnique));
%                             for iCon = 1:5
%                                 plot(plotHandles2(1,iCon),cValsUnique,alphaData(iCon,:),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 plot(plotHandles3(1,iCon),cValsUnique2,alphaData(:,iCon),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 plot(hRowCRF(1,1),cValsUnique,alphaData(iCon,:),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 hold(hRowCRF(1,1),'on')
%                                 plot(hColumnCRF(1,1),cValsUnique2,alphaData(:,iCon),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 hold(hColumnCRF(1,1),'on')
%                             end
%                          hold(hRowCRF(1,1),'off'); hold(hColumnCRF(1,1),'off')
%                          
%                         elseif analysisType == 5
%                         %color Coded Matrix of firing Rate
%                         if DataSize(1)==1
%                             gammaData = squeeze(mean(gammaData(:,1,:,:,:),5));
%                         elseif DataSize(1)>1
%                             gammaData = squeeze(mean(squeeze(mean(gammaData(:,1,:,:,:),5)),1));
%                         end
%                         gammaDataFlipped = ...
%                             [gammaData(5,:);gammaData(4,:);...
%                             gammaData(3,:);gammaData(2,:);gammaData(1,:);];
%                         imagesc(gammaDataFlipped,'parent',hNeuralMeasureColorMatrix);
%                         colorbar(hNeuralMeasureColorMatrix);
%                     
%                         % CRF Row-wise & Column-wise
%                          CRFColors = jet(length(cValsUnique));
%                             for iCon = 1:5
%                                 plot(plotHandles2(1,iCon),cValsUnique,gammaData(iCon,:),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 plot(plotHandles3(1,iCon),cValsUnique2,gammaData(:,iCon),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 plot(hRowCRF(1,1),cValsUnique,gammaData(iCon,:),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 hold(hRowCRF(1,1),'on')
%                                 plot(hColumnCRF(1,1),cValsUnique2,gammaData(:,iCon),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 hold(hColumnCRF(1,1),'on')
%                             end
%                          hold(hRowCRF(1,1),'off'); hold(hColumnCRF(1,1),'off')                            
%                             
%                         end
%                         
%                    elseif analysisType == 8
%                         if size(fftDataBL)==size(fftDataST)
%                             DataSize = size(fftDataST);
%                         else
%                             error('Size of fftDataBL and fftDataST do not match!')
%                         end
%                         
%                         if DataSize(1)==1
%                         fftDataBL = squeeze(fftDataBL(:,2,:,:,:));
%                         fftDataST = squeeze(fftDataST(:,2,:,:,:));
%                         elseif DataSize(1)>1
%                         fftDataBL = squeeze(mean(squeeze(fftDataBL(:,2,:,:,:)),1));
%                         fftDataST = squeeze(mean(squeeze(fftDataST(:,2,:,:,:)),1));
%                         end
%                         
%                         for c1 = 1:5
%                             for c2=1:5
%                                 plot(plotHandles(c1+PlotConstant(c1),c2),...
%                                     freqVals,squeeze(fftDataBL(c1,c2,:)),'g');
%                                 hold(plotHandles(c1+PlotConstant(c1),c2),'on')
%                                 plot(plotHandles(c1+PlotConstant(c1),c2),...
%                                     freqVals,squeeze(fftDataST(c1,c2,:)),'k');
%                                 hold(plotHandles(c1+PlotConstant(c1),c2),'off')
%                             end
%                         end
%                         
%                         
%                         %color Coded Matrix of firing Rate
%                         if DataSize(1)==1
%                             ssvepData = squeeze(mean(ssvepData(:,2,:,:,:),5));
%                         elseif DataSize(1)>1
%                             ssvepData = squeeze(mean(squeeze(mean(ssvepData(:,2,:,:,:),5)),1));
%                         end
%                         ssvepDataFlipped = ...
%                             [ssvepData(5,:);ssvepData(4,:);...
%                             ssvepData(3,:);ssvepData(2,:);ssvepData(1,:);];
%                         imagesc(ssvepDataFlipped,'parent',hNeuralMeasureColorMatrix);
%                         colorbar(hNeuralMeasureColorMatrix);
%                     
%                         % CRF Row-wise & Column-wise
%                          CRFColors = jet(length(cValsUnique));
%                             for iCon = 1:5
%                                 plot(plotHandles2(1,iCon),cValsUnique,ssvepData(iCon,:),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 plot(plotHandles3(1,iCon),cValsUnique2,ssvepData(:,iCon),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 plot(hRowCRF(1,1),cValsUnique,ssvepData(iCon,:),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 hold(hRowCRF(1,1),'on')
%                                 plot(hColumnCRF(1,1),cValsUnique2,ssvepData(:,iCon),...
%                                     'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
%                                 hold(hColumnCRF(1,1),'on')
%                             end
%                          hold(hRowCRF(1,1),'off'); hold(hColumnCRF(1,1),'off')
%                         
%                 end
            end
            
            if analysisMeasure<=3 %ERP or spikes
                xMin = str2double(get(hStimMin,'String'));
                xMax = str2double(get(hStimMax,'String'));
            elseif analysisType <=6  % LFP fft analysis
                xMin = str2double(get(hFFTMin,'String'));
                xMax = str2double(get(hFFTMax,'String'));                
            else
                xMin = str2double(get(hSTAMin,'String'));
                xMax = str2double(get(hSTAMax,'String'));
            end
            
            YLabel = get(hNeuralMeasureColorMatrix,'YTickLabel');
            set(hNeuralMeasureColorMatrix,'YTickLabel',flipud(YLabel));
            rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
            rescaleData(plotHandles2,0,100,getYLims(plotHandles2));
            rescaleData(plotHandles3,0,100,getYLims(plotHandles3));
            rescaleData(hRowCRF,0,100,getYLims(hRowCRF));
            rescaleData(hColumnCRF,0,100,getYLims(hColumnCRF));
        end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleXY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType==6
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        else    
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(plotHandles2,0,100,getYLims(plotHandles2));
        rescaleData(plotHandles3,0,100,getYLims(plotHandles3));

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
        
        holdOnGivenPlotHandle(plotHandles,holdOnState);
        holdOnGivenPlotHandle(plotHandles2,holdOnState);
        holdOnGivenPlotHandle(plotHandles3,holdOnState);
        
        if holdOnState
            set(hElectrodes,'Nextplot','add');
        else
            set(hElectrodes,'Nextplot','replace');
        end

        function holdOnGivenPlotHandle(plotHandles,holdOnState)
            
            [numRows,numCols] = size(plotHandles);
            if holdOnState
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','add');

                    end
                end
            else
                for i=1:numRows
                    for j=1:numCols
                        set(plotHandles(i,j),'Nextplot','replace');
                    end
                end
            end
        end 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        
        claGivenPlotHandle(plotHandles);
        claGivenPlotHandle(plotHandles2);
        claGivenPlotHandle(plotHandles3);
        claGivenPlotHandle(hRowCRF);
        claGivenPlotHandle(hColumnCRF);
        claGivenPlotHandle(hNeuralMeasureColorMatrix);
        claGivenPlotHandle(hRFCentresandStim)
       
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[erpData,firingRateData,fftData,energyData,electrodeArray] = ...
    getData(folderSourceString,fileNameStringTMP,ElectrodeListTMP,...
    erpRange,blRange,stRange,freqRanges)

numDatasets = length(fileNameStringTMP);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[erpData,firingRateData,fftData,energyData,electrodeArray]...
= getDataSingleSession(folderSourceString,fileNameStringTMP,...
ElectrodeListTMP{1},erpRange,blRange,stRange,freqRanges); 

if length(fileNameStringTMP)>1
    for i=2:numDatasets
        if isempty(ElectrodeListTMP{i})
           continue
        end
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [erpDataTMP,firingRateDataTMP,fftDataTMP,energyDataTMP,~] = getDataSingleSession(folderSourceString,fileNameStringTMP,...
            ElectrodeListTMP{i},erpRange,blRange,stRange,freqRanges);
        
        erpData.data = cat(1,erpData.data,erpData.dataTMP);
        erpData.analysisData = cat(1,erpData.analysisData,erpDataTMP.analysisData);
        
        firingRateData.data = cat(1,firingRateData.data,firingRateDataTMP.data);
        firingRateData.analysisDataBL = cat(1,firingRateData.analysisDataBL,firingRateDataTMP.analysisDataBL);
        firingRateData.analysisDataST = cat(1,firingRateData.analysisDataST,firingRateDataTMP.analysisDataST);
        
        fftData.dataBL = cat(1,fftData.dataBL,fftDataTMP.dataBL);
        fftData.dataST = cat(1,fftData.dataST,fftDataTMP.dataST);
        fftData.analysisdataBL = cat(1,fftData.analysisdataBL,fftDataTMP.analysisdataBL);
        fftData.analysisdataST = cat(1,fftData.analysisdataST,fftDataTMP.analysisdataST);
        
        energyData.dataBL = cat(1,energyData.dataBL,energyDataTMP.dataBL);
        energyData.dataST = cat(1,energyData.dataST,energyDataTMP.dataST);
        energyData.analysisdataBL = cat(1,energyData.analysisdataBL,energyDataTMP.analysisdataBL);
        energyData.analysisdataST = cat(1,energyData.analysisdataST,energyDataTMP.analysisdataST);
        
        electrodeArray = [];
%         psthData = cat(1,psthData,psthDataTMP);
%         firingRates = cat(1,firingRates,firingRatesTMP);
%         erpData = cat(1,erpData,erpDataTMP);
%         fftDataBL = cat(1,fftDataBL,fftDataBLTMP);
%         fftDataST = cat(1,fftDataST,fftDataSTTMP);
%         alphaData = cat(1,alphaData,alphaDataTMP);
%         gammaData = cat(1,gammaData,gammaDataTMP);
%         ssvepData = cat(1,ssvepData,ssvepDataTMP);
%         RMSvalsERP = cat(1,RMSvalsERP,RMSvalsERPTMP);
        
%         disp(size(firingRates));
    end
end
end



function [erpData,firingRateData,fftData,energyData,electrodeArray] = ...
    getDataSingleSession(folderSourceString,fileNameStringTMP,...
    ElectrodeListTMP,erpRange,blRange,stRange,freqRanges)

if strcmp(fileNameStringTMP{1}(1:5),'alpaH')       
    monkeyName = 'alpaH'; 
    expDate = fileNameStringTMP{1}(6:11); 
    protocolName = fileNameStringTMP{1}(12:end); 
elseif strcmp(fileNameStringTMP{1}(1:7),'kesariH')
    monkeyName = 'kesariH';
    expDate = fileNameStringTMP{1}(8:13); 
    protocolName = fileNameStringTMP{1}(14:end);
end
gridType = 'microelectrode';
folderName = fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');

% Get Combinations
[parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,~,~,~,cValsUnique,tValsUnique, ...
    aValsUnique2,eValsUnique2,~,~,~,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract);

if aValsUnique ~= aValsUnique2 || eValsUnique ~= eValsUnique2
    error('Azimuths and/or elevations do not match!');
end                                          
a=1; e=1; s=1; f=1; o=1; 
if tValsUnique ~= tValsUnique2
    error('Azimuths and/or elevations do not match!');
else
    tList = 1:length(tValsUnique);
end

c1List = 1:length(cValsUnique);
c2List = 1:length(cValsUnique2);


% TimeVals info
[~,timeVals,~,~] = loadlfpInfo(folderLFP);
if iscell(ElectrodeListTMP)
ElectrodeList = cell2mat(ElectrodeListTMP);
else
    ElectrodeList = ElectrodeListTMP;
end

electrodeArray = ElectrodeList;
Fs = round(1/(timeVals(2)-timeVals(1)));
range = blRange;
rangePos = round(diff(range)*Fs);
erpRangePos = round(diff(erpRange)*Fs);
blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
erpPos = find(timeVals>=erpRange(1),1)+ (1:erpRangePos);
freqVals = 0:1/diff(range):Fs-1/diff(range);
numFreqs = length(freqRanges);
% alphaRangeHz = [8 12];
% gammaRangeHz = [30 60];
% ssvepFreqHz = 16;
% alphaPos = intersect(find(freqVals>=alphaRangeHz(1)),find(freqVals<=alphaRangeHz(2)));
% gammaPos = intersect(find(freqVals>=gammaRangeHz(1)),find(freqVals<=gammaRangeHz(2)));
% ssvepPos = find(freqVals==ssvepFreqHz);
unitID = 0;
for iElec = 1:length(ElectrodeList)
    % Get LFP data
    clear analogData
    load(fullfile(folderLFP,['elec' num2str(ElectrodeList(iElec)) '.mat']));
    % Get Spike data
    clear spikeData
    load(fullfile(folderSpikes,['elec' num2str(ElectrodeList(iElec)) '_SID' num2str(unitID) '.mat']));

%     Get bad trials
    badTrialFile = fullfile(folderSegment,'badTrials.mat');
    if ~exist(badTrialFile,'file')
        disp('Bad trial file does not exist...');
        badTrials=[];
    else
        badTrials = loadBadTrials(badTrialFile);
        disp([num2str(length(badTrials)) ' bad trials']);
    end
    % Set up MT
%     Fs              = round(1/(timeVals(2)-timeVals(1)));
    params.tapers   = [1 1];
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 100];
    params.trialave = 1;
    for t = 1:length(tList)
        for c1=1:length(c1List)
            for c2=1:length(c2List)

                clear goodPos fftBL fftST erp dataBL dataST
                goodPos = parameterCombinations{a,e,s,f,o,c1List(c1),tList(t)};
                goodPos2 = parameterCombinations2{a,e,s,f,o,c2List(c2),tList(t)};
                goodPos = intersect(goodPos,goodPos2);
                goodPos = setdiff(goodPos,badTrials);

                if isempty(goodPos)
                    disp('No entries for this combination..');
                else
                    disp(['pos=(' num2str(c1) ',' num2str(c2) ') ,n=' num2str(length(goodPos))]);
                    N(iElec,t,c1,c2) = length(goodPos);
                    
                    if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
                        disp('baseline and stimulus ranges are not the same');
                    else
                       [psthData(iElec,t,c1,c2,:),xsFR] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
                       erp = mean(analogData(goodPos,:),1); %#ok<NODEF>
                       erpDataTMP(iElec,t,c1,c2,:) = erp;
                       RMSvalsBL(iElec,t,c1,c2) = rms(erp(blPos));
                       RMSvalsERP(iElec,t,c1,c2) = rms(erp(erpPos));
                       
                       firingRatesST(iElec,t,c1,c2) = mean(getSpikeCounts(spikeData(goodPos),stRange))/diff(stRange);
                       firingRatesBL(iElec,t,c1,c2) = mean(getSpikeCounts(spikeData(goodPos),blRange))/diff(blRange);
                       
                       fftBL = log10(squeeze(mean(abs(fft(analogData(goodPos,blPos),[],2)))));
                       fftST = log10(squeeze(mean(abs(fft(analogData(goodPos,stPos),[],2)))));
                       
                       fftDataBL(iElec,t,c1,c2,:) = fftBL;
                       fftDataST(iElec,t,c1,c2,:) = fftST;
                       
                       %Power by MT method
                       dataBL = analogData(goodPos,blPos)';
                       [tmpEBL,freqValsBL] = mtspectrumc(dataBL,params);
                       dataST = analogData(goodPos,stPos)';
                       [tmpEST,freqValsST] = mtspectrumc(dataST,params);
                       
                       if isequal(freqValsBL,freqValsST)
                           freqValsMT = freqValsST;
                       end
                       
                       mEnergyVsFreqBL(iElec,t,c1,c2,:) = conv2Log(tmpEBL);
                       mEnergyVsFreqST(iElec,t,c1,c2,:) = conv2Log(tmpEST);
                       
                       for i=1:numFreqs
                           fftAmpST{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(fftST(:),freqVals,freqRanges{i}));
                           fftAmpBL{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(fftBL(:),freqVals,freqRanges{i}));
                           energyValsST{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),freqValsMT,freqRanges{i}));
                           energyValsBL{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),freqValsMT,freqRanges{i}));
                       end
                       
                       
%                        alphaData(iElec,t,c1,c2,:) = fftST(alphaPos); %alphaPowerBL(iElec,t,c1,c2,:) = fftBL(alphaPos);
%                        gammaData(iElec,t,c1,c2,:) = fftST(gammaPos); %gammaPowerBL(iElec,t,c1,c2,:) = fftBL(gammaPos);
% 
%                        ssvepData(iElec,t,c1,c2) = fftST(ssvepPos);    %#ok<FNDSB> %ssvepPowerBL(iElec,t,c1,c2) = fftBL(ssvepPos);
                    end
                end
            end
        end
    end
end

erpData.data = erpDataTMP;
erpData.analysisDataBL = RMSvalsBL;
erpData.analysisDataST = RMSvalsERP;
erpData.timeVals = timeVals;
erpData.N = N;

firingRateData.data = psthData;
firingRateData.analysisDataBL = firingRatesBL;
firingRateData.analysisDataST = firingRatesST;
firingRateData.timeVals = xsFR;
firingRateData.N = N;

fftData.dataBL = fftDataST;
fftData.dataST = fftDataBL;
fftData.analysisDataBL = fftAmpBL;
fftData.analysisDataST = fftAmpST;
fftData.freqVals = freqVals;
fftData.N = N;

energyData.dataBL=mEnergyVsFreqBL;
energyData.dataST=mEnergyVsFreqST;
energyData.analysisDataBL = energyValsBL;
energyData.analysisDataST = energyValsST;
energyData.freqVals = freqValsMT;
energyData.N = N;








end

% Accessory Functions
function [parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,...
    cValsUnique,tValsUnique,aValsUnique2,eValsUnique2,sValsUnique2,...
    fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;            end

end
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end

function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
function [fileNameStringAll,fileNameStringListArray] = getFileNameStringList

[tmpFileNameStringList,monkeyNameList] = getNormalizationExperimentDetails;

fileNameStringAll = ''; pos=1;
clear fileNameStringListArray

for i=1:length(monkeyNameList)
    for j=1:length(tmpFileNameStringList{i})
        fileNameStringAll = [cat(2,fileNameStringAll,tmpFileNameStringList{i}{j}) '|'];
        fileNameStringListArray{pos} = tmpFileNameStringList{i}(j); %#ok<*AGROW>
        pos=pos+1;
    end
end

allNames = [];
for i=1:length(monkeyNameList)
    fileNameStringAll = [cat(2,fileNameStringAll,monkeyNameList{i}) ...
        ' (N=' num2str(length(tmpFileNameStringList{i})) ')|'];
    fileNameStringListArray{pos} = tmpFileNameStringList{i};
    allNames = cat(2,allNames,tmpFileNameStringList{i});
    pos=pos+1;
end

fileNameStringAll = cat(2,fileNameStringAll,['all (N=' num2str(length(allNames)) ')']);
fileNameStringListArray{pos} = allNames;
end

function [ElectrodeStringListAll,ElectrodeArrayListAll] = getElectrodesList(folderSourceString)

[tmpElectrodeStringList,tmpElectrodeArrayList,allElecs,monkeyNameList] = getGoodElectrodesDetails(folderSourceString);

ElectrodeStringListAll = ''; 
ElectrodeArrayListAll = [];

for i=1:length(monkeyNameList)
        ElectrodeStringListAll = cat(2,ElectrodeStringListAll,tmpElectrodeStringList{i});
        ElectrodeArrayListAll = cat(2,ElectrodeArrayListAll,tmpElectrodeArrayList{i});
         
end

Sessions = length(ElectrodeArrayListAll);
pos = Sessions+1;
for i=1:length(monkeyNameList)
    clear j
    for j = 1:length(tmpElectrodeArrayList{i})
        ElectrodeArrayListAll{1,pos}{1,j} = tmpElectrodeArrayList{1,i}{1,j}{1,end};
    end
    ElectrodeStringListAll{1,pos} = ['all (N=' num2str(allElecs(i)) ')'];
    pos=pos+1;
end

AllSessions = length(ElectrodeArrayListAll)+1;

for k=1:length(tmpElectrodeArrayList{1})+length(tmpElectrodeArrayList{2})
    ElectrodeArrayListAll{1,AllSessions}{1,k} = ElectrodeArrayListAll{1,k}{1,end};
    ElectrodeStringListAll{1,AllSessions} = ['all (N=' num2str(sum(allElecs)) ')'];
end
    


end
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange)

posToAverage = intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2)));
eValue   = mean(mEnergy(posToAverage));
end
function normData = normalizeData(x)
for iElec = 1:size(x,1)
    for t = 1:size(x,2)
%         normData = x(iElec
    end
end
end
function plotData(hPlot1,hPlot2,hPlot3,hPlot4,hPlot5,hPlot6,xs,data,colorName,analysisMeasure,AbsoluteMeasuresFlag)

% Main 5x5 plot for Neural Measure
if analysisMeasure == 1 || analysisMeasure == 2
    dataSize = size(data.data);
    if dataSize(1) == 1
        dataPlot = squeeze(squeeze(data.data(:,1,:,:,:)));
    elseif dataSize(1) >1
        dataPlot = squeeze(mean(squeeze(data.data(:,1,:,:,:)),1));
    end
elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6
    if size(data.dataBL) == size(data.dataST) 
        dataSize = size(data.dataST);
    else
        error('Size of fftDataBL and fftDataST do not match!')
    end
    if dataSize(1) == 1
        dataPlotBL = squeeze(data.dataBL(:,1,:,:,:));
        dataPlotST = squeeze(data.dataST(:,1,:,:,:));
        if analysisMeasure == 6
            dataPlotBL = squeeze(data.dataBL(:,2,:,:,:));
            dataPlotST = squeeze(data.dataST(:,2,:,:,:));    
        end
    elseif dataSize(1) >1
        dataPlotBL = squeeze(mean(squeeze(data.dataBL(:,1,:,:,:)),1));
        dataPlotST = squeeze(mean(squeeze(data.dataST(:,1,:,:,:)),1));
        if analysisMeasure == 6
            dataPlotBL = squeeze(mean(squeeze(data.dataBL(:,1,:,:,:)),1));
            dataPlotST = squeeze(mean(squeeze(data.dataST(:,1,:,:,:)),1));
        end
    end
end
    
PlotConstant = [4 2 0 -2 -4];
for c1 = 1:5
    for c2 = 1:5
        if analysisMeasure == 1 || analysisMeasure == 2
            plot(hPlot1(c1+PlotConstant(c1),c2),xs,squeeze(dataPlot(c1,c2,:)),'color',colorName);
        elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6
            plot(hPlot1(c1+PlotConstant(c1),c2),xs,squeeze(dataPlotBL(c1,c2,:)),'g');
            hold(hPlot1(c1+PlotConstant(c1),c2),'on')
            plot(hPlot1(c1+PlotConstant(c1),c2),xs,squeeze(dataPlotST(c1,c2,:)),'k');
            hold(hPlot1(c1+PlotConstant(c1),c2),'off')
        end
    end
end

% Color coded Matrix of Neural Measure

if dataSize(1)==1
    if analysisMeasure == 1 || analysisMeasure == 2
        analysisData = squeeze(data.analysisDataST(:,1,:,:));
    elseif analysisMeasure == 4 
        analysisData = squeeze(data.analysisDataST{1}(:,1,:,:));
    elseif analysisMeasure == 5
        analysisData = squeeze(data.analysisDataST{2}(:,1,:,:));
    elseif analysisMeasure == 6
        analysisData = squeeze(data.analysisDataST{3}(:,2,:,:));
    end
elseif dataSize(1)>1
    if analysisMeasure == 1 || analysisMeasure == 2
        analysisData = squeeze(mean(squeeze(data.analysisDataST(:,1,:,:)),1));
    elseif analysisMeasure == 4
        analysisData = squeeze(mean(squeeze(data.analysisDataST{1}(:,1,:,:)),1));
    elseif analysisMeasure == 5
        analysisData = squeeze(mean(squeeze(data.analysisDataST{2}(:,1,:,:)),1));
    elseif analysisMeasure == 6
        analysisData = squeeze(mean(squeeze(data.analysisDataST{3}(:,2,:,:)),1));
    end
end
   
imagesc(flip(analysisData,1),'parent',hPlot2);colorbar(hPlot2);

% Contrast Response curves Row-Wise & Column-wise
cValsUnique = [0 12.5 25 50 100];
cValsUnique2 = [0 12.5 25 50 100];
CRFColors = jet(length(cValsUnique));
for iCon = 1:5
    plot(hPlot3(1,iCon),cValsUnique,analysisData(iCon,:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    plot(hPlot4(1,iCon),cValsUnique2,analysisData(:,iCon),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    plot(hPlot5(1,1),cValsUnique,analysisData(iCon,:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlot5(1,1),'on')
    plot(hPlot6(1,1),cValsUnique2,analysisData(:,iCon),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlot6(1,1),'on')
end
hold(hPlot5(1,1),'off'); hold(hPlot6(1,1),'off')
end
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored); %#ok<NODEF>
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

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
function rescaleData(plotHandles,xMin,xMax,yLims)

[numRows,numCols] = size(plotHandles);
labelSize=11;
for i=1:numRows
    for j=1:numCols
        axis(plotHandles(i,j),[xMin xMax yLims]);
        if (i==numRows && rem(j,2)==1)
            if j==1
                set(plotHandles(i,j),'fontSize',labelSize);
            elseif j~=1
                set(plotHandles(i,j),'YTickLabel',[],'fontSize',labelSize);
            end
        elseif (rem(i,2)==0 && j==1)
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);

        else 
            set(plotHandles(i,j),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
        end
    end
end

% Remove Labels on the four corners
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




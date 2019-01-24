function displayNormalizationData(folderSourceString)

if ~exist('folderSourceString','var');      folderSourceString = 'M:\data\PlaidNorm\'; end

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.72; backgroundColor = 'w'; % Panels

hFigure = figure(1);
set(hFigure,'units','normalized','outerposition',[0 0 1 1])
hFigure2 = figure(2);
set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Session(s) panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% FileNameString
[fileNameStringAll,fileNameStringListArray] = getFileNameStringList;

figure(1);
hSessionPanel= uipanel('Title','Session(s)','titleposition','centertop',...
                'fontSize',fontSizeLarge,'Unit','Normalized','Position',...
                [0.05 0.925 0.9 0.08]);
            
hSession = uicontrol('Parent',hSessionPanel,'Unit','Normalized', ...
                    'BackgroundColor', backgroundColor, 'Position', ...
                    [0.4 0.83 0.2 0.16], 'Style','popup','String',...
                    fileNameStringAll,'FontSize',fontSizeLarge);
                
hOriTunedCheckbox = uicontrol('Parent',hSessionPanel,'Unit','Normalized',...
    'Position',[0.65 0.1 0.1 0.6],'Style','checkbox','String','Ori-tuned Elecs',...
    'FontSize',fontSizeMedium);
uicontrol('Parent',hSessionPanel,'Unit','Normalized',...
    'Position',[0.75 0.1 0.1 0.8],'Style','pushbutton','String','Select Session',...
    'FontSize',fontSizeMedium,'Callback',{@selectSession_Callback});

    function selectSession_Callback(~,~)
        sessionNum = get(hSession,'val'); 
        fileNameStringTMP = fileNameStringListArray{sessionNum};
        if sessionNum <=12 || sessionNum == 23 || sessionNum == 25
            monkeyName = fileNameStringTMP{1}(1:5);
            expDate = fileNameStringTMP{1}(6:11);
            protocolName = fileNameStringTMP{1}(12:end);
        elseif sessionNum >12 && sessionNum <= 22 || sessionNum == 24
            monkeyName = fileNameStringTMP{1}(1:7);
            expDate = fileNameStringTMP{1}(8:13);
            protocolName = fileNameStringTMP{1}(14:end);

            if sessionNum>=12 && sessionNum<=22
                sessionNum = sessionNum-12; % SessionNums for monkey: kesariH
            end
            
        end
        gridType = 'Microelectrode';
        
        oriSelectiveFlag = get(hOriTunedCheckbox,'val');

        % Show electrodes on Grid
        electrodeGridPos = [0.05 panelStartHeight 0.2 panelHeight];
        hElectrodesonGrid = showElectrodeLocations(electrodeGridPos,[], ...
        [],[],1,0,gridType,monkeyName); %#ok<NASGU>
        
        if sessionNum == 25
            monkeyName = 'all';
        end
            

        %%%%%%%%%%%%%%%%%%%%%%%%%% Find Good Electrodes %%%%%%%%%%%%%%%%%%%
        [ElectrodeStringListAll,ElectrodeArrayListAll]= ...
            getElectrodesList(monkeyName,sessionNum,oriSelectiveFlag,folderSourceString);
        if sessionNum <=22 % Single Session for either monkey
            ElectrodeListSession = ElectrodeArrayListAll{1};
            ElectrodeStringSession = ElectrodeStringListAll{1};
        elseif sessionNum == 23 % all Sessions for alpaH
            ElectrodeListSession = ElectrodeArrayListAll{13};
            ElectrodeStringSession = ElectrodeStringListAll{13};
        elseif sessionNum == 24 % all Sessions for kesariH
            ElectrodeListSession = ElectrodeArrayListAll{11};
            ElectrodeStringSession = ElectrodeStringListAll{11};
        elseif sessionNum == 25 % all Sessions for both monkeys combined
            ElectrodeListSession = ElectrodeArrayListAll{25};
            ElectrodeStringSession = ElectrodeStringListAll{25};
        end


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
            'Style','popup','String',ElectrodeStringSession,...
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
        
        hRelativeMeasures = uicontrol('Parent',hParameterPanel,...
            'Unit','Normalized', ...
            'Position',[0 1-5*paramsHeight 1 paramsHeight], ...
            'Style','togglebutton',...
            'String','Show Relative Measures (Stimulus - Baseline)',...
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
        
        cValsUnique = [0 12.5 25 50 100]./2;
        cValsUnique2 = [0 12.5 25 50 100]./2;
        numRows = length(cValsUnique); numCols = length(cValsUnique2);
        gridPos=[0.02+startXPos startYPos mainTabWidth mainTabHeight];
        gap = 0.002;
        
        figure(1);
        plotHandles = getPlotHandles(numRows,numCols,gridPos,gap);
        
        % Orientation Tuning for electrodes of single session
        hOriTuning = ...
            getPlotHandles(1,1,[0.52 0.38 0.14 0.22],0.05,0.05,1);
        
        % Color Matrix of Neural measures in a 5x5 contrast conditions
        hNeuralMeasureColorMatrix = ...
            getPlotHandles(1,1,[0.52 0.05 0.14 0.25],0.001,0.001,1); 
        
        % CRF plots for preferred Orientation (Row-wise)
        plotHandles2= getPlotHandles(1,5,[0.7 0.5 0.25 0.1],0.001,0.001,1);
        
        % CRF plots for null Orientation (column-wise)
        plotHandles3= getPlotHandles(1,5,[0.7 0.35 0.25 0.1],0.001,0.001,1);
        
        % CRF plots for preferred Orientation (Overlaid)
        hRowCRF = getPlotHandles(1,1,[0.85 0.2 0.1 0.1],0.001,0.001,1);
        
        % CRF plots for null Orientation (Overlaid)
        hColumnCRF = getPlotHandles(1,1,[0.85 0.05 0.1 0.1],0.001,0.001,1);
        
        % Population histogram of Normalization Index
        hNormIndex = getPlotHandles(1,1,[0.7 0.05 0.12 0.25],0.001,0.001,1);
        
        % Defining text in new axes for Orientation Type (pref/null) [0/90; 22.5/112.5; 45/135; 67.5/157.5] 
        textH1 = getPlotHandles(1,1,[0.2 0.65 0.01 0.01]); set(textH1,'Visible','Off');
        textH2 = getPlotHandles(1,1,[0.02 0.25 0.01 0.01]); set(textH2,'Visible','Off');
        
        % plots for 2nd annual Work Presentation (3rd Year)  
        figure(2);
        hPlotPreferred = getPlotHandles(1,5,[0.1 0.7 0.5 0.15],0.01,0.01,1); linkaxes(hPlotPreferred);
        hOtherMeaures = getPlotHandles(1,3,[0.1 0.3 0.5 0.25],0.05,0.05,0);
        
        % Defining Text  in new axes for Figure 2
        textH3 = getPlotHandles(1,1,[0.2 0.9 0.01 0.01]); set(textH3,'Visible','Off');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        freqRanges{1} = [8 12]; % alpha
        freqRanges{2} = [30 80]; % gamma
        freqRanges{3} = [16 16];  % SSVEP
        
%         freqRangeStr = {'alpha','gamma','SSVEP'};
%         numFreqRanges = length(freqRanges);        
        
        % Plotting Functions
        function plotData_Callback(~,~)
            
            %%%%%%%%%%%%%%%%%%%%%%%% Read values %%%%%%%%%%%%%%%%%%%%%%%%%%
            electrodeString = get(hElectrode,'val'); 
            if length(fileNameStringTMP) ==1
                ElectrodeListTMP = ElectrodeListSession(electrodeString);
                if isempty(ElectrodeListTMP{1})
                    error(['No electrode found for analysis!'...
                            'Please try another session!'])
                end
                folderName = fullfile(folderSourceString,'data',monkeyName,gridType,expDate,protocolName);
                folderExtract = fullfile(folderName,'extractedData');
                [~,~,~,~,~,~,oValsUnique,~,~,~,~,~,~,oValsUnique2,~,~] = loadParameterCombinations(folderExtract);

            elseif length(fileNameStringTMP)>1
                ElectrodeListTMP = ElectrodeListSession;
            end
  
            analysisMethod = get(hAnalysisMethod,'val');
            analysisMeasure = get(hAnalysisType,'val');
            NormalizeDataFlag = get(hNormalizeData,'val');
            relativeMeasuresFlag = get(hRelativeMeasures,'val');
            
            erpRange = [str2double(get(hERPMin,'String'))...
                        str2double(get(hERPMax,'String'))];
            blRange = [str2double(get(hBaselineMin,'String'))...
                       str2double(get(hBaselineMax,'String'))];
            stRange = [str2double(get(hStimPeriodMin,'String'))...
                       str2double(get(hStimPeriodMax,'String'))];

            plotColor = colorNames(get(hChooseColor,'val'));
            holdOnState = get(hHoldOn,'val'); %#ok<NASGU>
            
            ColorNeuralMeasures = jet(5);

            % get Data for Selected Session & Parameters
            [erpData,firingRateData,fftData,energyData,oriTuningData,~] = getData(folderSourceString,...
             fileNameStringTMP,ElectrodeListTMP,erpRange,blRange,...
             stRange,freqRanges); 
             
            if analysisMeasure == 1 % computing ERP
                plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,hNormIndex,hOriTuning,sessionNum,erpData.timeVals,erpData,oriTuningData,plotColor,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
                plotDataFig2(hPlotPreferred,hOtherMeaures,sessionNum,erpData.timeVals,erpData,ColorNeuralMeasures,analysisMethod,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
            elseif analysisMeasure == 2 % computing Firing rate
                plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,hNormIndex,hOriTuning,sessionNum,firingRateData.timeVals,firingRateData,oriTuningData,plotColor,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
                plotDataFig2(hPlotPreferred,hOtherMeaures,sessionNum,firingRateData.timeVals,firingRateData,ColorNeuralMeasures,analysisMethod,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
            elseif analysisMeasure == 3 % computing Raster Plot from spike data
                error('Still working on raster data!')
            elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6 % computing alpha
                if analysisMethod == 1
                    plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,hNormIndex,hOriTuning,sessionNum,fftData.freqVals,fftData,oriTuningData,plotColor,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
                    plotDataFig2(hPlotPreferred,hOtherMeaures,sessionNum,fftData.freqVals,fftData,ColorNeuralMeasures,analysisMethod,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
                elseif analysisMethod ==2 
                    plotData(plotHandles,hNeuralMeasureColorMatrix,plotHandles2,plotHandles3,hRowCRF,hColumnCRF,hNormIndex,hOriTuning,sessionNum,energyData.freqVals,energyData,oriTuningData,plotColor,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
                    plotDataFig2(hPlotPreferred,hOtherMeaures,sessionNum,energyData.freqVals,energyData,ColorNeuralMeasures,analysisMethod,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
                end
            elseif analysisType == 7 % need to work on STA!
                error('STA computation method not found') 
            end
            
            figure(1);
            % Text for Orientation for main 5x5 Neural measure matrix
            textH1 = getPlotHandles(1,1,[0.2 0.65 0.01 0.01]);
            set(textH1,'Visible','Off');
            textH2 = getPlotHandles(1,1,[0.02 0.25 0.01 0.01]);
            set(textH2,'Visible','Off');
                
            if length(fileNameStringTMP) ==1 && length(ElectrodeListTMP{1})==1    
                text(0.35,1.15,['Null Orientation: ' num2str(oValsUnique)],'unit','normalized','fontsize',20,'fontweight','bold','rotation',90,'parent',textH2);
                text(0.35,1.15,['Preferred Orientation: ' num2str(oValsUnique2)],'unit','normalized','fontsize',20,'fontweight','bold','parent',textH1);
            else
                text(0.35,1.15,'Null Orientation','unit','normalized','fontsize',20,'fontweight','bold','rotation',90,'parent',textH2);
                text(0.35,1.15,'Preferred Orientation' ,'unit','normalized','fontsize',20,'fontweight','bold','parent',textH1);
            end
            
            flippedcValsUnique2 = flip(cValsUnique2);
            for c = 1:5
            title(plotHandles(1,c),[num2str(cValsUnique(c)) ' %']);
            ylabel(plotHandles(c,1),[num2str(flippedcValsUnique2(c)) ' %'],'fontWeight','bold');
            end
            
            % Setting xlabels and ylabels for Color matrix of neural
            % measures
            set(hNeuralMeasureColorMatrix,'XTickLabel',cValsUnique);
            set(hNeuralMeasureColorMatrix,'YTickLabel',flip(cValsUnique2));
            
            % Setting title,xlabel and ylabel for OriTuning Plot
            title(hOriTuning,'Ori Tuning for Single Session (Spike Data)');
%             set(hOriTuning,'XTickLabel',oValsUnique_Tuning);


            % Setting Plot Ranges
            if analysisMeasure<=3 %ERP or spikes
                xMin = str2double(get(hStimMin,'String'));
                xMax = str2double(get(hStimMax,'String'));
            elseif analysisMeasure <=6  % LFP fft analysis
                xMin = str2double(get(hFFTMin,'String'));
                xMax = str2double(get(hFFTMax,'String'));                
            else
                xMin = str2double(get(hSTAMin,'String'));
                xMax = str2double(get(hSTAMax,'String'));
            end
            
            figure(2)
            % Text for Orientation for main 5x5 Neural measure matrix
%             textH3 = getPlotHandles(1,1,[0.2 0.65 0.01 0.01]);
%             set(textH3,'Visible','Off');
            text(0.35,1.15,'Contrast along Preferred Orientation','unit','normalized','fontsize',20,'fontweight','bold','parent',textH3);

            rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
            rescaleData(plotHandles2,0,50,getYLims(plotHandles2));
            rescaleData(plotHandles3,0,50,getYLims(plotHandles3));
            rescaleData(hRowCRF,0,50,getYLims(hRowCRF));
            rescaleData(hColumnCRF,0,50,getYLims(hColumnCRF));
            rescaleData(hOriTuning,0,160,getYLims(hOriTuning));
            rescaleData(hPlotPreferred,xMin,xMax,getYLims(hPlotPreferred));
            rescaleData(hOtherMeaures(1),0,50,getYLims(hOtherMeaures(1)));
            

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
        rescaleData(plotHandles2,0,50,getYLims(plotHandles2));
        rescaleData(plotHandles3,0,50,getYLims(plotHandles3));

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
        claGivenPlotHandle(hOriTuning);
        claGivenPlotHandle(hNormIndex);
        delete(findobj(textH1,'type','text'));
        delete(findobj(textH2,'type','text'));
        delete(findobj(textH3,'type','text'));
        
        claGivenPlotHandle(hPlotPreferred);
        claGivenPlotHandle(hOtherMeaures);
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
function[erpData,firingRateData,fftData,energyData,oriTuningData,electrodeArray] = ...
    getData(folderSourceString,fileNameStringTMP,ElectrodeListTMP,...
    erpRange,blRange,stRange,freqRanges)

numDatasets = length(fileNameStringTMP);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[erpData,firingRateData,fftData,energyData,oriTuningData,electrodeArray]...
= getDataSingleSession(folderSourceString,fileNameStringTMP{1},...
ElectrodeListTMP{1},erpRange,blRange,stRange,freqRanges); 

if length(fileNameStringTMP)>1
    for i=2:numDatasets
        if isempty(ElectrodeListTMP{i})
           continue
        end
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [erpDataTMP,firingRateDataTMP,fftDataTMP,energyDataTMP,~,~] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},...
            ElectrodeListTMP{i},erpRange,blRange,stRange,freqRanges);
        
        erpData.data = cat(1,erpData.data,erpDataTMP.data);
        erpData.analysisDataBL = cat(1,erpData.analysisDataBL,erpDataTMP.analysisDataBL);
        erpData.analysisDataST = cat(1,erpData.analysisDataST,erpDataTMP.analysisDataST);
        
        firingRateData.data = cat(1,firingRateData.data,firingRateDataTMP.data);
        firingRateData.analysisDataBL = cat(1,firingRateData.analysisDataBL,firingRateDataTMP.analysisDataBL);
        firingRateData.analysisDataST = cat(1,firingRateData.analysisDataST,firingRateDataTMP.analysisDataST);
        
        fftData.dataBL = cat(1,fftData.dataBL,fftDataTMP.dataBL);
        fftData.dataST = cat(1,fftData.dataST,fftDataTMP.dataST);
        for j = 1:3
            fftData.analysisDataBL{j} = cat(1,fftData.analysisDataBL{j},fftDataTMP.analysisDataBL{j});
            fftData.analysisDataST{j} = cat(1,fftData.analysisDataST{j},fftDataTMP.analysisDataST{j});
        end
        
        energyData.dataBL = cat(1,energyData.dataBL,energyDataTMP.dataBL);
        energyData.dataST = cat(1,energyData.dataST,energyDataTMP.dataST);
        for j =1:3
            energyData.analysisDataBL{j} = cat(1,energyData.analysisDataBL{j},energyDataTMP.analysisDataBL{j});
            energyData.analysisDataST{j} = cat(1,energyData.analysisDataST{j},energyDataTMP.analysisDataST{j});
        end

% Combining OriData across sessions need to be done!        
%         oriTuningDataTMP.PO = 
%         oriTuningDataTMP.OS
%         oriTuningDataTMP.FR
        
        electrodeArray = [];
    end
end
end



function [erpData,firingRateData,fftData,energyData,oriTuningData,electrodeArray] = ...
    getDataSingleSession(folderSourceString,fileNameStringTMP,...
    ElectrodeListTMP,erpRange,blRange,stRange,freqRanges)

if strcmp(fileNameStringTMP(1:5),'alpaH')       
    monkeyName = 'alpaH'; 
    expDate = fileNameStringTMP(6:11); 
    protocolName = fileNameStringTMP(12:end);
    oriTuning_protocolName = ['GRF_00' num2str(str2double(protocolName(5:end))-1)]; % The protocol Number is just the immediate precedent of the main protocol 
elseif strcmp(fileNameStringTMP(1:7),'kesariH')
    monkeyName = 'kesariH';
    expDate = fileNameStringTMP(8:13); 
    protocolName = fileNameStringTMP(14:end);
    oriTuning_protocolName = ['GRF_00' num2str(str2double(protocolName(5:end))-1)]; % The protocol Number is just the immediate precedent of the main protocol 
end
gridType = 'microelectrode';


folderName = fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,protocolName);
tuningProtocol_folderName =  fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,oriTuning_protocolName);                           

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');

folderSave = fullfile(tuningProtocol_folderName,'savedData');
if ~exist(folderSave,'dir')
    mkdir(folderSave);
end
fileToSave = fullfile(folderSave,['oriTuningData_' num2str(1000*stRange(1)) 'ms_' num2str(1000*stRange(2)) 'ms.mat']);

if exist(fileToSave,'file')
    disp(['Loading file ' fileToSave]);
    load(fileToSave);
else
    % Get OrientationTuning Data
    [computationVals,PO,OS] = savePrefOriAndOriSelectivitySpikes(monkeyName,expDate,oriTuning_protocolName,folderSourceString,gridType);
end

    oriTuningData.PO = PO(ElectrodeListTMP);
    oriTuningData.OS = OS(ElectrodeListTMP);
    oriTuningData.FR = computationVals(ElectrodeListTMP,:);

    % Get Combinations
    [parameterCombinations,parameterCombinations2,...
        aValsUnique,eValsUnique,~,~,oValsUnique,cValsUnique,tValsUnique, ...
        aValsUnique2,eValsUnique2,~,~,oValsUnique2,cValsUnique2,tValsUnique2] = ...
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

    c1List = length(cValsUnique):-1:1; %  Flipping data Row-wise so that positive x-axis and positive y-axis denotes increase in Contrast
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
                        disp(['pos=(' num2str(c1List(c1)) ',' num2str(c2List(c2)) ') ,n=' num2str(length(goodPos))]);
                        N(iElec,t,c1List(c1),c2List(c2)) = length(goodPos);

                        if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
                            disp('baseline and stimulus ranges are not the same');
                        else
                           
                           erp = mean(analogData(goodPos,:),1); %#ok<NODEF>
                           erpDataTMP(iElec,t,c1,c2,:) = erp;
                           RMSvalsBL(iElec,t,c1,c2) = rms(erp(blPos));
                           RMSvalsERP(iElec,t,c1,c2) = rms(erp(erpPos));
                            
                           [psthData(iElec,t,c1,c2,:),xsFR] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
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

    fftData.dataBL = fftDataBL;
    fftData.dataST = fftDataST;
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
    
    elecs_neededtoFlipped = find(abs(oriTuningData.PO-oValsUnique)<abs(oriTuningData.PO-oValsUnique2));
    erpData = segregate_Pref_Null_data(erpData,elecs_neededtoFlipped);
    firingRateData = segregate_Pref_Null_data(firingRateData,elecs_neededtoFlipped);
    fftData = segregate_Pref_Null_data(fftData,elecs_neededtoFlipped);
    energyData = segregate_Pref_Null_data(energyData,elecs_neededtoFlipped);
    

%     % Save Data for particular session
%     save(fileToSave,'oriTuningData');
% end
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

function [ElectrodeStringListAll,ElectrodeArrayListAll] = getElectrodesList(monkeyName,sessionNum,oriSelectiveFlag,folderSourceString)

[tmpElectrodeStringList,tmpElectrodeArrayList,allElecs,monkeyNameList] = getGoodElectrodesDetails(monkeyName,sessionNum,oriSelectiveFlag,folderSourceString);

    ElectrodeStringListAll = ''; 
    ElectrodeArrayListAll = [];
    
    for i=1:length(monkeyNameList)
            ElectrodeStringListAll = cat(2,ElectrodeStringListAll,tmpElectrodeStringList{i});
            ElectrodeArrayListAll = cat(2,ElectrodeArrayListAll,tmpElectrodeArrayList{i});

    end
    if sessionNum >=23
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
    end
    AllSessions = length(ElectrodeArrayListAll)+1;
    
    if sessionNum == 25
        for k=1:length(tmpElectrodeArrayList{1})+length(tmpElectrodeArrayList{2})
            ElectrodeArrayListAll{1,AllSessions}{1,k} = ElectrodeArrayListAll{1,k}{1,end};
            ElectrodeStringListAll{1,AllSessions} = ['all (N=' num2str(sum(allElecs)) ')'];
        end
    end
    
end


function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange)

posToAverage = intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2)));
eValue   = mean(mEnergy(posToAverage));
end
function normData = normalizeData(x)
for iElec = 1:size(x.data,1)
    for t = 1:size(x.data,2)
        normData.data(iElec,t,:,:,:) = x.data(iElec,t,:,:,:)./max(max(max(abs(x.data(iElec,t,:,:,:)))));
        normData.analysisDataBL(iElec,t,:,:) = x.analysisDataBL(iElec,t,:,:)./max(max(abs(x.analysisDataBL(iElec,t,:,:))));
        normData.analysisDataST(iElec,t,:,:) = x.analysisDataST(iElec,t,:,:)./max(max(abs(x.analysisDataST(iElec,t,:,:))));
        normData.timeVals = x.timeVals;
        normData.N = x.N;
    end
end
end
function plotData(hPlot1,hPlot2,hPlot3,hPlot4,hPlot5,hPlot6,hPlot7,hPlot8,sessionNum,xs,data,oriTuningData,colorName,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)

% Main 5x5 plot for Neural Measure
if analysisMeasure == 1 || analysisMeasure == 2
    if NormalizeDataFlag
    % Normalize ERP and spike Data (fft or energy data need not be normalized
    % as they are expressed in log units)
    data = normalizeData(data);
    end
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
            dataPlotBL = squeeze(mean(squeeze(data.dataBL(:,2,:,:,:)),1));
            dataPlotST = squeeze(mean(squeeze(data.dataST(:,2,:,:,:)),1));
        end
    end
    % When Change in neural measures are to be plotted
    if relativeMeasuresFlag
        if analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6
            dataPlotdiffSTvsBL = 10*(dataPlotST-dataPlotBL); % Change in Power expressed in deciBel
        else
            dataPlotdiffSTvsBL = dataPlotST-dataPlotBL; 
        end
    end
end
    
% Plotting 5x5 plots for raw Neural Measures
for c1 = 1:5
    for c2 = 1:5
        if analysisMeasure == 1 || analysisMeasure == 2
            plot(hPlot1(c1,c2),xs,squeeze(dataPlot(c1,c2,:)),'color',colorName);
        elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6
            if ~relativeMeasuresFlag
            plot(hPlot1(c1,c2),xs,squeeze(dataPlotBL(c1,c2,:)),'g');
            hold(hPlot1(c1,c2),'on')
            plot(hPlot1(c1,c2),xs,squeeze(dataPlotST(c1,c2,:)),'k');
            hold(hPlot1(c1,c2),'off')
            else
                plot(hPlot1(c1,c2),xs,squeeze(dataPlotdiffSTvsBL(c1,c2,:)),'b');
            end
    
        end
    end
end

% Color coded 5x5 Matrix of raw Neural Measures

if dataSize(1)==1
    if analysisMeasure == 1 || analysisMeasure == 2
        analysisData = squeeze(data.analysisDataST(:,1,:,:));
        analysisDataBL = squeeze(data.analysisDataBL(:,1,:,:));
    elseif analysisMeasure == 4 
        analysisData = squeeze(data.analysisDataST{1}(:,1,:,:));
        analysisDataBL = squeeze(data.analysisDataBL{1}(:,1,:,:));
    elseif analysisMeasure == 5
        analysisData = squeeze(data.analysisDataST{2}(:,1,:,:));
        analysisDataBL = squeeze(data.analysisDataBL{2}(:,1,:,:));
    elseif analysisMeasure == 6
        analysisData = squeeze(data.analysisDataST{3}(:,2,:,:));
        analysisDataBL = squeeze(data.analysisDataBL{3}(:,2,:,:));
    end
elseif dataSize(1)>1
    if analysisMeasure == 1 || analysisMeasure == 2
        analysisData = squeeze(mean(squeeze(data.analysisDataST(:,1,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL(:,1,:,:)),1));
        for iElec= 1:size(data.analysisDataST,1)
            clear electrodeVals
            electrodeVals =  squeeze(data.analysisDataST(iElec,1,:,:));
            NI_population(iElec) = abs((electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5));
        end
    elseif analysisMeasure == 4
        analysisData = squeeze(mean(squeeze(data.analysisDataST{1}(:,1,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL{1}(:,1,:,:)),1));
        for iElec= 1:size(data.analysisDataST{1},1)
            clear electrodeVals
            electrodeVals =  squeeze(data.analysisDataST{1}(iElec,1,:,:));
            NI_population(iElec) = abs((electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5));
        end        
    elseif analysisMeasure == 5
        analysisData = squeeze(mean(squeeze(data.analysisDataST{2}(:,1,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL{2}(:,1,:,:)),1));
        for iElec= 1:size(data.analysisDataST{2},1)
            clear electrodeVals
            electrodeVals =  squeeze(data.analysisDataST{2}(iElec,1,:,:));
            NI_population(iElec) = abs((electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5));
        end        
    elseif analysisMeasure == 6
        analysisData = squeeze(mean(squeeze(data.analysisDataST{3}(:,2,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL{3}(:,2,:,:)),1));
        for iElec= 1:size(data.analysisDataST{3},1)
            clear electrodeVals
            electrodeVals =  squeeze(data.analysisDataST{3}(iElec,2,:,:));
            NI_population(iElec) = abs((electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5));
        end        
    end
end
if relativeMeasuresFlag
    if analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6
        analysisData = 10*(analysisData-analysisDataBL); % Change in power expressed in deciBel
    else
        analysisData = analysisData-analysisDataBL; % No need to do this for ERP & firing Rate
    end
end

imagesc(analysisData,'parent',hPlot2);colorbar(hPlot2);set(hPlot2,'Position',[0.52 0.05 0.12 0.25]);
NIAnalysisData = analysisData;
NormIndex = (NIAnalysisData(1,1)+ NIAnalysisData(5,5))/NIAnalysisData(1,5);
title(hPlot2,['NI: ',num2str(NormIndex)],'fontWeight','bold');

if size(data.analysisDataST,1) == 1
    if sessionNum <=12
        colorNamesOriTuning = hsv(30);
        oValsUnique_Tuning = [0 22.5 45 67.5 90 112.5 135 157.5];
        for index = 1:size(data.analysisDataST,1)
        plot(hPlot8,oValsUnique_Tuning,oriTuningData.FR(index,:),'Marker','o','color',colorNamesOriTuning(index,:,:));
        text(0.05,index*0.1+0.5,['PO: ' num2str(oriTuningData.PO(index)) ',OS: ' num2str(oriTuningData.OS(index))],'color',colorNamesOriTuning(index,:,:),'unit','normalized','parent',hPlot8);
        hold(hPlot8,'on');
        end
        hold(hPlot8,'off');
    else
    end
end

if sessionNum>12
% NI population histogram
histogram(hPlot7,NI_population);
end

% Contrast Response curves Row-Wise & Column-wise
cValsUnique = [0 12.5 25 50 100]./2;
cValsUnique2 = [0 12.5 25 50 100]./2;
conFlipped = 5:-1:1;
CRFColors = jet(length(cValsUnique));
for iCon = 1:5
    plot(hPlot3(1,iCon),cValsUnique,analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    plot(hPlot4(1,iCon),cValsUnique2,analysisData(:,conFlipped(iCon)),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    plot(hPlot5(1,1),cValsUnique,analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlot5(1,1),'on')
    plot(hPlot6(1,1),cValsUnique2,analysisData(:,conFlipped(iCon)),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlot6(1,1),'on')
end
hold(hPlot5(1,1),'off'); hold(hPlot6(1,1),'off')
end


function plotDataFig2(hPlot1,hPlot2,sessionNum,xs,data,colors,analysisMethod,analysisMeasure,relativeMeasuresFlag,NormalizeDataFlag)
cValsUnique = [0 12.5 25 50 100]./2;
% Main 5x5 plot for Neural Measure
if analysisMeasure == 1 || analysisMeasure == 2
    if NormalizeDataFlag
    % Normalize ERP and spike Data (fft or energy data need not be normalized
    % as they are expressed in log units)
    data = normalizeData(data);
    end
    dataSize = size(data.data);
    if dataSize(1) == 1
        dataPlot = squeeze(squeeze(data.data(:,1,:,:,:)));
    elseif dataSize(1) >1
        dataPlot = squeeze(mean(squeeze(data.data(:,1,:,:,:)),1));
    end
elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6
    if size(data.dataBL) == size(data.dataST) 
        dataSize = size(data.dataST);
%         data.dataBL = 
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
        dataPlotBL = squeeze(mean(mean(mean(squeeze(data.dataBL(:,1,:,:,:)),3),2),1)); %CommonBaseline
        dataPlotBL = reshape(repmat(dataPlotBL',[25 1]),[5 5 size(dataPlotBL)]);
        dataPlotST = squeeze(mean(squeeze(data.dataST(:,1,:,:,:)),1));
        if analysisMeasure == 6
            dataPlotBL = squeeze(mean(mean(mean(squeeze(data.dataBL(:,2,:,:,:)),3),2),1)); %CommonBaseline
            dataPlotBL = reshape(repmat(dataPlotBL',[25 1]),[5 5 size(dataPlotBL)]);
            dataPlotST = squeeze(mean(squeeze(data.dataST(:,2,:,:,:)),1));
        end
    end
    % When Change in neural measures are to be plotted
    if relativeMeasuresFlag
        if analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6
            dataPlotdiffSTvsBL = dataPlotST-dataPlotBL; % Change in Power expressed in deciBel
        else
            dataPlotdiffSTvsBL = dataPlotST-dataPlotBL; 
        end
    end
end

% Plotting 5x5 plots for raw Neural Measures
cFlipped = 5:-1:1;
for c1 = 1:5
    for c2 = 1:5
        if analysisMeasure == 1 || analysisMeasure == 2
            plot(hPlot1(1,c2),xs,squeeze(dataPlot(cFlipped(c1),c2,:)),'color',colors(c1,:,:),'LineWidth',2);
            
            hold(hPlot1(1,c2),'on');
            
        elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6
            if ~relativeMeasuresFlag
            plot(hPlot1(1,c2),xs,squeeze(dataPlotBL(cFlipped(c1),c2,:)),'k');
            hold(hPlot1(1,c2),'on')
            plot(hPlot1(1,c2),xs,squeeze(dataPlotST(cFlipped(c1),c2,:)),'color',colors(c1,:,:),'LineWidth',2);
            
            else
                plot(hPlot1(1,c2),xs,squeeze(dataPlotBL(cFlipped(c1),c2,:))-squeeze(dataPlotBL(cFlipped(c1),c2,:)),'k');
                plot(hPlot1(1,c2),xs,squeeze(dataPlotdiffSTvsBL(cFlipped(c1),c2,:)),'color',colors(c1,:,:),'LineWidth',2);
            end
            
        end
        
%         text(0.05,c1*0.12+0.3,['Null: ',num2str(cValsUnique(c1)) ' %'],'color',colors(c1,:,:),'unit','normalized','parent',hPlot1(1))
        title(hPlot1(1,c2),[num2str(cValsUnique(c2)) ' %'])
    end
end
hold(hPlot1(1,c2),'off');
% % Color coded 5x5 Matrix of raw Neural Measures

if dataSize(1)==1
    if analysisMeasure == 1 || analysisMeasure == 2
        analysisData = squeeze(data.analysisDataST(:,1,:,:));
        analysisDataBL = squeeze(data.analysisDataBL(:,1,:,:));
    elseif analysisMeasure == 4 
        analysisData = squeeze(data.analysisDataST{1}(:,1,:,:));
        analysisDataBL = squeeze(data.analysisDataBL{1}(:,1,:,:));
    elseif analysisMeasure == 5
        analysisData = squeeze(data.analysisDataST{2}(:,1,:,:));
        analysisDataBL = squeeze(data.analysisDataBL{2}(:,1,:,:));
    elseif analysisMeasure == 6
        analysisData = squeeze(data.analysisDataST{3}(:,2,:,:));
        analysisDataBL = squeeze(data.analysisDataBL{3}(:,2,:,:));
    end
elseif dataSize(1)>1
    if analysisMeasure == 1 || analysisMeasure == 2
        analysisData = squeeze(mean(squeeze(data.analysisDataST(:,1,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL(:,1,:,:)),1));
        for iElec= 1:size(data.analysisDataST,1)
            clear electrodeVals
            if ~relativeMeasuresFlag
                electrodeVals =  squeeze(data.analysisDataST(iElec,1,:,:));
            else
                electrodeVals =  squeeze(data.analysisDataST(iElec,1,:,:))-squeeze(data.analysisDataBL(iElec,1,:,:));
            end
            NI_population(iElec) = (electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5);
        end
    elseif analysisMeasure == 4
        analysisData = squeeze(mean(squeeze(data.analysisDataST{1}(:,1,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL{1}(:,1,:,:)),1));
        for iElec= 1:size(data.analysisDataST{1},1)
            clear electrodeVals
            if ~relativeMeasuresFlag
                electrodeVals =  squeeze(data.analysisDataST{1}(iElec,1,:,:));
            else
                electrodeVals =  squeeze(data.analysisDataST{1}(iElec,1,:,:))-squeeze(data.analysisDataBL{1}(iElec,1,:,:));
            end
            NI_population(iElec) = (electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5);
        end        
    elseif analysisMeasure == 5
        analysisData = squeeze(mean(squeeze(data.analysisDataST{2}(:,1,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL{2}(:,1,:,:)),1));
        for iElec= 1:size(data.analysisDataST{2},1)
            clear electrodeVals
            if ~relativeMeasuresFlag
                electrodeVals =  squeeze(data.analysisDataST{2}(iElec,1,:,:));
            else
                electrodeVals =  squeeze(data.analysisDataST{2}(iElec,1,:,:))-squeeze(data.analysisDataBL{2}(iElec,1,:,:));

            end
            NI_population(iElec) = (electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5);
        end        
    elseif analysisMeasure == 6
        analysisData = squeeze(mean(squeeze(data.analysisDataST{3}(:,2,:,:)),1));
        analysisDataBL = squeeze(mean(squeeze(data.analysisDataBL{3}(:,2,:,:)),1));
        for iElec= 1:size(data.analysisDataST{3},1)
            clear electrodeVals
            if ~relativeMeasuresFlag
                electrodeVals =  squeeze(data.analysisDataST{3}(iElec,2,:,:));
            else
                electrodeVals =  squeeze(data.analysisDataST{3}(iElec,2,:,:))-squeeze(data.analysisDataBL{3}(iElec,2,:,:));
            end
            NI_population(iElec) = (electrodeVals(1,1)+electrodeVals(5,5))/electrodeVals(1,5);
        end        
    end
end
if relativeMeasuresFlag
    if analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6
        analysisData = analysisData-analysisDataBL; % Change in power expressed in deciBel
    else
        analysisData = analysisData-analysisDataBL; % No need to do this for ERP & firing Rate
    end
end

imagesc(analysisData,'parent',hPlot2(3));
% grid on;
color_Bar = colorbar(hPlot2(3));% 
colorYlabelHandle = get(color_Bar,'Ylabel');
if analysisMeasure==2
    YlabelString = 'Spikes/s';
elseif analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6
    if ~relativeMeasuresFlag
        YlabelString = 'log_1_0(FFT Amplitude)';
    else
        YlabelString = 'log_1_0(\Delta FFT Amplitude)';
    end
    
end

set(colorYlabelHandle,'String',YlabelString,'fontSize',14);

plotPos = get(hPlot2(3),'Position');
set(hPlot2(3),'Position',[plotPos(1) plotPos(2) 0.12 plotPos(4)]);
NIAnalysisData = analysisData;
NormIndex = (NIAnalysisData(1,1)+ NIAnalysisData(5,5))/NIAnalysisData(1,5);
title(hPlot2(3),['NI: ',num2str(NormIndex)],'fontWeight','bold');
% 
% if size(data.analysisDataST,1) == 1
%     if sessionNum <=12
%         colorNamesOriTuning = hsv(30);
%         oValsUnique_Tuning = [0 22.5 45 67.5 90 112.5 135 157.5];
%         for index = 1:size(data.analysisDataST,1)
%         plot(hPlot8,oValsUnique_Tuning,oriTuningData.FR(index,:),'Marker','o','color',colorNamesOriTuning(index,:,:));
%         text(0.4,index*0.07+0.5,['PO: ' num2str(oriTuningData.PO(index)) ',OS: ' num2str(oriTuningData.OS(index))],'color',colorNamesOriTuning(index,:,:),'unit','normalized','parent',hPlot8);
%         hold(hPlot8,'on');
%         end
%         hold(hPlot8,'off');
%     else
%     end
% end
% 
if sessionNum>12
% NI population histogram
histogram(hPlot2(2),NI_population);
end

% Contrast Response curves Row-Wise & Column-wise

% flippedcVals = flip(cValsUnique);
% cValsUnique2 = [0 12.5 25 50 100]./2;
conFlipped = 5:-1:1;
CRFColors = jet(length(cValsUnique));
for iCon = 1:5
    plot(hPlot2(1),cValsUnique,analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlot2(1),'on')
end
hold(hPlot2(1),'off');

if analysisMeasure == 1||analysisMeasure == 2
    displayRange(hPlot1,[0.2 0.4],getYLims(hPlot1),'k');
elseif analysisMeasure ==4 
    displayRange(hPlot1,[8 12],getYLims(hPlot1),'k');
elseif analysisMeasure == 5  
    displayRange(hPlot1,[30 80],getYLims(hPlot1),'k');
elseif analysisMeasure == 6     
    displayRange(hPlot1,[16 16],getYLims(hPlot1),'k');
end

tickLengthPlot = 2*get(hPlot2(1),'TickLength');

for idx =1:5
set(hPlot1(idx),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
if analysisMeasure == 1||analysisMeasure == 2
    xlabel(hPlot1(1),'Time(ms)');
elseif analysisMeasure == 4||analysisMeasure == 5|| analysisMeasure == 6
    xlabel(hPlot1(1),'Frequency(Hz)');
end

if analysisMeasure == 1
    ylabel(hPlot1(1),'ERP (\mu V)');ylabel(hPlot2(1),'RMS value of ERP');
elseif analysisMeasure == 2 
    ylabel(hPlot1(1),'Spikes/s');ylabel(hPlot2(1),'spikes/s');
elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6
    if analysisMethod ==1
        ylabel(hPlot1(1),'log_1_0 (FFT Amplitude)'); ylabel(hPlot2(1),'log_1_0 (FFT Amplitude)');
        if relativeMeasuresFlag
            ylabel(hPlot1(1),'log_1_0 (\Delta FFT Amplitude)'); ylabel(hPlot2(1),'log_1_0 (\Delta FFT Amplitude)')
        end

    elseif analysisMethod ==2
        ylabel(hPlot1(1),'log_1_0 (Power)');ylabel(hPlot2(1),'log_1_0 (Power)')
        if relativeMeasuresFlag
            ylabel(hPlot1(1),'Change in Power (dB)');  ylabel(hPlot2(1),'Change in Power (dB)');
        end
    end

end



%             set(hOtherMeaures(1),'XScale','linear')
text(0.5,0.3,['N = ' num2str(dataSize(1))],'color','k','unit','normalized','fontSize',14,'fontWeight','bold','parent',hPlot2(1))
set(hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot2(1),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlot2(1),'Contrast (%)');
title(hPlot2(1),'CRF along Pref Ori')


set(hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
xlabel(hPlot2(2),'Normalization Index'); ylabel(hPlot2(2),'No. of Electrodes');
title(hPlot2(2),'Population NI Histogram')

set(hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlot2(3),'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlot2(3),'Contrast (%)');ylabel(hPlot2(3),'Contrast (%)');

end

function data = segregate_Pref_Null_data(data,elecs_neededtoFlipped)
    for iElec = 1:length(elecs_neededtoFlipped)
        for iTF = 1:2
            disp ([num2str(iElec), ' ' num2str(iTF)])
            if numel(fieldnames(data)) == 5
                data.data(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.data(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
            elseif numel(fieldnames(data)) == 6
                data.dataBL(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.dataBL(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
                data.dataST(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.dataST(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
            end
            
            if ~iscell(data.analysisDataBL) && ~iscell(data.analysisDataST)
                data.analysisDataBL(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataBL(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                data.analysisDataST(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataST(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
            elseif iscell(data.analysisDataBL) && iscell(data.analysisDataST)
                for iCell = 1:length(data.analysisDataBL)
                    data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                    data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                end
            end
        end
    end
  
end
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored); %#ok<NODEF>
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

function displayRange(plotHandles,range,yLims,colorName)
[nX,nY] = size(plotHandles);
%yLims = getYLims(plotHandles);

yVals = yLims(1):(yLims(2)-yLims(1))/100:yLims(2);
xVals1 = range(1) + zeros(1,length(yVals));
xVals2 = range(2) + zeros(1,length(yVals));

for i=1:nX
    for j=1:nY
        hold(plotHandles(i,j),'on');
        plot(plotHandles(i,j),xVals1,yVals,'color',colorName);
        plot(plotHandles(i,j),xVals2,yVals,'color',colorName);
    end
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
labelSize=14;
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




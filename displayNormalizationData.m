function displayNormalizationData(folderSourceString)

if ~exist('folderSourceString','var') 
   folderSourceString = 'M:\data\PlaidNorm\';
end

close all; % closes any open figure to avoid any overlaying issues

% Display Options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16; % Fonts
panelHeight = 0.2; panelStartHeight = 0.68; backgroundColor = 'w'; % Panels

hFigure = figure(1);
set(hFigure,'units','normalized','outerposition',[0 0 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Session(s) panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% FileNameString
[fileNameStringAll,~,fileNameStringListArray] = getFileNameStringList;

figure(1);
hSessionPanel= uipanel('Title','Session Panel','titleposition','centertop',...
                'fontSize',fontSizeLarge,'Unit','Normalized','Position',...
                [0.05 0.89 0.225 0.11]);
            
hSession = uicontrol('Parent',hSessionPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [0 0.9 1 0.12], 'HorizontalAlignment','Center','Style',...
        'popup','String',...
        fileNameStringAll,'FontSize',fontSizeLarge);
    
% To align the strings inside the pop-up menu, a separate function jobj 
% is required which uses JAVA graphics properties   
% jobj = findjobj(hSession);
% renderer = jobj.getRenderer();
% renderer.setHorizontalAlignment(renderer.CENTER);
    
hOriTunedCheckbox = uicontrol('Parent',hSessionPanel,'Unit','Normalized',...
'Position',[0.18 0.1 0.8 0.3],'Style','checkbox',...
'String','Orientation tuned Electrodes','FontSize',fontSizeMedium);
                           
hTimingPanel = uipanel('Title','Timing Panel','titleposition','centertop',...
                'fontSize',fontSizeLarge,'Unit','Normalized','Position',...
                [0.275 0.89 0.225 0.11]);
            
timingTextWidth = 0.5; timingBoxWidth = 0.2;
timingHeight = 1/3; 

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight 0.5 timingHeight],...
    'Style','text','String','Parameter','FontSize',fontSizeSmall);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth...
    timingHeight], 'Style','text','String','Min',...
    'FontSize',fontSizeSmall);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight ...
    timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeSmall);
                
% Baseline Range
baseline = [-0.25 0];
uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
    'Position',[0 1-2*timingHeight 0.5 timingHeight]...
    ,'Style','text','String','Baseline (s)',...
    'FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,...
    'Unit','Normalized','BackgroundColor', backgroundColor,...
    'Position',[timingTextWidth 1-2*timingHeight...
    timingBoxWidth timingHeight],'Style','edit',...
    'String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,...
    'Unit','Normalized','BackgroundColor', backgroundColor,...
    'Position',[timingTextWidth+timingBoxWidth 1-2*timingHeight...
    timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),...
    'FontSize',fontSizeSmall);

% Stimulus Range
stimPeriod = [0.15 0.4];
uicontrol('Parent',hTimingPanel,'Unit','Normalized',...
    'Position',[0 1-3*timingHeight 0.5 timingHeight]...
    ,'Style','text','String','Stim Period (s)',...
    'FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,...
    'Unit','Normalized','BackgroundColor', backgroundColor,...
    'Position',[timingTextWidth 1-3*timingHeight...
    timingBoxWidth timingHeight],'Style','edit',...
    'String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,...
    'Unit','Normalized','BackgroundColor', backgroundColor,...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight...
    timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),...
    'FontSize',fontSizeSmall);  
            
            
hMTPanel = uipanel('Title','Multi-taper Panel','titleposition','centertop',...
                'fontSize',fontSizeLarge,'Unit','Normalized','Position',...
                [0.5 0.89 0.225 0.11]);
            
mtTextWidth = 0.5; mtBoxWidth = 0.2;
mtHeight = 1/3; 

            
% Time Bandwidth and Tapers
uicontrol('Parent',hMTPanel,'Unit','Normalized', ...
    'Position',[0 1-mtHeight mtTextWidth mtHeight], ...
    'Style','text','String','MT Parameter','FontSize',fontSizeSmall);

uicontrol('Parent',hMTPanel,'Unit','Normalized', ...
    'Position',[mtTextWidth 1-mtHeight mtBoxWidth mtHeight], ...
    'Style','text','String','TW','FontSize',fontSizeSmall);

uicontrol('Parent',hMTPanel,'Unit','Normalized', ...
    'Position',[mtTextWidth+mtBoxWidth 1-mtHeight mtBoxWidth mtHeight], ...
    'Style','text','String','K','FontSize',fontSizeSmall);

taperParameters = [1 1];
uicontrol('Parent',hMTPanel,'Unit','Normalized', ...
    'Position',[0 1-2*mtHeight mtTextWidth mtHeight], ...
    'Style','text','String','tapers','FontSize',fontSizeSmall);
hTimeBandWidthProduct = uicontrol('Parent',hMTPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[mtTextWidth 1-2*mtHeight mtBoxWidth mtHeight], ...
    'Style','edit','String',num2str(taperParameters(1)),'FontSize',fontSizeSmall);
hTapers = uicontrol('Parent',hMTPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[mtTextWidth+mtBoxWidth 1-2*mtHeight mtBoxWidth mtHeight], ...
    'Style','edit','String',num2str(taperParameters(2)),'FontSize',fontSizeSmall);
            
hLoadDataPanel = uipanel('Title','Data Processing Panel',...
                         'titleposition','centertop','fontSize',fontSizeLarge,...
                         'Unit','Normalized','Position',...
                [0.725 0.89 0.225 0.11]); 
            
hLFPResponsePanel = uibuttongroup('Parent',hLoadDataPanel,...
                    'Unit','Normalized','Position',[0 0.5 1 0.5]);
uicontrol('parent',hLFPResponsePanel,'Style','radiobutton', ...
          'String','Evoked Response','Unit','Normalized',...
          'Position',[0.1 0 0.5 1],'FontSize',fontSizeSmall);
              
uicontrol('parent',hLFPResponsePanel,'Style','radiobutton',...
  'String','Induced Response','Unit','Normalized',...
  'Position',[0.5 0 0.5 1],'FontSize',fontSizeSmall);

uicontrol('Parent',hLoadDataPanel,'Unit','Normalized',...
'Position',[0.25 0 0.5 0.5],'Style','pushbutton','String','Load DATA',...
'FontSize',fontSizeMedium,'Callback',{@processData_Callback});


    function processData_Callback(~,~)
        
        try
            % We turn the interface off for processing.
            InterfaceObj_DataProcessing = findobj(hFigure,'Enable','on');
            set(InterfaceObj_DataProcessing,'Enable','off');
            
            % get Session Information
            sessionNum = get(hSession,'val'); 
            fileNameStringTMP = fileNameStringListArray{sessionNum};
            blRange = [str2double(get(hBaselineMin,'String'))...
               str2double(get(hBaselineMax,'String'))];
            stRange = [str2double(get(hStimPeriodMin,'String'))...
               str2double(get(hStimPeriodMax,'String'))];

            gridType = 'Microelectrode';
            oriSelectiveFlag = get(hOriTunedCheckbox,'val');
            LFPdataProcessingMethod = get(get(hLFPResponsePanel,'SelectedObject'),'String');

            % parameters for MT analysis
            TW = str2double(get(hTimeBandWidthProduct,'String'));
            K = str2double(get(hTapers,'String'));
            if K~=2*TW-1
                error('Set TW and K as K = 2TW-1')
            elseif K == 2*TW-1
                tapers_MT = [TW K];
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%% Find Good Electrodes %%%%%%%%%%%%%%%%%%%
            [ElectrodeStringListAll,ElectrodeArrayListAll]= getElectrodesList(fileNameStringTMP,oriSelectiveFlag,folderSourceString);

            % Show electrodes on Grid
            electrodeGridPos = [0.05 panelStartHeight 0.22 panelHeight];
            hElectrodesonGrid = showElectrodeLocations(electrodeGridPos,[], ...
            [],[],1,0,gridType,'alpaH');% Electrode grid Layout are similar for both alpaH and kesariH hybrid grid


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freqRanges{1} = [8 12]; % alpha
            freqRanges{2} = [30 80]; % gamma
            freqRanges{3} = [104 250]; % hi-gamma
            freqRanges{4} = [16 16];  % SSVEP

            dataParameters.blRange = blRange;
            dataParameters.stRange = stRange;
            dataParameters.erpRange = [0.05 0.2];


           % get Data for Selected Session & Parameters
            [erpData,firingRateData,fftData,energyData,oriTuningData,~] = getData(folderSourceString,...
             fileNameStringTMP,ElectrodeArrayListAll,dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); 
%             freqRangeStr = {'alpha','gamma','SSVEP'};
%             numFreqRanges = length(freqRanges);       

            figure(1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters panel %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

            hParameterPanel = uipanel('Title','Parameters','titleposition','centertop','fontSize', ...
                fontSizeLarge,'Unit','Normalized','Position',...
                [0.275 panelStartHeight 0.225 panelHeight]);
            paramsHeight=1/6;

            % ElectrodeString
            uicontrol('Parent',hParameterPanel,'Unit','Normalized', ...
                'Position',[0 1-paramsHeight 0.5 paramsHeight],...
                'Style','text','String','Electrode','FontSize',fontSizeMedium);
            hElectrode = uicontrol('Parent',hParameterPanel,...
                'Unit','Normalized','BackgroundColor', backgroundColor, ...
                'Position',[0.5 1-paramsHeight 0.5 paramsHeight],...
                'Style','popup','String',ElectrodeStringListAll{1},...
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
            ['ERP|Firing Rate|Raster|Alpha [8-12 Hz]|Gamma Power [30-80 Hz]| Hi-Gamma Power [104-250 Hz]|'...
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
            axesLimitsTextWidth = 0.5; axesLimitsBoxWidth = 0.25;
            hAxesLimitsPanel = uipanel('Title','X and Y Limits','titleposition','centertop',...
                'fontSize', fontSizeLarge,'Unit','Normalized',...
                'Position',[0.5 panelStartHeight 0.225 panelHeight]);
            axesLimitsHeight = 1/6; 

            signalRange = [-0.1 0.5]; 
    %         erpPeriod = [0.05 0.2]; 
            fftRange = [0 100];

            % Signal Range
            uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized', ...
                'Position',[0 1-axesLimitsHeight axesLimitsTextWidth axesLimitsHeight],...
                'Style','text','String','Parameter','FontSize',fontSizeMedium);

            uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized', ...
                'Position',[axesLimitsTextWidth 1-axesLimitsHeight axesLimitsBoxWidth...
                axesLimitsHeight], 'Style','text','String','Min',...
                'FontSize',fontSizeMedium);

            uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized', ...
                'Position',[axesLimitsTextWidth+axesLimitsBoxWidth 1-axesLimitsHeight ...
                axesLimitsBoxWidth axesLimitsHeight], ...
                'Style','text','String','Max','FontSize',fontSizeMedium);

            % Stim Range
            uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized', ...
                'Position',[0 1-2*axesLimitsHeight axesLimitsTextWidth axesLimitsHeight]...
                ,'Style','text','String','Stim Range (s)',...
                'FontSize',fontSizeSmall);
            hStimMin = uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized', ...
                'BackgroundColor', backgroundColor, ...
                'Position',[axesLimitsTextWidth 1-2*axesLimitsHeight axesLimitsBoxWidth...
                axesLimitsHeight], ...
                'Style','edit','String',num2str(signalRange(1)),...
                'FontSize',fontSizeSmall);
            hStimMax = uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized',...
                'BackgroundColor', backgroundColor, ...
                'Position',[axesLimitsTextWidth+axesLimitsBoxWidth 1-2*axesLimitsHeight...
                axesLimitsBoxWidth axesLimitsHeight],'Style','edit',...
                'String',num2str(signalRange(2)),'FontSize',fontSizeSmall);


            % Frequency Range
            uicontrol('Parent',hAxesLimitsPanel,'Unit','Normalized',...
                'Position',...
                [0 1-3*axesLimitsHeight axesLimitsTextWidth axesLimitsHeight]...
                ,'Style','text','String','Frequency Range (Hz)',...
                'FontSize',fontSizeSmall);
            hFFTMin = uicontrol('Parent',hAxesLimitsPanel,...
                'Unit','Normalized', ...
                'BackgroundColor', backgroundColor, ...
                'Position',[axesLimitsTextWidth 1-3*axesLimitsHeight...
                axesLimitsBoxWidth axesLimitsHeight], ...
                'Style','edit','String',num2str(fftRange(1)),...
                'FontSize',fontSizeSmall);
            hFFTMax = uicontrol('Parent',hAxesLimitsPanel,...
                'Unit','Normalized',...
                'BackgroundColor', backgroundColor, ...
                'Position',...
                [axesLimitsTextWidth+axesLimitsBoxWidth 1-3*axesLimitsHeight...
                axesLimitsBoxWidth axesLimitsHeight],'Style','edit',...
                'String',num2str(fftRange(2)),'FontSize',fontSizeSmall);



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hPlotOptionsPanel = uipanel('Title','Plotting Options',...
                'titleposition','centertop','fontSize', fontSizeLarge,...
                'Unit','Normalized',...
                'Position',[0.725 panelStartHeight 0.225 panelHeight]);
            plotOptionsHeight = 1/6;

            % Button for Plotting
            [colorString, colorNames] = getColorString;
            uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
                'Position',...
                [0 5*plotOptionsHeight 0.6 plotOptionsHeight],...
                'Style','text','String','Color','FontSize',fontSizeSmall)
            hChooseColor = uicontrol('Parent',hPlotOptionsPanel,...
                'Unit','Normalized','BackgroundColor', backgroundColor,...
                'Position',...
                [0.6 5*plotOptionsHeight 0.4 plotOptionsHeight],...
                'Style','popup','String',colorString,...
                'FontSize',fontSizeSmall);

            uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
                'Position',...
                [0 3*plotOptionsHeight 0.5 plotOptionsHeight],...
                'Style','pushbutton','String','clear Plot',...
                'FontSize',fontSizeMedium,'Callback',{@cla_Callback});
            
            uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
                'Position',...
                [0.5 3*plotOptionsHeight 0.5 plotOptionsHeight],...
                'Style','pushbutton','String','Select New Session',...
                'BackgroundColor','w','ForegroundColor','r','FontSize',...
                fontSizeMedium,'Callback',{@ChooseNewSession_Callback});
            
            hHoldOn = uicontrol('Parent',hPlotOptionsPanel,...
                'Unit','Normalized','Position',[0 2*plotOptionsHeight ...
                1 plotOptionsHeight],'Style','togglebutton',...
                'String','hold on','FontSize',fontSizeMedium, ...
                'Callback',{@holdOn_Callback});

            uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
                'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
                'Style','pushbutton','String','Rescale',...
                'FontSize',fontSizeMedium,'Callback',{@rescaleXY_Callback})

            uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized',...
                'Position',[0 0 1 plotOptionsHeight],...
                'Style','pushbutton','String','plot',...
                'FontSize',fontSizeMedium,'Callback',{@plotData_Callback})

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            hFigure2 = figure(2); hFigure3 = figure(3);
            set(hFigure2,'units','normalized','outerposition',[0 0 1 1])
            set(hFigure3,'units','normalized','outerposition',[0 0 1 1])
            figure(2);
            hPlotPreferred = getPlotHandles(1,5,[0.1 0.7 0.5 0.15],0.01,0.01,1); linkaxes(hPlotPreferred);
            hOtherMeaures = getPlotHandles(1,3,[0.1 0.3 0.5 0.25],0.05,0.05,0);

            % Defining Text  in new axes for Figure 2
            textH3 = getPlotHandles(1,1,[0.2 0.9 0.01 0.01]); set(textH3,'Visible','Off');

            % more plots for AWP 
            figure(3);
            hEnergyPlots = getPlotHandles(1,2,[0.1 0.55 0.5 0.35],0.03,0.05,0); linkaxes(hEnergyPlots);
            hOtherPanels = getPlotHandles(1,3,[0.1 0.15 0.5 0.25],0.05,0.05,0);

        catch
            % We turn back on the interface
            set(InterfaceObj_DataProcessing,'Enable','on');
        end
        
        % Plotting Functions
        function plotData_Callback(~,~)
            
            try
                % We turn the interface off for processing.
                InterfaceObj_DataPlotting = findobj(hFigure,'Enable','on');
                set(InterfaceObj_DataPlotting,'Enable','off');

                %%%%%%%%%%%%%%%%%%%%%%%% Read values %%%%%%%%%%%%%%%%%%%%%%%%%%
                electrodeString = get(hElectrode,'val'); 
                if length(ElectrodeArrayListAll)==1
                    if isempty(ElectrodeArrayListAll{1}{electrodeString})
                        error('No electrode found for analysis! Pease try another session!')
                    elseif length(ElectrodeArrayListAll{1}{electrodeString})==1
                        electrodeNum = electrodeString;
                        disp('performing analysis on the selected single electrode');
                    elseif length(ElectrodeArrayListAll{1}{electrodeString})>1
                        disp('performing analysis on all electrodes for the selected single session')
                        electrodeNum = 'all';
                    end
                else
                   disp('performing analysis on all electrodes across sessions');
                   electrodeNum = 'all';
                end


                analysisMethod = get(hAnalysisMethod,'val');
                analysisMeasure = get(hAnalysisType,'val');
                NormalizeDataFlag = get(hNormalizeData,'val');
                relativeMeasuresFlag = get(hRelativeMeasures,'val');

                plotColor = colorNames(get(hChooseColor,'val'));
                holdOnState = get(hHoldOn,'val'); %#ok<NASGU>

    %             ColorNeuralMeasures = jet(5);

                hPlots_Fig1.hPlot1 = plotHandles;
                hPlots_Fig1.hPlot2 = hNeuralMeasureColorMatrix;
                hPlots_Fig1.hPlot3 = plotHandles2;
                hPlots_Fig1.hPlot4 = plotHandles3;
                hPlots_Fig1.hPlot5 = hRowCRF;
                hPlots_Fig1.hPlot6 = hColumnCRF;
                hPlots_Fig1.hPlot7 = hNormIndex;
                hPlots_Fig1.hPlot8 = hOriTuning;

                hPlots_Fig2.hPlot1 = hPlotPreferred;
                hPlots_Fig2.hPlot2 = hOtherMeaures; 

                hPlots_Fig3.hPlot1 = hEnergyPlots;
                hPlots_Fig3.hPlot2 = hOtherPanels;            


                if analysisMeasure == 1 % computing ERP
                    plotData(hPlots_Fig1,hPlots_Fig2,hPlots_Fig3,sessionNum,erpData.timeVals,erpData,oriTuningData,plotColor,analysisMethod,analysisMeasure,LFPdataProcessingMethod,relativeMeasuresFlag,NormalizeDataFlag,electrodeNum)
                elseif analysisMeasure == 2 % computing Firing rate
                    plotData(hPlots_Fig1,hPlots_Fig2,hPlots_Fig3,sessionNum,firingRateData.timeVals,firingRateData,oriTuningData,plotColor,analysisMethod,analysisMeasure,LFPdataProcessingMethod,relativeMeasuresFlag,NormalizeDataFlag,electrodeNum)
                elseif analysisMeasure == 3 % computing Raster Plot from spike data
                    error('Still working on raster data!')
                elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6|| analysisMeasure == 7% computing alpha
                    if analysisMethod == 1
                        plotData(hPlots_Fig1,hPlots_Fig2,hPlots_Fig3,sessionNum,fftData.freqVals,fftData,oriTuningData,plotColor,analysisMethod,analysisMeasure,LFPdataProcessingMethod,relativeMeasuresFlag,NormalizeDataFlag,electrodeNum)
                    elseif analysisMethod ==2 
                        plotData(hPlots_Fig1,hPlots_Fig2,hPlots_Fig3,sessionNum,energyData.freqVals,energyData,oriTuningData,plotColor,analysisMethod,analysisMeasure,LFPdataProcessingMethod,relativeMeasuresFlag,NormalizeDataFlag,electrodeNum)
                    end
                elseif analysisMeasure == 8 % need to work on STA!
                    error('STA computation method not found') 
                end

                figure(1);
                % Text for Orientation for main 5x5 Neural measure matrix
                textH1 = getPlotHandles(1,1,[0.2 0.65 0.01 0.01]);
                set(textH1,'Visible','Off');
                textH2 = getPlotHandles(1,1,[0.02 0.25 0.01 0.01]);
                set(textH2,'Visible','Off');

    %             if length(fileNameStringTMP) ==1 && length(ElectrodeListTMP{1})==1    
    %                 text(0.35,1.15,['Null Orientation: ' num2str(oValsUnique)],'unit','normalized','fontsize',20,'fontweight','bold','rotation',90,'parent',textH2);
    %                 text(0.35,1.15,['Preferred Orientation: ' num2str(oValsUnique2)],'unit','normalized','fontsize',20,'fontweight','bold','parent',textH1);
    %             else
    %                 text(0.35,1.15,'Null Orientation','unit','normalized','fontsize',20,'fontweight','bold','rotation',90,'parent',textH2);
    %                 text(0.35,1.15,'Preferred Orientation' ,'unit','normalized','fontsize',20,'fontweight','bold','parent',textH1);
    %             end

                % Setting Plot Ranges
                if analysisMeasure<=3 %ERP or spikes
                    xMin = str2double(get(hStimMin,'String'));
                    xMax = str2double(get(hStimMax,'String'));
                elseif analysisMeasure <=7  % LFP fft analysis
                    xMin = str2double(get(hFFTMin,'String'));
                    xMax = str2double(get(hFFTMax,'String'));                
                else
                    xMin = str2double(get(hSTAMin,'String'));
                    xMax = str2double(get(hSTAMax,'String'));
                end

                figure(2)
                % Text for Orientation for main 5x5 Neural measure matrix
                text(0.35,1.15,'Contrast along Orientation 1','unit','normalized','fontsize',20,'fontweight','bold','parent',textH3);

                rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
                rescaleData(plotHandles2,0,50,getYLims(plotHandles2));
                rescaleData(plotHandles3,0,50,getYLims(plotHandles3));
                rescaleData(hRowCRF,0,50,getYLims(hRowCRF));
                rescaleData(hColumnCRF,0,50,getYLims(hColumnCRF));
    %             rescaleData(hOriTuning,0,160,getYLims(hOriTuning));
                rescaleData(hPlotPreferred,xMin,xMax,getYLims(hPlotPreferred));
                rescaleData(hOtherMeaures(1),0,50,getYLims(hOtherMeaures(1)));
                rescaleData(hEnergyPlots,xMin,xMax,getYLims(hEnergyPlots));
                rescaleData(hOtherPanels(1),0,50,getYLims(hOtherPanels(1)));
                
                % We turn back on the Plotting interface
                set(InterfaceObj_DataPlotting,'Enable','on');
         
            catch e
            	% We turn back on the entire interface (plotting and data
            	% processing)
%                 fprintf(1,'The identifier was:\n%s',e.identifier);
%                 fprintf(1,'There was an error! The message was:\n%s',e.message);
                set(InterfaceObj_DataPlotting,'Enable','on');
                rethrow(e)
            end
        end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleXY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');

        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType==8
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
        
        claGivenPlotHandle(hEnergyPlots);
        claGivenPlotHandle(hOtherPanels);
        
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for i=1:numRows
                for j=1:numCols
                    cla(plotHandles(i,j));
                end
            end
        end

    end
    
    function ChooseNewSession_Callback(~,~)
        % We turn back on the interface
        set(InterfaceObj_DataProcessing,'Enable','on');
        delete(findobj(hElectrodesonGrid));
        delete(findobj(hParameterPanel));
        delete(findobj(hAxesLimitsPanel));
        delete(findobj(hPlotOptionsPanel));
        delete(findobj(plotHandles));
        delete(findobj(plotHandles2));
        delete(findobj(plotHandles3));
        delete(findobj(hOriTuning));
        delete(findobj(hNeuralMeasureColorMatrix));
        delete(findobj(hRowCRF));
        delete(findobj(hColumnCRF));
        delete(findobj(hNormIndex));
        close(hFigure2)
        close(hFigure3)
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,electrodeArray] = ...
    getData(folderSourceString,fileNameStringTMP,ElectrodeListTMP,dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod)

numDatasets = length(fileNameStringTMP);
disp(['Working on dataset 1 of ' num2str(numDatasets)]);
[erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,electrodeArray]...
= getDataSingleSession(folderSourceString,fileNameStringTMP{1},...
ElectrodeListTMP{1},dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod); 

if length(fileNameStringTMP)>1
    for i=2:numDatasets
        if isempty(ElectrodeListTMP{i}{end})
           continue
        end
        disp(['Working on dataset ' num2str(i) ' of ' num2str(length(fileNameStringTMP))]);
        [erpDataTMP,firingRateDataTMP,fftDataTMP,energyDataTMP,~,NI_DataTMP,~] = getDataSingleSession(folderSourceString,fileNameStringTMP{i},...
            ElectrodeListTMP{i},dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod);
        
        erpData.data = cat(1,erpData.data,erpDataTMP.data);
        erpData.analysisDataBL = cat(1,erpData.analysisDataBL,erpDataTMP.analysisDataBL);
        erpData.analysisData_cBL = cat(1,erpData.analysisData_cBL,erpDataTMP.analysisData_cBL);
        erpData.analysisDataST = cat(1,erpData.analysisDataST,erpDataTMP.analysisDataST);
        
        firingRateData.data = cat(1,firingRateData.data,firingRateDataTMP.data);
        firingRateData.analysisDataBL = cat(1,firingRateData.analysisDataBL,firingRateDataTMP.analysisDataBL);
        firingRateData.analysisData_cBL = cat(1,firingRateData.analysisData_cBL,firingRateDataTMP.analysisData_cBL);
        firingRateData.analysisDataST = cat(1,firingRateData.analysisDataST,firingRateDataTMP.analysisDataST);
        
        fftData.dataBL = cat(1,fftData.dataBL,fftDataTMP.dataBL);
        fftData.data_cBL = cat(1,fftData.data_cBL,fftDataTMP.data_cBL);
        fftData.dataST = cat(1,fftData.dataST,fftDataTMP.dataST);
        for j = 1:length(fftData.analysisDataBL)
            fftData.analysisDataBL{j} = cat(1,fftData.analysisDataBL{j},fftDataTMP.analysisDataBL{j});
            fftData.analysisData_cBL{j} = cat(1,fftData.analysisData_cBL{j},fftDataTMP.analysisData_cBL{j});
            fftData.analysisDataST{j} = cat(1,fftData.analysisDataST{j},fftDataTMP.analysisDataST{j});
        end
        
        energyData.dataBL = cat(1,energyData.dataBL,energyDataTMP.dataBL);
        energyData.data_cBL = cat(1,energyData.data_cBL,energyDataTMP.data_cBL);
        energyData.dataST = cat(1,energyData.dataST,energyDataTMP.dataST);
        for j =1:length(energyData.analysisDataBL)
            energyData.analysisDataBL{j} = cat(1,energyData.analysisDataBL{j},energyDataTMP.analysisDataBL{j});
            energyData.analysisData_cBL{j} = cat(1,energyData.analysisData_cBL{j},energyDataTMP.analysisData_cBL{j});
            energyData.analysisDataST{j} = cat(1,energyData.analysisDataST{j},energyDataTMP.analysisDataST{j});
        end
        
        NI_Data.erp = cat(1, NI_Data.erp, NI_DataTMP.erp);
        NI_Data.firingRate = cat(1, NI_Data.firingRate, NI_DataTMP.firingRate);
        for j =1:length(NI_Data.fft)
            NI_Data.fft{j} = cat(1, NI_Data.fft{j}, NI_DataTMP.fft{j});
            NI_Data.energy{j} = cat(1, NI_Data.energy{j}, NI_DataTMP.energy{j});
        end

        % Combining OriData across sessions need to be done! Right now, there is no requirement!        
        electrodeArray = [];
    end
end
end



function [erpData,firingRateData,fftData,energyData,oriTuningData,NI_Data,electrodeArray] = ...
    getDataSingleSession(folderSourceString,fileNameStringTMP,...
    ElectrodeListTMP,dataParameters,tapers_MT,freqRanges,oriSelectiveFlag,LFPdataProcessingMethod)

gridType = 'microelectrode';
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


folderName = fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,protocolName);
tuningProtocol_folderName =  fullfile(folderSourceString,'data',...
                        monkeyName,gridType,expDate,oriTuning_protocolName);                           

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');
folderSave = fullfile(folderSourceString,'Projects\PlaidNormalizationProject\savedData');
folderSave_oriTuning = fullfile(tuningProtocol_folderName,'savedData');

if ~exist(folderSave,'dir')
    mkdir(folderSave);
end
if ~exist(folderSave_oriTuning,'dir')
    mkdir(folderSave_oriTuning);
end

% Load Orientation Tuning dataFile for ori Tuning protocol for all elecs
oriTuningDataFile = fullfile(folderSave_oriTuning,['oriTuningData_' num2str(1000*dataParameters.stRange(1)) 'ms_' num2str(1000*dataParameters.stRange(2)) 'ms.mat']);
if exist(oriTuningDataFile,'file')
    disp(['Loading file ' oriTuningDataFile]);
    load(oriTuningDataFile);
else
    % Get OrientationTuning Data
    [computationVals,PO,OS] = savePrefOriAndOriSelectivitySpikes(monkeyName,expDate,oriTuning_protocolName,folderSourceString,gridType);
end

oriTuningData.PO = PO(ElectrodeListTMP{end});
oriTuningData.OS = OS(ElectrodeListTMP{end});
oriTuningData.FR = computationVals(ElectrodeListTMP{end},:);

if oriSelectiveFlag 
    fileToSave = fullfile(folderSave,[fileNameStringTMP '_OriTunedElecData_StimPeriod_' num2str(1000*dataParameters.stRange(1)) '_' num2str(1000*dataParameters.stRange(2)) 'ms_tapers' num2str(tapers_MT(1)) '_' num2str(tapers_MT(2)) '_' strtok(LFPdataProcessingMethod) '.mat']);
else
    fileToSave = fullfile(folderSave,[fileNameStringTMP '_allElecData_StimPeriod_' num2str(1000*dataParameters.stRange(1)) '_' num2str(1000*dataParameters.stRange(2)) 'ms_tapers' num2str(tapers_MT(1)) '_' num2str(tapers_MT(2)) '_' strtok(LFPdataProcessingMethod) '.mat']);
end

if exist(fileToSave,'file')
    disp(['Loading file ' fileToSave]);
    load(fileToSave)
else
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

    % electrode info
    ElectrodeList = ElectrodeListTMP{end};
    electrodeArray = ElectrodeList;

    % TimeVals info
    [~,timeVals,~,~] = loadlfpInfo(folderLFP);

    % Set up fft
    Fs = round(1/(timeVals(2)-timeVals(1)));
    range = dataParameters.blRange;
    rangePos = round(diff(range)*Fs);
    erpRangePos = round(diff(dataParameters.erpRange)*Fs);
    blPos = find(timeVals>=dataParameters.blRange(1),1)+ (1:rangePos);
    stPos = find(timeVals>=dataParameters.stRange(1),1)+ (1:rangePos);
    erpPos = find(timeVals>=dataParameters.erpRange(1),1)+ (1:erpRangePos);
    freqVals = 0:1/diff(range):Fs-1/diff(range);
    numFreqs = length(freqRanges);
    unitID = 0;

    % Set up params for MT
    params.tapers   = tapers_MT;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 250];
    params.trialave = 1;

    c1List = length(cValsUnique):-1:1; %  Flipping data Row-wise so that positive x-axis and positive y-axis denotes increase in Contrast
    c2List = 1:length(cValsUnique2);

    % Main Loop (Stores data in elec x tempFreq x Contrast of Ori 1 x
    % Contrast of Ori 2 x dataPoints)

    for iElec = 1:length(ElectrodeList)
        % Get LFP data
        clear analogData
        load(fullfile(folderLFP,['elec' num2str(ElectrodeList(iElec)) '.mat']));
        % Get Spike data
        clear spikeData
        load(fullfile(folderSpikes,['elec' num2str(ElectrodeList(iElec)) '_SID' num2str(unitID) '.mat']));

        % Get bad trials
        badTrialFile = fullfile(folderSegment,'badTrials.mat');
        if ~exist(badTrialFile,'file')
            disp('Bad trial file does not exist...');
            badTrials=[];
        else
            badTrials = loadBadTrials(badTrialFile);
            disp([num2str(length(badTrials)) ' bad trials']);
        end


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
                        if round(diff(dataParameters.blRange)*Fs) ~= round(diff(dataParameters.stRange)*Fs)
                            disp('baseline and stimulus ranges are not the same');
                        else
                            
                           % Event-related potential
                           erp = mean(analogData(goodPos,:),1); %#ok<NODEF>
                           erpDataTMP(iElec,t,c1,c2,:) = erp;
                           RMSvalsBL(iElec,t,c1,c2) = rms(erp(blPos));
                           RMSvalsERP(iElec,t,c1,c2) = rms(erp(erpPos));

                           % PSTH & firing rate 
                           [psthData(iElec,t,c1,c2,:),xsFR] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
                           firingRatesBL(iElec,t,c1,c2) = mean(getSpikeCounts(spikeData(goodPos),dataParameters.blRange))/diff(dataParameters.blRange);
                           firingRatesST(iElec,t,c1,c2) = mean(getSpikeCounts(spikeData(goodPos),dataParameters.stRange))/diff(dataParameters.stRange);

                           % fft data
                           if strcmp(strtok(LFPdataProcessingMethod),'Evoked')
                               fftBL = squeeze(mean(abs(fft(analogData(goodPos,blPos),[],2))));
                               fftST = squeeze(mean(abs(fft(analogData(goodPos,stPos),[],2))));
                           elseif strcmp(strtok(LFPdataProcessingMethod),'Induced')
                               fftBL = squeeze(mean(abs(fft(removeERP(analogData(goodPos,blPos)),[],2))));
                               fftST = squeeze(mean(abs(fft(removeERP(analogData(goodPos,stPos)),[],2))));
                           end
                           fftDataBL(iElec,t,c1,c2,:) = conv2Log(fftBL);
                           fftDataST(iElec,t,c1,c2,:) = conv2Log(fftST);
                           
                           % Power Estimation by MT method
                           if strcmp(strtok(LFPdataProcessingMethod),'Evoked')
                               dataBL = analogData(goodPos,blPos)';
                               dataST = analogData(goodPos,stPos)';
                           elseif strcmp(strtok(LFPdataProcessingMethod),'Induced')
                               dataBL = removeERP(analogData(goodPos,blPos))';
                               dataST = removeERP(analogData(goodPos,stPos))';                           
                           end
                           [tmpEBL,freqValsBL] = mtspectrumc(dataBL,params);
                           [tmpEST,freqValsST] = mtspectrumc(dataST,params);
                           
                           if isequal(freqValsBL,freqValsST)
                               freqValsMT = freqValsST;
                           end
                           mEnergyVsFreqBL(iElec,t,c1,c2,:) = conv2Log(tmpEBL);
                           mEnergyVsFreqST(iElec,t,c1,c2,:) = conv2Log(tmpEST);

                           % computing analysis Data for particular
                           % frequency band
                           for i=1:numFreqs
                               fftAmpBL{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(fftBL(:),freqVals,freqRanges{i}));
                               fftAmpST{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(fftST(:),freqVals,freqRanges{i}));
                               energyValsBL{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),freqValsMT,freqRanges{i}));
                               energyValsST{i}(iElec,t,c1,c2,:) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),freqValsMT,freqRanges{i}));
                           end
                        end
                    end
                end
            end
        end
    end
    
    % Time-Domain data
    erpData.data = erpDataTMP;
    erpData.analysisDataBL = RMSvalsBL;
    erpData.analysisData_cBL = getCommonBaseline(RMSvalsBL); % gets common (mean baseline) for all 5 x 5 contrast conditions 
    erpData.analysisDataST = RMSvalsERP;
    erpData.timeVals = timeVals;
    erpData.N = N;

    firingRateData.data = psthData;
    firingRateData.analysisDataBL = firingRatesBL;
    firingRateData.analysisData_cBL = getCommonBaseline(firingRatesBL); % gets common (mean baseline) for all 5 x 5 contrast conditions 
    firingRateData.analysisDataST = firingRatesST;
    firingRateData.timeVals = xsFR;
    firingRateData.N = N;

    % Freq-Domain data
    fftData.dataBL = fftDataBL;
    fftData.data_cBL = getCommonBaseline(fftDataBL);% gets common (mean baseline) for all 5 x 5 contrast conditions 
    fftData.dataST = fftDataST;
    fftData.analysisDataBL = fftAmpBL;
    fftData.analysisData_cBL = getCommonBaseline(fftAmpBL);% gets fft amplitude of common (mean baseline) for all 5 x 5 contrast conditions for selected freqRanges
    fftData.analysisDataST = fftAmpST;
    fftData.freqVals = freqVals;
    fftData.N = N;

    energyData.dataBL = mEnergyVsFreqBL;
    energyData.data_cBL = getCommonBaseline(mEnergyVsFreqBL);% gets common (mean baseline) for all 5 x 5 contrast conditions 
    energyData.dataST = mEnergyVsFreqST;
    energyData.analysisDataBL = energyValsBL;
    energyData.analysisData_cBL = getCommonBaseline(energyValsBL);% gets energyData of  common (mean baseline) for all 5 x 5 contrast conditions for selected freqRanges
    energyData.analysisDataST = energyValsST;
    energyData.freqVals = freqValsMT;
    energyData.N = N;

    % Segregation into Preferred-null axis is done when analysis is being
    % done for orientation selective electrodes
    if oriSelectiveFlag
        elecs_neededtoFlipped = find(abs(oriTuningData.PO-oValsUnique)<abs(oriTuningData.PO-oValsUnique2));
        erpData = segregate_Pref_Null_data(erpData,elecs_neededtoFlipped);
        firingRateData = segregate_Pref_Null_data(firingRateData,elecs_neededtoFlipped);
        fftData = segregate_Pref_Null_data(fftData,elecs_neededtoFlipped);
        energyData = segregate_Pref_Null_data(energyData,elecs_neededtoFlipped);
    end

    
    % Get Normalization Indices
    NI_Data.erp = getNI(erpData);
    NI_Data.firingRate = getNI(fftData);
    NI_Data.fft = getNI(fftData);
    NI_Data.energy = getNI(energyData);

    % Save Data for particular session
    save(fileToSave,'erpData','firingRateData','fftData','energyData','oriTuningData','NI_Data','electrodeArray');
end
end

function plotData(hPlots_Fig1,hPlots_Fig2,hPlots_Fig3,sessionNum,xs,data,oriTuningData,colorName,analysisMethod,analysisMeasure,LFPdataProcessingMethod,relativeMeasuresFlag,NormalizeDataFlag,electrodeNum)

if strcmp(strtok(LFPdataProcessingMethod),'Induced') && analysisMeasure == 7
    error('SSVEP is a time locked Signal. Subtracting ERP from each trials will wipe out the SSVEP signal at 2nd Harmonic. Please use evoked SSVEP response')
end

cValsUnique = [0 12.5 25 50 100]./2;
cValsUnique2 = [0 12.5 25 50 100]./2;

if isnumeric(electrodeNum)
    data = getDataSingleElec(data,electrodeNum,analysisMeasure);
elseif strcmp(electrodeNum,'all')
end

% Main 5x5 plot for Neural Measure
if analysisMeasure == 1 || analysisMeasure == 2
    if NormalizeDataFlag
    % Normalize ERP and spike Data (fft or energy data need not be normalized
    % as they are expressed in log units)
    data = normalizeData(data);
    end
    dataPlot = squeeze(mean(data.data(:,1,:,:,:),1)); % Handles single electrode or multi-electrode data similarly
    
elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
    if size(data.dataBL) ~= size(data.dataST) 
       error('Size of fftDataBL and fftDataST do not match!')
    end
    
    % Absolute BL & ST neural measure
    dataPlotBL = squeeze(mean(data.data_cBL(:,1,:,:,:),1));
    dataPlotST = squeeze(mean(data.dataST(:,1,:,:,:),1));
    if analysisMeasure == 7
        dataPlotBL = squeeze(mean(data.data_cBL(:,2,:,:,:),1));
        dataPlotST = squeeze(mean(data.dataST(:,2,:,:,:),1));
    end

    % When Change in neural measures are to be plotted
    if relativeMeasuresFlag
        dataPlotdiffSTvsBL = dataPlotST-dataPlotBL;
        if analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
            if analysisMethod == 2
                dataPlotdiffSTvsBL = 10*(dataPlotST-dataPlotBL);  % Change in Power expressed in deciBel for MT method
            end
        end
    end
end
    
% Plotting 5x5 plots for raw Neural Measures (Figure 1)
for c1 = 1:5
    for c2 = 1:5
        if analysisMeasure == 1 || analysisMeasure == 2
            plot(hPlots_Fig1.hPlot1(c1,c2),xs,squeeze(dataPlot(c1,c2,:)),'color',colorName);
        elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
            if ~relativeMeasuresFlag
            plot(hPlots_Fig1.hPlot1(c1,c2),xs,squeeze(dataPlotBL(c1,c2,:)),'g');
            hold(hPlots_Fig1.hPlot1(c1,c2),'on')
            plot(hPlots_Fig1.hPlot1(c1,c2),xs,squeeze(dataPlotST(c1,c2,:)),'k');
            hold(hPlots_Fig1.hPlot1(c1,c2),'off')
            else
                plot(hPlots_Fig1.hPlot1(c1,c2),xs,squeeze(dataPlotdiffSTvsBL(c1,c2,:)),'b');
            end
        end
    end
end

conFlipped = 5:-1:1;
for c = 1:5
title(hPlots_Fig1.hPlot1(1,c),[num2str(cValsUnique(c)) ' %']);
ylabel(hPlots_Fig1.hPlot1(c,1),[num2str(cValsUnique(conFlipped(c))) ' %'],'fontWeight','bold');
end

colors = jet(5);
% Plot Figure 2 neural response plots
for c1 = 1:5
    for c2 = 1:5
        if analysisMeasure == 1 || analysisMeasure == 2
            plot(hPlots_Fig2.hPlot1(1,c2),xs,squeeze(dataPlot(conFlipped(c1),c2,:)),'color',colors(c1,:,:),'LineWidth',2);
            hold(hPlots_Fig2.hPlot1(1,c2),'on');
        elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
            if ~relativeMeasuresFlag
                plot(hPlots_Fig2.hPlot1(1,c2),xs,squeeze(dataPlotBL(conFlipped(c1),c2,:)),'k');
                hold(hPlots_Fig2.hPlot1(1,c2),'on')
                plot(hPlots_Fig2.hPlot1(1,c2),xs,squeeze(dataPlotST(conFlipped(c1),c2,:)),'color',colors(c1,:,:),'LineWidth',2);
            else
                plot(hPlots_Fig2.hPlot1(1,c2),xs,squeeze(dataPlotBL(conFlipped(c1),c2,:))-squeeze(dataPlotBL(conFlipped(c1),c2,:)),'k');
                hold(hPlots_Fig2.hPlot1(1,c2),'on')
                plot(hPlots_Fig2.hPlot1(1,c2),xs,squeeze(dataPlotdiffSTvsBL(conFlipped(c1),c2,:)),'color',colors(c1,:,:),'LineWidth',2);
            end
        end
        title(hPlots_Fig2.hPlot1(1,c2),[num2str(cValsUnique(c2)) ' %'])
    end
end

Norm_colorNames(1,:) = colors(5,:);
Norm_colorNames(2:5,:) = flip(parula(4),1);
% Plot Figure 3 neural response plots
for c = 1:5
    if analysisMeasure == 1 || analysisMeasure == 2
        plot(hPlots_Fig3.hPlot1(1),xs,squeeze(dataPlot(conFlipped(1),c,:)),'color',colors(c,:,:),'LineWidth',2);
        hold(hPlots_Fig3.hPlot1(1),'on');
        plot(hPlots_Fig3.hPlot1(2),xs,squeeze(dataPlot(c,conFlipped(1),:)),'color',Norm_colorNames(conFlipped(c),:,:),'LineWidth',2);
        hold(hPlots_Fig3.hPlot1(2),'on');
    elseif analysisMeasure == 4 || analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
        if ~relativeMeasuresFlag
        plot(hPlots_Fig3.hPlot1(1),xs,squeeze(dataPlotBL(conFlipped(1),c,:)),'k');
        hold(hPlots_Fig3.hPlot1(1),'on')
        plot(hPlots_Fig3.hPlot1(1),xs,squeeze(dataPlotST(conFlipped(1),c,:)),'color',colors(c,:,:),'LineWidth',2);
        
        plot(hPlots_Fig3.hPlot1(2),xs,squeeze(dataPlotBL(conFlipped(c),conFlipped(1),:)),'k');
        hold(hPlots_Fig3.hPlot1(2),'on')
        plot(hPlots_Fig3.hPlot1(2),xs,squeeze(dataPlotST(conFlipped(c),conFlipped(1),:)),'color',Norm_colorNames(c,:,:),'LineWidth',2);      
        else
            plot(hPlots_Fig3.hPlot1(1),xs,squeeze(dataPlotBL(conFlipped(1),c,:))-squeeze(dataPlotBL(conFlipped(1),c,:)),'k');
            hold(hPlots_Fig3.hPlot1(1),'on')
            plot(hPlots_Fig3.hPlot1(1),xs,squeeze(dataPlotdiffSTvsBL(conFlipped(1),c,:)),'color',colors(c,:,:),'LineWidth',2);
            
            plot(hPlots_Fig3.hPlot1(2),xs,squeeze(dataPlotBL(conFlipped(c),conFlipped(1),:))-squeeze(dataPlotBL(conFlipped(c),conFlipped(1),:)),'k');
            hold(hPlots_Fig3.hPlot1(2),'on')
            if analysisMeasure == 4||analysisMeasure ==5||analysisMeasure ==6
                plot(hPlots_Fig3.hPlot1(2),xs,squeeze(dataPlotdiffSTvsBL(conFlipped(c),conFlipped(1),:)),'color',Norm_colorNames(c,:,:),'LineWidth',2);
            elseif analysisMeasure == 7
                plot(hPlots_Fig3.hPlot1(2),xs,squeeze(dataPlotdiffSTvsBL(c,conFlipped(1),:)),'color',Norm_colorNames(conFlipped(c),:,:),'LineWidth',2);
            end
        end
    end
    
    if analysisMeasure == 2
        text(0.55,0.45+c*0.08,[num2str(cValsUnique(c)) ' %'],'color',colors(c,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(1))
        text(0.55,0.45+c*0.08,[num2str(cValsUnique(c)) ' %'],'color',Norm_colorNames(c,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(2))
        text(0.55,0.45+6*0.08,'Ori:1','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(1))
        text(0.55,0.45+6*0.08,'Ori:2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(2))

    elseif analysisMeasure == 5
        text(0.05,0.45+c*0.08,[num2str(cValsUnique(c)) ' %'],'color',colors(c,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(1))
        text(0.05,0.45+6*0.08,'Ori:1','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(1))

        text(0.05,0.45+c*0.08,[num2str(cValsUnique(conFlipped(c))) ' %'],'color',Norm_colorNames(conFlipped(c),:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(2))
        text(0.05,0.45+6*0.08,'Ori:2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(2))
    elseif analysisMeasure == 6||analysisMeasure == 7
        text(0.1+c*0.15,0.15,[num2str(cValsUnique(c)) ' %'],'color',colors(c,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(1))
        if c==1
            text(0.1,0.15,'Ori 1:','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(1))
        end

        text(0.1+c*0.15,0.15,[num2str(cValsUnique(c)) ' %'],'color',Norm_colorNames(c,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(2))
        if c==1
            text(0.1,0.15,'Ori 2:','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot1(2))
        end
    end

end
if analysisMeasure == 2
    title(hPlots_Fig3.hPlot1(1),'Spike Response for Orientation 1')
    title(hPlots_Fig3.hPlot1(2),{'Spike Response at 50% contrast of Ori 1','when Ori 2 is added'},'HorizontalAlignment','Center')
else
    title(hPlots_Fig3.hPlot1(1),'Change in PSD for Orientation 1')
    title(hPlots_Fig3.hPlot1(2),{'Change in PSD at 50% contrast of Ori 1','when Ori 2 is added'},'HorizontalAlignment','Center')
end
hold(hPlots_Fig3.hPlot1(1),'off');
hold(hPlots_Fig3.hPlot1(2),'off');
% Figure 1
% Color coded 5x5 Matrix of raw Neural Measures
if analysisMeasure == 1 || analysisMeasure == 2
    analysisData = squeeze(mean(data.analysisDataST(:,1,:,:),1));
    sem_analysisData = squeeze(std(squeeze(data.analysisDataST(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST,1)));
    analysisDataBL = squeeze(mean(data.analysisData_cBL(:,1,:,:),1));
    for iElec= 1:size(data.analysisDataST,1)
        clear electrodeVals
        if ~relativeMeasuresFlag
            electrodeVals =  squeeze(data.analysisDataST(iElec,1,:,:));
        else
            electrodeVals =  squeeze(data.analysisDataST(iElec,1,:,:))-squeeze(data.analysisDataBL(iElec,1,:,:));
        end
        NI_population(iElec) = electrodeVals(1,5)/(((electrodeVals(1,1)+electrodeVals(5,5)))/2)-1;
    end
elseif analysisMeasure == 4 
    analysisData = squeeze(mean(data.analysisDataST{1}(:,1,:,:),1));
    sem_analysisData = squeeze(std(squeeze(data.analysisDataST{1}(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST{1},1)));
    analysisDataBL = squeeze(mean(data.analysisData_cBL{1}(:,1,:,:),1));
    for iElec= 1:size(data.analysisDataST{1},1)
        clear electrodeVals
        if ~relativeMeasuresFlag
            electrodeVals =  squeeze(data.analysisDataST{1}(iElec,1,:,:));
        else
            electrodeVals =  squeeze(data.analysisDataST{1}(iElec,1,:,:))-squeeze(data.analysisDataBL{1}(iElec,1,:,:));
        end
        NI_population(iElec) = electrodeVals(1,5)/(((electrodeVals(1,1)+electrodeVals(5,5)))/2)-1;
    end        
elseif analysisMeasure == 5
    analysisData = squeeze(mean(data.analysisDataST{2}(:,1,:,:),1));
    sem_analysisData = squeeze(std(squeeze(data.analysisDataST{2}(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST{2},1)));
    analysisDataBL = squeeze(mean(data.analysisData_cBL{2}(:,1,:,:),1));
    for iElec= 1:size(data.analysisDataST{2},1)
        clear electrodeVals
        if ~relativeMeasuresFlag
            electrodeVals =  squeeze(data.analysisDataST{2}(iElec,1,:,:));
        else
            electrodeVals =  squeeze(data.analysisDataST{2}(iElec,1,:,:))-squeeze(data.analysisDataBL{2}(iElec,1,:,:));
        end
        NI_population(iElec) = electrodeVals(1,5)/(((electrodeVals(1,1)+electrodeVals(5,5)))/2)-1;
    end
elseif analysisMeasure == 6
    analysisData = squeeze(mean(data.analysisDataST{3}(:,1,:,:),1));
    sem_analysisData = squeeze(std(squeeze(data.analysisDataST{3}(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST{3},1)));
    analysisDataBL = squeeze(mean(data.analysisData_cBL{3}(:,1,:,:),1));
    for iElec= 1:size(data.analysisDataST{3},1)
        clear electrodeVals
        if ~relativeMeasuresFlag
            electrodeVals =  squeeze(data.analysisDataST{3}(iElec,1,:,:));
            NI_population(iElec) = abs(electrodeVals(1,5)/(((electrodeVals(1,1)+electrodeVals(5,5)))/2)-1); % log of low power values are negative!
        else
            electrodeVals =  squeeze(data.analysisDataST{3}(iElec,1,:,:))-squeeze(data.analysisDataBL{3}(iElec,1,:,:));
            NI_population(iElec) = electrodeVals(1,5)/(((electrodeVals(1,1)+electrodeVals(5,5)))/2)-1;
        end
    end 
elseif analysisMeasure == 7
    analysisData = squeeze(mean(data.analysisDataST{4}(:,2,:,:),1));
    sem_analysisData = squeeze(std(squeeze(data.analysisDataST{4}(:,2,:,:)),[],1)./sqrt(size(data.analysisDataST{4},1)));
    analysisDataBL = squeeze(mean(data.analysisData_cBL{4}(:,2,:,:),1));
    for iElec= 1:size(data.analysisDataST{4},1)
        clear electrodeVals
        if ~relativeMeasuresFlag
            electrodeVals =  squeeze(data.analysisDataST{4}(iElec,2,:,:));
        else
            electrodeVals =  squeeze(data.analysisDataST{4}(iElec,2,:,:))-squeeze(data.analysisDataBL{4}(iElec,2,:,:));
        end
        NI_population(iElec) = electrodeVals(1,5)/(((electrodeVals(1,1)+electrodeVals(5,5)))/2)-1;
    end
end

if relativeMeasuresFlag
    if analysisMeasure == 1||analysisMeasure == 2
            analysisData = analysisData-analysisDataBL;
            sem_analysisData =  squeeze(std(squeeze(data.analysisDataST(:,1,:,:))-squeeze(data.analysisDataBL(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST,1)));
    elseif analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
        if analysisMethod == 2
            analysisData = 10*(analysisData-analysisDataBL);% Change in power expressed in deciBel
        else
            analysisData = (analysisData-analysisDataBL);
        end
        
        if analysisMeasure == 4
            sem_analysisData = squeeze(std(squeeze(data.analysisDataST{1}(:,1,:,:))-squeeze(data.analysisDataBL{1}(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST{1},1)));
        elseif analysisMeasure == 5
            sem_analysisData = squeeze(std(squeeze(data.analysisDataST{2}(:,1,:,:))-squeeze(data.analysisDataBL{2}(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST{2},1)));
        elseif analysisMeasure == 6
            sem_analysisData = squeeze(std(squeeze(data.analysisDataST{3}(:,1,:,:))-squeeze(data.analysisDataBL{3}(:,1,:,:)),[],1)./sqrt(size(data.analysisDataST{3},1)));
        elseif analysisMeasure == 7
            sem_analysisData = squeeze(std(squeeze(data.analysisDataST{4}(:,2,:,:))-squeeze(data.analysisDataBL{4}(:,2,:,:)),[],1)./sqrt(size(data.analysisDataST{4},1)));
        end
    end
end

% Figure 1
imagesc(analysisData,'parent',hPlots_Fig1.hPlot2);colorbar(hPlots_Fig1.hPlot2);set(hPlots_Fig1.hPlot2,'Position',[0.52 0.05 0.12 0.25]);
title(hPlots_Fig1.hPlot2,['Mean NI: ',num2str(round(mean(NI_population),2))],'fontWeight','bold');

% Setting xlabels and ylabels for Color matrix of neural measures
set(hPlots_Fig1.hPlot2,'XTickLabel',cValsUnique);
set(hPlots_Fig1.hPlot2,'YTickLabel',flip(cValsUnique2));

% Ori Tuning data shown only for single session (single or all electrodes)
if analysisMeasure == 1||analysisMeasure == 2
    num_Electrodes = size(data.data,1);
else
    num_Electrodes = size(data.dataST,1);    
end
colorNamesOriTuning = hsv(num_Electrodes);
oValsUnique_Tuning = [0 22.5 45 67.5 90 112.5 135 157.5];

if sessionNum <=22 % Single Session
    if num_Electrodes == 1
        plot(hPlots_Fig1.hPlot8,oValsUnique_Tuning,oriTuningData.FR(electrodeNum,:),'Marker','o','color',colorNamesOriTuning);
        text(0.7,0.5,['PO: ' num2str(oriTuningData.PO(electrodeNum)) ',OS: ' num2str(round(oriTuningData.OS(electrodeNum),2))],'color',colorNamesOriTuning,'unit','normalized','fontSize',6,'parent',hPlots_Fig1.hPlot8);
    else
        for iElec = 1:num_Electrodes
            plot(hPlots_Fig1.hPlot8,oValsUnique_Tuning,oriTuningData.FR(iElec,:),'Marker','o','color',colorNamesOriTuning(iElec,:,:));
            text(0.7,iElec*0.05+0.1,['PO: ' num2str(oriTuningData.PO(iElec)) ',OS: ' num2str(round(oriTuningData.OS(iElec),2))],'color',colorNamesOriTuning(iElec,:,:),'unit','normalized','fontSize',6,'parent',hPlots_Fig1.hPlot8);
            hold(hPlots_Fig1.hPlot8,'on');
        end
        hold(hPlots_Fig1.hPlot8,'off');
    end
else
    text(0.5,0.5,{'Orientation Tuning data' 'not shown for' 'multiple sessions'},'unit','normalized','HorizontalAlignment','center','fontSize',10,'parent',hPlots_Fig1.hPlot8);
end

% Setting title,xlabel and ylabel for OriTuning Plot
title(hPlots_Fig1.hPlot8,'Ori Tuning for Single Session (Spike Data)','fontSize',10);
xlim(hPlots_Fig1.hPlot8,[0 250]);
set(hPlots_Fig1.hPlot8,'XTick',oValsUnique_Tuning); set(hPlots_Fig1.hPlot8,'XTickLabelRotation',90);

% NI population histogram
if num_Electrodes>1
histogram(hPlots_Fig1.hPlot7,NI_population,-2:0.2:2);
end

% Figure 2
imagesc(analysisData,'parent',hPlots_Fig2.hPlot2(3));
color_Bar = colorbar(hPlots_Fig2.hPlot2(3)); 
colorYlabelHandle = get(color_Bar,'Ylabel');
if analysisMeasure==1
    YlabelString = 'Potential';
elseif analysisMeasure==2
    YlabelString = 'Spikes/s';
elseif analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6||analysisMeasure == 7
    if ~relativeMeasuresFlag
        YlabelString = 'log_1_0(FFT Amplitude)';
    else
        YlabelString = 'log_1_0(\Delta FFT Amplitude)';
        if analysisMethod == 2
            YlabelString = 'Change in Power (dB)';
        end
    end
end

set(colorYlabelHandle,'String',YlabelString,'fontSize',14);
plotPos = get(hPlots_Fig2.hPlot2(3),'Position');
set(hPlots_Fig2.hPlot2(3),'Position',[plotPos(1) plotPos(2) 0.12 plotPos(4)]);
title(hPlots_Fig2.hPlot2(3),['Mean NI: ',num2str(round(mean(NI_population),2))],'fontWeight','bold');

if num_Electrodes>1
% NI population histogram
histogram(hPlots_Fig2.hPlot2(2),NI_population,-2:0.2:2);
end
tickLengthPlot = 2*get(hPlots_Fig2.hPlot2(1),'TickLength');
set(hPlots_Fig2.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
xlabel(hPlots_Fig2.hPlot2(2),'Normalization Index'); ylabel(hPlots_Fig2.hPlot2(2),'No. of Electrodes');
title(hPlots_Fig2.hPlot2(2),'Population NI Histogram')

set(hPlots_Fig2.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlots_Fig2.hPlot2(3),'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlots_Fig2.hPlot2(3),'Contrast (%)');ylabel(hPlots_Fig2.hPlot2(3),'Contrast (%)');

% Figure 3
imagesc(analysisData,'parent',hPlots_Fig3.hPlot2(3));
% grid on;
color_Bar = colorbar(hPlots_Fig3.hPlot2(3));% 
colorYlabelHandle = get(color_Bar,'Ylabel');
if analysisMeasure==1
    YlabelString = 'Potential';
elseif analysisMeasure==2
    YlabelString = 'Spikes/s';
elseif analysisMeasure == 4||analysisMeasure == 5||analysisMeasure == 6
    if ~relativeMeasuresFlag
        YlabelString = 'log_1_0(FFT Amplitude)';
    else
        YlabelString = 'log_1_0(\Delta FFT Amplitude)';
        if analysisMethod == 2
            YlabelString = 'Change in Power (dB)';
        end
    end
end

set(colorYlabelHandle,'String',YlabelString,'fontSize',14);
plotPos = get(hPlots_Fig3.hPlot2(3),'Position');
set(hPlots_Fig3.hPlot2(3),'Position',[plotPos(1) plotPos(2) 0.12 plotPos(4)]);
title(hPlots_Fig3.hPlot2(3),['Mean NI: ',num2str(round(mean(NI_population),2))],'fontWeight','bold');

if num_Electrodes>1
% NI population histogram
histogram(hPlots_Fig3.hPlot2(2),NI_population,-2:0.2:2);
end

set(hPlots_Fig3.hPlot2(2),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
xlabel(hPlots_Fig3.hPlot2(2),'Normalization Index'); ylabel(hPlots_Fig3.hPlot2(2),'No. of Electrodes');
title(hPlots_Fig3.hPlot2(2),{'Population NI Histogram,' ,['Median NI: ' num2str(round(median(NI_population),2))]})
line([round(median(NI_population),2) round(median(NI_population),2)],getYLims(hPlots_Fig3.hPlot2(2)),'color','k','lineWidth',2,'parent',hPlots_Fig3.hPlot2(2))

set(hPlots_Fig3.hPlot2(3),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlots_Fig3.hPlot2(3),'XTick',1:5,'XTickLabelRotation',90,'XTickLabel',cValsUnique,'YTickLabel',flip(cValsUnique));
xlabel(hPlots_Fig3.hPlot2(3),'Ori 1 Contrast (%)');ylabel(hPlots_Fig3.hPlot2(3),'Ori 2 Contrast (%)');

% Figure 1
% Contrast Response curves Row-Wise & Column-wise
CRFColors = jet(length(cValsUnique));
for iCon = 1:5
    plot(hPlots_Fig1.hPlot3(1,iCon),cValsUnique,analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    plot(hPlots_Fig1.hPlot4(1,iCon),cValsUnique2,flip(analysisData(:,iCon),1),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    plot(hPlots_Fig1.hPlot5(1,1),cValsUnique,analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlots_Fig1.hPlot5(1,1),'on')
    plot(hPlots_Fig1.hPlot6(1,1),cValsUnique2,flip(analysisData(:,iCon),1),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlots_Fig1.hPlot6(1,1),'on')
end
hold(hPlots_Fig1.hPlot5(1,1),'off'); hold(hPlots_Fig1.hPlot6(1,1),'off')

% Figure 2
CRFColors = jet(length(cValsUnique));
for iCon = 1:5
    plot(hPlots_Fig2.hPlot2(1),cValsUnique,analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
    hold(hPlots_Fig2.hPlot2(1),'on')
end
hold(hPlots_Fig2.hPlot2(1),'off');

if analysisMeasure == 1||analysisMeasure == 2
    displayRange(hPlots_Fig2.hPlot1,[0.2 0.4],getYLims(hPlots_Fig2.hPlot1),'k');
elseif analysisMeasure ==4 
    displayRange(hPlots_Fig2.hPlot1,[8 12],getYLims(hPlots_Fig2.hPlot1),'k');
elseif analysisMeasure == 5  
    displayRange(hPlots_Fig2.hPlot1,[30 80],getYLims(hPlots_Fig2.hPlot1),'k');
elseif analysisMeasure == 6     
    displayRange(hPlots_Fig2.hPlot1,[16 16],getYLims(hPlots_Fig2.hPlot1),'k');
end

tickLengthPlot = 2*get(hPlots_Fig2.hPlot2(1),'TickLength');

for idx =1:5
set(hPlots_Fig2.hPlot1(idx),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
if analysisMeasure == 1||analysisMeasure == 2
    xlabel(hPlots_Fig2.hPlot1(1),'Time(ms)');
elseif analysisMeasure == 4||analysisMeasure == 5|| analysisMeasure == 6 ||analysisMeasure == 7
    xlabel(hPlots_Fig2.hPlot1(1),'Frequency(Hz)');
end

if analysisMeasure == 1
    ylabel(hPlots_Fig2.hPlot1(1),'ERP (\mu V)');ylabel(hPlots_Fig2.hPlot2(1),'RMS value of ERP');
elseif analysisMeasure == 2 
    ylabel(hPlots_Fig2.hPlot1(1),'Spikes/s');ylabel(hPlots_Fig2.hPlot2(1),'spikes/s');
elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6||analysisMeasure == 7
    if analysisMethod ==1
        ylabel(hPlots_Fig2.hPlot1(1),'log_1_0 (FFT Amplitude)'); ylabel(hPlots_Fig2.hPlot2(1),'log_1_0 (FFT Amplitude)');
        if relativeMeasuresFlag
            ylabel(hPlots_Fig2.hPlot1(1),'log_1_0 (\Delta FFT Amplitude)'); ylabel(hPlots_Fig2.hPlot2(1),'log_1_0 (\Delta FFT Amplitude)')
        end

    elseif analysisMethod ==2
        ylabel(hPlots_Fig2.hPlot1(1),'log_1_0 (Power)');ylabel(hPlots_Fig2.hPlot2(1),'log_1_0 (Power)')
        if relativeMeasuresFlag
            ylabel(hPlots_Fig2.hPlot1(1),'Change in Power (dB)');  ylabel(hPlots_Fig2.hPlot2(1),'Change in Power (dB)');
        end
    end

end

text(0.5,0.3,['N = ' num2str(num_Electrodes)],'color','k','unit','normalized','fontSize',14,'fontWeight','bold','parent',hPlots_Fig2.hPlot2(1))
set(hPlots_Fig2.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlots_Fig2.hPlot2(1),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlots_Fig2.hPlot2(1),'Contrast (%)');
title(hPlots_Fig2.hPlot2(1),'CRF along Orientation 1')

% Figure 3
% errorbar(cValsUnique,analysisData(conFlipped(1),:),sem_analysisData(conFlipped(1),:),...
%     'Marker','o','LineWidth',2,'color',Norm_colorNames(1,:,:),'parent',hPlots_Fig3.hPlot2(1))
%     hold(hPlots_Fig3.hPlot2(1),'on');
%     
% errorbar(cValsUnique,diag(flipud(analysisData)),diag(flipud(sem_analysisData)),'Marker','o','LineWidth',2,'color','k','parent',hPlots_Fig3.hPlot2(1));
%     
for iCon = 1:5
    errorbar(cValsUnique,analysisData(conFlipped(iCon),:),sem_analysisData(conFlipped(iCon),:),...
        'Marker','o','LineWidth',2,'color',Norm_colorNames(iCon,:,:),'parent',hPlots_Fig3.hPlot2(1))
    hold(hPlots_Fig3.hPlot2(1),'on')
    text(0.5,iCon*0.07+0.05,['Ori 2: ' num2str(cValsUnique(iCon)) ' %'],'color',Norm_colorNames(iCon,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot2(1))
end
hold(hPlots_Fig3.hPlot2(1),'off');
% text(0.5,0.2,'Ori 2: 0%','color',colors(5,:,:),'fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot2(1))
% text(0.5,0.1,'Ori 1 = Ori 2','color','k','fontWeight','bold','fontSize',14,'unit','normalized','parent',hPlots_Fig3.hPlot2(1))


if analysisMeasure == 1||analysisMeasure == 2
    displayRange(hPlots_Fig3.hPlot1,[0.2 0.4],getYLims(hPlots_Fig3.hPlot1),'k');
elseif analysisMeasure ==4 
    displayRange(hPlots_Fig3.hPlot1,[108 248],getYLims(hPlots_Fig3.hPlot1),'k');
elseif analysisMeasure == 5  
    displayRange(hPlots_Fig3.hPlot1,[30 80],getYLims(hPlots_Fig3.hPlot1),'k');
elseif analysisMeasure == 6     
    displayRange(hPlots_Fig3.hPlot1,[104 250],getYLims(hPlots_Fig3.hPlot1),'k');
elseif analysisMeasure == 7    
displayRange(hPlots_Fig3.hPlot1,[16 16],getYLims(hPlots_Fig3.hPlot1),'k');
end

tickLengthPlot = 2*get(hPlots_Fig3.hPlot2(1),'TickLength');

for idx =1:2
    set(hPlots_Fig3.hPlot1(idx),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
end
if analysisMeasure == 1||analysisMeasure == 2
    xlabel(hPlots_Fig3.hPlot1(1),'Time(ms)');
elseif analysisMeasure == 4||analysisMeasure == 5|| analysisMeasure == 6|| analysisMeasure == 7
    xlabel(hPlots_Fig3.hPlot1(1),'Frequency(Hz)');
end

if analysisMeasure == 1
    ylabel(hPlots_Fig3.hPlot1(1),'ERP (\mu V)');ylabel(hPlots_Fig3.hPlot2(1),'RMS value of ERP');
elseif analysisMeasure == 2 
    ylabel(hPlots_Fig3.hPlot1(1),'Spikes/s');ylabel(hPlots_Fig3.hPlot2(1),'spikes/s');
elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6|| analysisMeasure == 7
    if analysisMethod ==1
        ylabel(hPlots_Fig3.hPlot1(1),'log_1_0 (FFT Amplitude)'); ylabel(hPlots_Fig3.hPlot2(1),'log_1_0 (FFT Amplitude)');
        if relativeMeasuresFlag
            ylabel(hPlots_Fig3.hPlot1(1),'log_1_0 (\Delta FFT Amplitude)'); ylabel(hPlots_Fig3.hPlot2(1),'log_1_0 (\Delta FFT Amplitude)')
        end

    elseif analysisMethod ==2
        ylabel(hPlots_Fig3.hPlot1(1),'log_1_0 (Power)');ylabel(hPlots_Fig3.hPlot2(1),'log_1_0 (Power)')
        if relativeMeasuresFlag
            ylabel(hPlots_Fig3.hPlot1(1),'Change in Power (dB)');  ylabel(hPlots_Fig3.hPlot2(1),'Change in Power (dB)');
        end
    end

end
%             set(hOtherMeaures(1),'XScale','linear')
% text(0.5,0.3,['N = ' num2str(dataSize(1))],'color','k','unit','normalized','fontSize',14,'fontWeight','bold','parent',hPlots_Fig3.hPlot2(1))
set(hPlots_Fig3.hPlot2(1),'fontSize',14,'TickDir','out','Ticklength',tickLengthPlot,'box','off')
set(hPlots_Fig3.hPlot2(1),'XTick',cValsUnique,'XTickLabelRotation',90,'XTickLabel',cValsUnique);
xlabel(hPlots_Fig3.hPlot2(1),'Ori 1 Contrast (%)');
title(hPlots_Fig3.hPlot2(1),'CRF along Ori 1')
end

% Accessory Functions

% Load LFP Info
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored); %#ok<NODEF>
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end

% Get parameter combinations
function [parameterCombinations,parameterCombinations2,...
    aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,...
    cValsUnique,tValsUnique,aValsUnique2,eValsUnique2,sValsUnique2,...
    fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = ...
    loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;            end

end

% Get Bad Trials
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end

% Get Color String
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end

% Get FileNamesList
function [fileNameStringAll,fileNameStringListAll,fileNameStringListArray] = getFileNameStringList

[tmpFileNameStringList,monkeyNameList] = getNormalizationExperimentDetails;

fileNameStringAll = ''; pos=1;
clear fileNameStringListArray

for i=1:length(monkeyNameList)
    for j=1:length(tmpFileNameStringList{i})
        fileNameStringAll = [cat(2,fileNameStringAll,tmpFileNameStringList{i}{j}) '|'];
        fileNameStringListAll{pos} = tmpFileNameStringList{i}(j);
        fileNameStringListArray{pos} = tmpFileNameStringList{i}(j); %#ok<*AGROW>
        pos=pos+1;
    end
end

allNames = [];
for i=1:length(monkeyNameList)
    fileNameStringAll = [cat(2,fileNameStringAll,monkeyNameList{i}) ' (N=' num2str(length(tmpFileNameStringList{i})) ')|'];
    fileNameStringListAll{pos} = {[monkeyNameList{i} ' (N=' num2str(length(tmpFileNameStringList{i})) ')']};
    fileNameStringListArray{pos} = tmpFileNameStringList{i};
    allNames = cat(2,allNames,tmpFileNameStringList{i});
    pos=pos+1;
end

fileNameStringAll = cat(2,fileNameStringAll,['all (N=' num2str(length(allNames)) ')']);
fileNameStringListAll{pos} = {['all (N=' num2str(length(allNames)) ')']};
fileNameStringListArray{pos} = allNames;
end

% Get ElectrodesList
function [ElectrodeStringListAll,ElectrodeArrayListAll] = getElectrodesList(fileNameStringTMP,oriSelectiveFlag,folderSourceString)

[tmpElectrodeStringList,tmpElectrodeArrayList,allElecs] = getGoodElectrodesDetails(fileNameStringTMP,oriSelectiveFlag,folderSourceString);

if length(tmpElectrodeStringList)> 1
   clear tmpElectrodeStringList
   tmpElectrodeStringList = {['all (N=' num2str(allElecs) ')']};
end

ElectrodeStringListAll = tmpElectrodeStringList;
ElectrodeArrayListAll = tmpElectrodeArrayList;
end

% Get Induced LFP data by subtracting trialaveraged ERP data from trialwise LFP Data
function Y = removeERP(X)
Y = X-repmat(mean(X,1),size(X,1),1);
end

% Get MeanEnergy for different frequency bands
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange)

posToAverage = intersect(find(freq>=freqRange(1)),find(freq<=freqRange(2)));
eValue   = mean(mEnergy(posToAverage));
end

% Normalize data for ERP and Spike data
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

% Get data for a single electrode in a selected session
function data = getDataSingleElec(data,electrodeNum,analysisMeasure)
    if analysisMeasure == 1 || analysisMeasure ==2
    data.data = data.data(electrodeNum,:,:,:,:);
    data.analysisDataBL = data.analysisDataBL(electrodeNum,:,:,:);
    data.analysisDataST = data.analysisDataST(electrodeNum,:,:,:);
    
    elseif analysisMeasure == 4 || analysisMeasure == 5 || analysisMeasure == 6||analysisMeasure == 7
        data.dataBL = data.dataBL(electrodeNum,:,:,:,:);
        data.dataST = data.dataST(electrodeNum,:,:,:,:);
        for i = 1:length(data.analysisDataST)
            data.analysisDataBL{i} = data.analysisDataBL{i}(electrodeNum,:,:,:);
            data.analysisDataST{i} = data.analysisDataST{i}(electrodeNum,:,:,:);
        end
    end
end

% get Common Baseline across all 5 (Ori 1) x 5 (Ori 2) contrast conditions
function data_BL = getCommonBaseline(data_BL)

if iscell(data_BL)
    size_data_BL = numel(size(data_BL{1}));
    num_con_Ori1 = size(data_BL{1},3);
    num_con_Ori2 = size(data_BL{1},4);
else
    size_data_BL = numel(size(data_BL));
    num_con_Ori1 = size(data_BL,3);
    num_con_Ori2 = size(data_BL,4);
end

if size_data_BL == 4 % baseline for analysis data (elec x TF x Num_Contrast_Ori1 x Num_Contrast_Ori2 x (analysisDataBL))
    if iscell(data_BL)
        for iElec = 1:size(data_BL{1},1)
            for iTF = 1:size(data_BL{1},2)
                for k = 1:length(data_BL)
                    data_BL{k}(iElec,iTF,:,:) = repmat(mean(mean(squeeze(data_BL{k}(iElec,iTF,:,:)),2),1),[num_con_Ori1 num_con_Ori2]);
                end
            end
        end
    else
        for iElec = 1:size(data_BL,1)
            for iTF = 1:size(data_BL,2)
               data_BL(iElec,iTF,:,:) = repmat(mean(mean(squeeze(data_BL(iElec,iTF,:,:)),2),1),[num_con_Ori1 num_con_Ori2]);
            end
        end
    end
elseif size_data_BL == 5 % baseline for timeSeries/PSD data (elec x TF x Num_Contrast_Ori1 x Num_Contrast_Ori2 x time/FreqVals x (dataBL))
    for iElec = 1:size(data_BL,1)
        for iTF = 1:size(data_BL,2)
           data_BL(iElec,iTF,:,:,:) = reshape(repmat(squeeze(mean(mean(squeeze(data_BL(iElec,iTF,:,:,:)),2),1))',[num_con_Ori1*num_con_Ori2 1]),[num_con_Ori1 num_con_Ori2 size(data_BL,5)]);
        end
    end    
end
end

% Segregate data into Preferred (positive x-axis -- increasing contrast) and Null (positive y axis -- increasing contrast)
function data = segregate_Pref_Null_data(data,elecs_neededtoFlipped)
    for iElec = 1:length(elecs_neededtoFlipped)
        for iTF = 1:2
            disp ([num2str(iElec), ' ' num2str(iTF)])
            if numel(fieldnames(data)) == 6
                data.data(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.data(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
            elseif numel(fieldnames(data)) == 8
                data.dataBL(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.dataBL(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
                data.data_cBL(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.data_cBL(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
                data.dataST(elecs_neededtoFlipped(iElec),iTF,:,:,:) = flip(flip(permute(squeeze(data.dataST(elecs_neededtoFlipped(iElec),iTF,:,:,:)),[2 1 3]),1),2);
            end
            
            if ~iscell(data.analysisDataBL) && ~iscell(data.analysisDataST)
                data.analysisDataBL(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataBL(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                data.analysisData_cBL(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisData_cBL(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                data.analysisDataST(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataST(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
            elseif iscell(data.analysisDataBL) && iscell(data.analysisDataST)
                for iCell = 1:length(data.analysisDataBL)
                    data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataBL{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                    data.analysisData_cBL{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisData_cBL{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                    data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:) = flip(flip(permute(squeeze(data.analysisDataST{iCell}(elecs_neededtoFlipped(iElec),iTF,:,:)),[2 1]),1),2);
                end
            end
        end
    end
end

% computing NI for different neural measures
function [NI_absolute,NI_relative] = getNI(data)
if iscell(data.analysisDataST)
    for iElec = 1:size(data.analysisDataST{1},1)
        for iTF = 1:size(data.analysisDataST{1},2)
            for k = 1:length(data.analysisDataST) 
                clear responseMatrix_elec
                responseMatrix_elec = squeeze(data.analysisDataST{k}(iElec,iTF,:,:));
                diff_responseMatrix_elec = squeeze(data.analysisDataST{k}(iElec,iTF,:,:))-squeeze(data.analysisDataBL{k}(iElec,1,:,:));
                NI_absolute{k}(iElec,iTF) = 2*responseMatrix_elec(1,5)/(responseMatrix_elec(1,1)+responseMatrix_elec(5,5))-1;
                NI_relative{k}(iElec,iTF) = 2*diff_responseMatrix_elec(1,5)/(diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))-1;
            end
        end
    end
else
    for iElec = 1:size(data.analysisDataST,1)
        for iTF = 1:size(data.analysisDataST,2)
            clear responseMatrix_elec
            responseMatrix_elec = squeeze(data.analysisDataST(iElec,iTF,:,:));
            diff_responseMatrix_elec = squeeze(data.analysisDataST(iElec,iTF,:,:))-squeeze(data.analysisDataBL(iElec,1,:,:));
            NI_absolute(iElec,iTF) = 2*responseMatrix_elec(1,5)/(responseMatrix_elec(1,1)+responseMatrix_elec(5,5))-1;
            NI_relative(iElec,iTF) = 2*diff_responseMatrix_elec(1,5)/(diff_responseMatrix_elec(1,1)+diff_responseMatrix_elec(5,5))-1;
        end
    end
end
end

% Draw lines for timeTange or FreqRange
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




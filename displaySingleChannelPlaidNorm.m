% Displays Plaid normalization data from a single electrode in Macaque V1

function displaySingleChannelPlaidNorm(subjectName,expDate,protocolName,folderSourceString,gridType)

if ~exist('folderSourceString','var');  folderSourceString='F:';        end
if ~exist('gridType','var');             gridType='ECoG';               end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');

% load LFP Information
[analogChannelsStored,timeVals,~,analogInputNums] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitIDs] = loadspikeInfo(folderSpikes);

% Get Combinations
[~,~,aValsUnique,eValsUnique,sValsUnique,fValsUnique,oValsUnique,cValsUnique,tValsUnique, ...
    aValsUnique2,eValsUnique2,sValsUnique2,fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = loadParameterCombinations(folderExtract);
if aValsUnique ~= aValsUnique2 || eValsUnique ~= eValsUnique2
    error('Azimuths and/or elevations do not match!');
end

% Get properties of the Stimulus
% stimResults = loadStimResults(folderExtract);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
% fonts
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Panels
panelHeight = 0.34; panelStartHeight = 0.61;
staticPanelWidth = 0.25; staticStartPos = 0.025;
dynamicPanelWidth = 0.25; dynamicStartPos = 0.275;
timingPanelWidth = 0.25; timingStartPos = 0.525;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.775;
backgroundColor = 'w';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.06; dynamicGap=0.015; dynamicTextWidth = 0.5;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Analog channel
[analogChannelStringList,analogChannelStringArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Analog Channel','FontSize',fontSizeSmall);
hAnalogChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analogChannelStringList,'FontSize',fontSizeSmall);

% Neural channel
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Neural Channel','FontSize',fontSizeSmall);
    
if ~isempty(neuralChannelsStored)
    neuralChannelString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs);
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','popup','String',neuralChannelString,'FontSize',fontSizeSmall);
else
    hNeuralChannel = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
        'Position', [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
        'Style','text','String','Not found','FontSize',fontSizeSmall);
end
% Sigma
sigmaString  = getStringFromValues(sValsUnique,1);
sigmaString2 = getStringFromValues(sValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Sigma (Deg)','FontSize',fontSizeSmall);
hSigma = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',sigmaString,'FontSize',fontSizeSmall);
hSigma2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-3*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',sigmaString2,'FontSize',fontSizeSmall);

% Spatial Frequency
spatialFreqString  = getStringFromValues(fValsUnique,1);
spatialFreqString2 = getStringFromValues(fValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Spatial Freq (CPD)','FontSize',fontSizeSmall);
hSpatialFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',spatialFreqString,'FontSize',fontSizeSmall);
hSpatialFreq2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-4*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',spatialFreqString2,'FontSize',fontSizeSmall);

% Orientation
orientationString  = getStringFromValues(oValsUnique,1);
orientationString2 = getStringFromValues(oValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Orientation (Deg)','FontSize',fontSizeSmall);
hOrientation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',orientationString,'FontSize',fontSizeSmall);
hOrientation2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-5*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',orientationString2,'FontSize',fontSizeSmall);

% Contrast
contrastString  = 'all'; %getStringFromValues(cValsUnique,1);
contrastString2 = 'all'; %getStringFromValues(cValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);
hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeSmall);
hContrast2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-6*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',contrastString2,'FontSize',fontSizeSmall);

% Temporal Frequency
temporalFreqString  = getStringFromValues(tValsUnique,1);
temporalFreqString2 = getStringFromValues(tValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-7*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeSmall);
hTemporalFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-7*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeSmall);
hTemporalFreq2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-7*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',temporalFreqString2,'FontSize',fontSizeSmall);

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|Raster|FFT & Alpha Power|FFT & Gamma Power|delta FFT|STA|SSVEP';
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-8*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-8*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);

% For all Plots
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-9.5*(dynamicHeight+dynamicGap) 1 dynamicHeight],...
    'Style','text','String','For all plots','FontSize',fontSizeSmall);
% Azimuth
azimuthString  = getStringFromValues(aValsUnique,1);
azimuthString2 = getStringFromValues(aValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-10.5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Azimuth (Deg)','FontSize',fontSizeSmall);
hAzimuth = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-10.5*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',azimuthString,'FontSize',fontSizeSmall);
hAzimuth2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-10.5*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',azimuthString2,'FontSize',fontSizeSmall);

% Elevation
elevationString  = getStringFromValues(eValsUnique,1);
elevationString2 = getStringFromValues(eValsUnique2,1);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-11.5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Elevation (Deg)','FontSize',fontSizeSmall);
hElevation = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-11.5*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',elevationString,'FontSize',fontSizeSmall);
hElevation2 = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth+(1-dynamicTextWidth)/2 1-11.5*(dynamicHeight+dynamicGap) (1-dynamicTextWidth)/2 dynamicHeight], ...
    'Style','popup','String',elevationString2,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [-0.1 0.5];
ERPRange = [0.05 0.2];
fftRange = [0 250];
baseline = [-0.2 0];
stimPeriod = [0.2 0.4];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% ERP Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-4*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','ERP Range (s)','FontSize',fontSizeSmall);
hERPMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(ERPRange(1)),'FontSize',fontSizeSmall);
hERPMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-4*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(ERPRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

% STA length
staLen = [-0.05 0.05]; 
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','STA len (s)','FontSize',fontSizeSmall);
hSTAMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(1)),'FontSize',fontSizeSmall);
hSTAMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(staLen(2)),'FontSize',fontSizeSmall);
hRemoveMeanSTA = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0 1-10*timingHeight 1 timingHeight], ...
    'Style','togglebutton','String','remove mean STA','FontSize',fontSizeMedium);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeSmall);

hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get plots and message handles

% Get electrode array information
electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,analogChannelsStored(get(hAnalogChannel,'val')), ...
    colorNames(get(hChooseColor,'val')),[],1,0,gridType,subjectName);

mapRatio = 1/2;
startXPos = staticStartPos; endXPos = 0.95; startYPos = 0.05; mainRFHeight = 0.48; centerGap = 0.05;
mainRFWidth = mapRatio*(endXPos-startXPos-centerGap);

% Main plot handles 
numRows = length(cValsUnique); numCols = length(cValsUnique2);
gridPos=[0.03+startXPos startYPos mainRFWidth mainRFHeight]; gap = 0.002;
plotHandles = getPlotHandles(numRows,numCols,gridPos,gap);

gridPos2 = [endXPos-mainRFWidth+0.02 startYPos+mainRFHeight*3/5 mainRFWidth 2*mainRFHeight/5]; gap = 0.002;
plotHandles2 = getPlotHandles(2,numCols,gridPos2,gap);

gridPos3 = [endXPos-mainRFWidth+0.02 startYPos mainRFWidth 2*mainRFHeight/5]; gap = 0.002;
plotHandles3 = getPlotHandles(1,2,gridPos3,gap);

uicontrol('Unit','Normalized','Position',[0 0.975 1 0.03],...
    'Style','text','String',[subjectName expDate protocolName],'FontSize',fontSizeLarge);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        a=[get(hAzimuth,'val') get(hAzimuth2,'val')];
        e=[get(hElevation,'val') get(hElevation2,'val')];
        s=[get(hSigma,'val') get(hSigma2,'val')];
        f=[get(hSpatialFreq,'val') get(hSpatialFreq2,'val')];
        o=[get(hOrientation,'val') get(hOrientation2,'val')];
        c=[get(hContrast,'val') get(hContrast2,'val')];
        t=[get(hTemporalFreq,'val') get(hTemporalFreq2,'val')];
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];
        ERPRange = [str2double(get(hERPMin,'String')) str2double(get(hERPMax,'String'))];
        staRange = [str2double(get(hSTAMin,'String')) str2double(get(hSTAMax,'String'))];
        holdOnState = get(hHoldOn,'val');
        removeMeanSTA = get(hRemoveMeanSTA,'val');

        if analysisType==7 % Spike triggered average
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};
            spikeChannelPos = get(hNeuralChannel,'val');
            spikeChannelNumber = neuralChannelsStored(spikeChannelPos);
            unitID = SourceUnitIDs(spikeChannelPos);
            
            plotColors{1} = 'g';
            plotColors{2} = 'k';
            plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                a,e,s,f,o,c,t,timeVals,plotColors,blRange,stRange,folderName,staRange,removeMeanSTA);
            
          
            if analogChannelPos<=length(analogChannelsStored)
                analogChannelNumber = analogChannelsStored(analogChannelPos);
            else
                analogChannelNumber = 0;
            end
            channelNumber = [analogChannelNumber spikeChannelNumber];
            
        elseif analysisType == 2 || analysisType == 3
            channelPos = get(hNeuralChannel,'val');
            channelNumber = neuralChannelsStored(channelPos);
            unitID = SourceUnitIDs(channelPos);
            plotSpikeData1Channel(plotHandles,plotHandles2,plotHandles3,channelNumber,a,e,s,f,o,c,t,folderSpikes,...
                analysisType,timeVals,plotColor,unitID,blRange,stRange,folderName);

           
        else
            analogChannelPos = get(hAnalogChannel,'val');
            analogChannelString = analogChannelStringArray{analogChannelPos};

            plotLFPData1Channel(plotHandles,plotHandles2,plotHandles3,analogChannelString,a,e,s,f,o,c,t,folderLFP,...
                analysisType,timeVals,plotColor,blRange,stRange,ERPRange,folderName);


            if analogChannelPos<=length(analogChannelsStored)
                channelNumber = analogChannelsStored(analogChannelPos);
            else
                channelNumber = 0;
            end
            

        end

        if analysisType<=3  % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        elseif analysisType <=6 || analysisType == 8 % LFP alpha power & gamma power or SSVEP for CP flicker stimuli
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        else
            xMin = str2double(get(hSTAMin,'String'));
            xMax = str2double(get(hSTAMax,'String'));
        end

        rescaleData(plotHandles,xMin,xMax,getYLims(plotHandles));
        rescaleData(plotHandles2,0,100,getYLims(plotHandles2));
        rescaleData(plotHandles3,0,100,getYLims(plotHandles3));

        showElectrodeLocations(electrodeGridPos,channelNumber,plotColor,hElectrodes,holdOnState,0,gridType,subjectName);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)

        analysisType = get(hAnalysisType,'val');
        
        if analysisType<=3 % ERP or spikes
            xMin = str2double(get(hStimMin,'String'));
            xMax = str2double(get(hStimMax,'String'));
        else
            xMin = str2double(get(hFFTMin,'String'));
            xMax = str2double(get(hFFTMax,'String'));
        end

        yLims = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,xMin,xMax,yLims);

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function plotLFPData1Channel(plotHandles,plotHandles2,plotHandles3,channelString,a,e,s,f,o,~,t,folderLFP,...
analysisType,timeVals,plotColor,blRange,stRange,ERPRange,folderName)

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

titleFontSize = 12;

[parameterCombinations,parameterCombinations2,~,~,~,~,oValsUnique,cValsUnique,tValsUnique, ...
                                              ~,~,~,~,oValsUnique2,cValsUnique2,tValsUnique2] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

c1List = 1:length(cValsUnique);
c2List = 1:length(cValsUnique2);

% Get the data
removeAvgRef = 0;
if removeAvgRef
    disp('Removing average reference');
    load(fullfile(folderLFP,'avgRef.mat'));
    avgRef = analogData;
end
clear signal analogData
load(fullfile(folderLFP,channelString));
if removeAvgRef
    analogData = analogData-avgRef;
end

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end



rfMapVals = zeros(numRows,numCols);
for c1=1:numRows
    for c2=1:numCols

        clear goodPos
        goodPos = parameterCombinations{a,e,s(1),f(1),o(1),c1List(c1),t(1)};
        goodPos2 = parameterCombinations2{a,e,s(2),f(2),o(2),c2List(c2),t(2)};
        goodPos = intersect(goodPos,goodPos2);
        goodPos = setdiff(goodPos,badTrials);
      
        if isempty(goodPos)
            disp('No entries for this combination..');
        else
            disp(['pos=(' num2str(c1) ',' num2str(c2) ') ,n=' num2str(length(goodPos))]);
    
            Fs = round(1/(timeVals(2)-timeVals(1)));
            if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
                disp('baseline and stimulus ranges are not the same');
            else
                range = blRange;
                rangePos = round(diff(range)*Fs);
                ERPRangePos = round(diff(ERPRange)*Fs);
                blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
                stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
                ERPPos = find(timeVals>=ERPRange(1),1)+ (1:ERPRangePos);
                xs = 0:1/diff(range):Fs-1/diff(range);
            end

            if analysisType == 1        % compute ERP
                clear erp
                erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
                plot(plotHandles(c1,c2),timeVals,erp,'color',plotColor);
                ERPVals(c1,c2) = rms(erp(ERPPos));
                rfMapVals(e,a) = rms(erp(stPos));

            elseif analysisType == 2  ||   analysisType == 3 % compute Firing rates
                disp('Use plotSpikeData instead of plotLFPData...');
            else
                fftBL = abs(fft(analogData(goodPos,blPos),[],2));
                fftST = abs(fft(analogData(goodPos,stPos),[],2));
                alphaRange = [8 12]; gammaRange = [30 60];
                alphaPos = find(xs >= alphaRange(1) & xs <= alphaRange(2));
                gammaPos = find(xs >= gammaRange(1) & xs <= gammaRange(2));
                                
                if analysisType == 4 || analysisType == 5 || analysisType == 8 % compute fft for power spectrum & power estimation
                    plot(plotHandles(c1,c2),xs,log10(mean(fftBL,1)),'g');
                    set(plotHandles(c1,c2),'Nextplot','add');
                    plot(plotHandles(c1,c2),xs,log10(mean(fftST,1)),'k');
                    set(plotHandles(c1,c2),'Nextplot','replace');
                    
                    if analysisType == 4 
                        %alpha Power Estimation
                        if length(alphaPos) > 1
                            alphaPowerBL(c1,c2) = log10(mean(mean(fftBL(:,alphaPos),2),1));
                            alphaPowerST(c1,c2) = log10(mean(mean(fftST(:,alphaPos),2),1));
                        elseif length(alphaPos) == 1
                            alphaPowerBL(c1,c2) = log10(mean(fftBL(:,alphaPos)));
                            alphaPowerST(c1,c2) = log10(mean(fftST(:,alphaPos)));
                        elseif alphaPos == 0
                            error('No alphaPos found for selected BL or ST Range!')
                        end
                        
                    elseif analysisType == 5
                        
                        if length(gammaPos) > 1
                            gammaPowerBL(c1,c2) = log10(mean(mean(fftBL(:,gammaPos),2),1));
                            gammaPowerST(c1,c2) = log10(mean(mean(fftST(:,gammaPos),2),1));
                        elseif length(gammaPos) == 1
                            gammaPowerBL(c1,c2) = log10(mean(fftBL(:,gammaPos)));
                            gammaPowerST(c1,c2) = log10(mean(fftST(:,gammaPos)));
                        elseif gammaPos == 0
                            error('No alphaPos found for selected BL or ST Range!')
                        end
                        
                    elseif analysisType == 8
                            
                        if t(1)>1 || t(2)>1
                            t1 = tValsUnique(t(1));
                            t2 = tValsUnique2(t(2));
                            tempFreqCP = max(t1,t2); %counterphase Flickering stimuli
                            ssvepFreq = 2*tempFreqCP;   % 2nd Harmonic
                            ssvepPos = find(xs == ssvepFreq);
                            if isempty(ssvepPos)
                                error('Adjust the BL and ST Range to incorporate SSVEP Freq')
                            end
                            ssvepPowerBL(c1,c2) = log10(mean(fftBL(:,ssvepPos)));
                            ssvepPowerST(c1,c2) = log10(mean(fftST(:,ssvepPos)));
                        end
                    end
                end
            end
        end

                if analysisType == 6
                    plot(plotHandles(c1,c2),xs,log10(mean(fftST,1))-log10(mean(fftBL,1)),'color',plotColor); % compute deltaFFT
                end
    end
end
 
CRFColors = jet(length(cValsUnique));
if analysisType==1
%     rfMapVals=[];
    for iCon = 1:length(cValsUnique2)
        plot(plotHandles2(1,iCon),cValsUnique,ERPVals(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles2(2,iCon),cValsUnique2,ERPVals(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles3(1,1),cValsUnique,ERPVals(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,1),'on')
        plot(plotHandles3(1,2),cValsUnique2,ERPVals(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,2),'on')
    end
    hold(plotHandles3(1,1),'off');hold(plotHandles3(1,2),'off')
end

if analysisType==4

    for iCon = 1:length(cValsUnique2)
        plot(plotHandles2(1,iCon),cValsUnique,alphaPowerST(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles2(2,iCon),cValsUnique2,alphaPowerST(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles3(1,1),cValsUnique,alphaPowerST(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,1),'on')
        plot(plotHandles3(1,2),cValsUnique2,alphaPowerST(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,2),'on')
    end
    hold(plotHandles3(1,1),'off');hold(plotHandles3(1,2),'off')
end

if analysisType==5

    for iCon = 1:length(cValsUnique2)
        plot(plotHandles2(1,iCon),cValsUnique,gammaPowerST(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles2(2,iCon),cValsUnique2,gammaPowerST(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles3(1,1),cValsUnique,gammaPowerST(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,1),'on')
        plot(plotHandles3(1,2),cValsUnique2,gammaPowerST(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,2),'on')
    end
    hold(plotHandles3(1,1),'off');hold(plotHandles3(1,2),'off')
end

if analysisType==8

    for iCon = 1:length(cValsUnique2)
        plot(plotHandles2(1,iCon),cValsUnique,ssvepPowerST(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles2(2,iCon),cValsUnique2,ssvepPowerST(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles3(1,1),cValsUnique,ssvepPowerST(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,1),'on')
        plot(plotHandles3(1,2),cValsUnique2,ssvepPowerST(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,2),'on')
    end
    hold(plotHandles3(1,1),'off');hold(plotHandles3(1,2),'off')
end

% Display Axis and titles
for c1 = 1:5
    for c2 = 1:5
        if c1 == 1
            title(plotHandles(c1,c2),[num2str(cValsUnique2(c2)) '%'],'FontSize',titleFontSize)
        end
        if c2 == 1
            ylabel(plotHandles(c1,c2),[num2str(cValsUnique(c1)) '%'],'FontSize',titleFontSize,'FontWeight','bold')
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSpikeData1Channel(plotHandles,plotHandles2,plotHandles3,channelNumber,a,e,s,f,o,~,t,folderSpikes,...
analysisType,timeVals,plotColor,unitID,blRange,stRange,folderName)

titleFontSize = 12;

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

[parameterCombinations,parameterCombinations2,~,~,~,~,oValsUnique,cValsUnique,~, ...
                                              ~,~,~,~,oValsUnique2,cValsUnique2,~] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

c1List = 1:length(cValsUnique);
c2List = 1:length(cValsUnique2);

% Get the data
clear signal spikeData
load(fullfile(folderSpikes,['elec' num2str(channelNumber) '_SID' num2str(unitID) '.mat']));

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

for c1=1:numRows
    for c2=1:numCols

        clear goodPos
        goodPos = parameterCombinations{a,e,s(1),f(1),o(1),c1List(c1),t(1)};
        goodPos2 = parameterCombinations2{a,e,s(2),f(2),o(2),c2List(c2),t(2)};
        goodPos = intersect(goodPos,goodPos2);
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            disp(['pos=(' num2str(c1) ',' num2str(c2) ') ,n=' num2str(length(goodPos))]);
            
            if analysisType == 2
                [psthVals,xs] = getPSTH(spikeData(goodPos),10,[timeVals(1) timeVals(end)]);
                plot(plotHandles(c1,c2),xs,psthVals,'color',plotColor);
                
                FRVals(c1,c2) = mean(getSpikeCounts(spikeData(goodPos),stRange))/diff(stRange);
            else
                X = spikeData(goodPos);
                axes(plotHandles(c1,c2)); %#ok<LAXES>
                rasterplot(X,1:length(X),plotColor);
            end
        end
    end
end

CRFColors = jet(length(cValsUnique));
if analysisType==2
%     rfMapVals=[];
    for iCon = 1:length(cValsUnique2)
        plot(plotHandles2(1,iCon),cValsUnique,FRVals(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles2(2,iCon),cValsUnique2,FRVals(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        plot(plotHandles3(1,1),cValsUnique,FRVals(iCon,:),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,1),'on')
        plot(plotHandles3(1,2),cValsUnique2,FRVals(:,iCon),'Marker','o','LineWidth',2,'color',CRFColors(iCon,:,:))
        hold(plotHandles3(1,2),'on')
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSTA1Channel(plotHandles,analogChannelString,spikeChannelNumber,unitID,folderLFP,folderSpikes,...
                s,f,o,c,t,timeVals,plotColors,blRange,stRange,folderName,staLen,removeMeanSTA)

titleFontSize = 12;

folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');

[parameterCombinations,parameterCombinations2,aValsUnique,eValsUnique] = loadParameterCombinations(folderExtract);
[numRows,numCols] = size(plotHandles);

% Get the analog data
clear signal analogData
load(fullfile(folderLFP,analogChannelString));

% Get the spike data
clear signal spikeData
load(fullfile(folderSpikes,['elec' num2str(spikeChannelNumber) '_SID' num2str(unitID) '.mat']));

% Get bad trials
badTrialFile = fullfile(folderSegment,'badTrials.mat');
if ~exist(badTrialFile,'file')
    disp('Bad trial file does not exist...');
    badTrials=[];
else
    badTrials = loadBadTrials(badTrialFile);
    disp([num2str(length(badTrials)) ' bad trials']);
end

staTimeLims{1} = blRange;
staTimeLims{2} = stRange;

for i=1:numRows
    e = numRows-i+1;
    for j=1:numCols
        a = j;
        clear goodPos
        goodPos = parameterCombinations{a,e,s(1),f(1),o(1),c(1),t(1)};
        goodPos2 = parameterCombinations2{a,e,s(2),f(2),o(2),c(2),t(2)};
        goodPos = intersect(goodPos,goodPos2);
        goodPos = setdiff(goodPos,badTrials);

        if isempty(goodPos)
            disp('No entries for this combination..')
        else
            goodSpikeData = spikeData(goodPos);
            goodAnalogSignal = analogData(goodPos,:);
            [staVals,numberOfSpikes,xsSTA] = getSTA(goodSpikeData,goodAnalogSignal,staTimeLims,timeVals,staLen,removeMeanSTA);
            
            disp([num2str(i) ' ' num2str(j) ', numStim: ' num2str(length(goodPos)) ', numSpikes: ' num2str(numberOfSpikes)]);
            if ~isempty(staVals{1})
                plot(plotHandles(i,j),xsSTA,staVals{1},'color',plotColors{1});
            end
            set(plotHandles(i,j),'Nextplot','add');
            if ~isempty(staVals{2})
                plot(plotHandles(i,j),xsSTA,staVals{2},'color',plotColors{2});
            end
            set(plotHandles(i,j),'Nextplot','replace');
        end
        
        % Display title
        if (i==1)
            if (j==1)
                title(plotHandles(i,j),['Azi: ' num2str(aValsUnique(a))],'FontSize',titleFontSize);
            else
                title(plotHandles(i,j),num2str(aValsUnique(a)),'FontSize',titleFontSize);
            end
        end

        if (j==numCols)
            if (i==1)
                title(plotHandles(i,j),[{'Ele'} {num2str(eValsUnique(e))}],'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            else
                title(plotHandles(i,j),num2str(eValsUnique(e)),'FontSize',titleFontSize,...
                    'Units','Normalized','Position',[1.25 0.5]);
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outString = getStringFromValues(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end
function [outString,outArray] = getAnalogStringFromValues(analogChannelsStored,analogInputNums)
outString='';
count=1;
for i=1:length(analogChannelsStored)
    outArray{count} = ['elec' num2str(analogChannelsStored(i))]; %#ok<AGROW>
    outString = cat(2,outString,[outArray{count} '|']);
    count=count+1;
end
if ~isempty(analogInputNums)
    for i=1:length(analogInputNums)
        outArray{count} = ['ainp' num2str(analogInputNums(i))]; %#ok<AGROW>
        outString = cat(2,outString,[outArray{count} '|']);
        count=count+1;
    end
end
end
function outString = getNeuralStringFromValues(neuralChannelsStored,SourceUnitIDs)
outString='';
for i=1:length(neuralChannelsStored)
    outString = cat(2,outString,[num2str(neuralChannelsStored(i)) ', SID ' num2str(SourceUnitIDs(i)) '|']);
end 
end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%
%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,timeVals,goodStimPos,analogInputNums] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
analogChannelsStored=sort(analogChannelsStored);
if ~exist('analogInputNums','var')
    analogInputNums=[];
end
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = fullfile(folderSpikes,'spikeInfo.mat');
if exist(fileName,'file')
    load(fileName);
    [neuralChannelsStored,I]=sort(neuralChannelsStored);
    SourceUnitID=SourceUnitID(I);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,parameterCombinations2,aValsUnique,eValsUnique,sValsUnique,...
    fValsUnique,oValsUnique,cValsUnique,tValsUnique,aValsUnique2,eValsUnique2,...
    sValsUnique2,fValsUnique2,oValsUnique2,cValsUnique2,tValsUnique2] = loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));

if ~exist('sValsUnique','var');    sValsUnique=rValsUnique;             end
if ~exist('cValsUnique','var');    cValsUnique=100;                     end
if ~exist('tValsUnique','var');    tValsUnique=0;                       end
end
function badTrials = loadBadTrials(badTrialFile)
load(badTrialFile);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prefOrientation,orientationSelectivity] = getOrientationTuning(computationVals,oValsUnique)
num=0;
den=0;

for j=1:length(oValsUnique)
    num = num+computationVals(j)*sind(2*oValsUnique(j));
    den = den+computationVals(j)*cosd(2*oValsUnique(j));
end

prefOrientation = 90*atan2(num,den)/pi;
orientationSelectivity = abs(den+1i*num)/sum(computationVals);

if prefOrientation<0
    prefOrientation = prefOrientation+180;
end
end
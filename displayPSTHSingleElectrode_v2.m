function displayPSTHSingleElectrode_v2(elecNum)
if strcmp(getenv('username'),'RayLabPC-Aritra') || strcmp(getenv('username'),'Lab Computer-Aritra')
    folderSourceString_Project = 'E:\';
elseif strcmp(getenv('username'),'Supratim Ray')
    folderSourceString_Project = 'M:\';
end
close all;
folderName = fullfile(folderSourceString_Project,'Projects\Aritra_PlaidNormalizationProject\savedData_Figures\');
fileName = fullfile(folderName,'all_N15_S2_allElecs_T150_400_d0_0.75_tapers1_removeERP0_cne0_gse1_Microelectrode_UnitID0.mat');

if exist(fileName,'file')
    load(fullfile(fileName)); %#ok<LOAD>
    disp(['Loading file ' fileName]);
end

% if strcmp(exampleElectrodeType,'Additive')
% elecList = [2 8 10 13 19 23 30 48 57 58 66 77 88 100 117 119 129 158 159 161 162 170 174 177 178 179 182 184]; 
% elseif strcmp(exampleElectrodeType,'Norm')
% elecList = [2 8 10 13 19 23 30 48 57 58 66 77 88 100 117 119 129 158 159 161 162 170 174 177 178 179 182 184]; 
% end

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

figure(4);


for iElec = elecNum: elecNum% 1:length(elecList) %1:size(firingRateData.data,1)
    hPlot = getPlotHandles(numRows,numCols,gridPos,gap); linkaxes(hPlot);

    for cOri2 = 1:5
        for cOri1 = 1:5
            plot(hPlot(cOri2,cOri1),firingRateData.timeVals,squeeze(firingRateData.data(iElec,1,cOri2,cOri1,:)))
            hold(hPlot(cOri2,cOri1),'on');

%             if cOri2 ==1 && cOri1 == 5
%                 plot(hPlot(cOri2,cOri1),firingRateData.timeVals,(squeeze(firingRateData.data(iElec,1,1,1,:))+squeeze(firingRateData.data(iElec,1,5,5,:))),'color',[0.8500, 0.3250, 0.0980],'LineStyle',':','LineWidth',2)
%             end

            ylim(getYLims(hPlot)+[0 20]);
            xlim([-0.3 0.5]);
%             rescaleData(hPlot,-0.1,0.5,getYLims(hPlot),12);
        end
    end
    title(hPlot(1,1),['elec: ' num2str(iElec)]);
    pause;
    hold off;
    clf;
    
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

%Rescale data
function rescaleData(plotHandles,xMin,xMax,yLims,labelSize) %#ok<DEFNU>

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
            set(plotHandles(i,j),'XTickLabel',[],'fontSize',labelSize);
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
%set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
%set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

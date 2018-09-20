function varargout = DynamicPoster_DC_SFN2017(varargin)
% DYNAMICPOSTER_DC_SFN2017 MATLAB code for DynamicPoster_DC_SFN2017.fig
%      DYNAMICPOSTER_DC_SFN2017, by itself, creates a new DYNAMICPOSTER_DC_SFN2017 or raises the existing
%      singleton*.
%
%      H = DYNAMICPOSTER_DC_SFN2017 returns the handle to a new DYNAMICPOSTER_DC_SFN2017 or the handle to
%      the existing singleton*.
%
%      DYNAMICPOSTER_DC_SFN2017('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DYNAMICPOSTER_DC_SFN2017.M with the given input arguments.
%
%      DYNAMICPOSTER_DC_SFN2017('Property','Value',...) creates a new DYNAMICPOSTER_DC_SFN2017 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DynamicPoster_DC_SFN2017_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DynamicPoster_DC_SFN2017_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DynamicPoster_DC_SFN2017

% Last Modified by GUIDE v2.5 07-Nov-2017 11:23:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DynamicPoster_DC_SFN2017_OpeningFcn, ...
                   'gui_OutputFcn',  @DynamicPoster_DC_SFN2017_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before DynamicPoster_DC_SFN2017 is made visible.
function DynamicPoster_DC_SFN2017_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DynamicPoster_DC_SFN2017 (see VARARGIN)

% Choose default command line output for DynamicPoster_DC_SFN2017
handles.output = hObject;

handles.dataFolderArray ={'L01','L02','L03','L04'};
handles.sliceThickness = 20;
handles.update_fts_legend = true;
% handles.heatColor = jet(10);

homepath = getenv('HOME');
fs = filesep;
handles.posterHome = [homepath,fs,'Dropbox',fs,'Talks',fs,'SFN2017_DynamicPoster'];
addpath(handles.posterHome);

% handles.plotFrequencyBins = linspace(0, 0.005,9); %defined in load
% dataset
% handles.plotAmplitudeBins = linspace(0, 0.08, 9);
% handles.plotCorrelationBins = linspace(0, 1, 9);
% handles.plotAgeBins = 0.5;

handles.stackAcqFreq = 1.24;

handles = loadDataset(handles,eventdata);

% handles.recruitMapData = load([handles.posterHome,fs,'Data',fs,'recruitmentMap.mat']);
% load([handles.posterHome,fs,'Data',fs,'grm_eae_dF0.mat']);
% handles.dF0 = dF0;
% handles.ftsRw = alxFTS; %could not save due to large cell size

% handles.recruitmentMap.Amplitude = handles.recruitMapData.rcAmp;
% handles.recruitmentMap.Frequency = handles.recruitMapData.rcFreq;

cornellLogo = imread([handles.posterHome,fs,'CornellLOGO.png'] );
bPx = find(cornellLogo == 0);
cornellLogo(bPx) = 255;
axes(handles.cornellLogo)
imshow(cornellLogo,'Parent',handles.cornellLogo);
axis off

% zfishSchematic = imread([handles.posterHome,fs,'zfish_schematic.png'] );
% axes(handles.zfishSchematic);
% image(zfishSchematic,'Parent', handles.zfishSchematic);
% axis off; axis square

axes(handles.static_txRed_AlxBacfill_Image);
alxBackfill = imread([handles.posterHome,fs,'AlxKaede_txRedBacfill_xuvStitch.png'] );
imshow(alxBackfill);
% xlabel('TxRed RS Backfill','color','r','FontUnits','normalized','FontSize',0.06,'Parent',handles.static_txRed_AlxBacfill_Image);
% ylabel('Chx10 Kaede','color','g','FontUnits','normalized','FontSize',0.06,'Parent',handles.static_txRed_AlxBacfill_Image);


alxColorChange = imread([handles.posterHome,fs,'AlxKaede_ColorChange.png'] );
imshow(alxColorChange,'InitialMagnification',100,'Parent',handles.colorChange);

%%%%heat map panel
stckRateHz = 1.24; 

load([handles.posterHome,fs,'triggPanelData.mat'] );
axes(handles.triggAvgHeatMap); cla;
imagesc(trgShck);
caxis([0 1]);
set(gca,'XTick',linspace(1,24,5),'XTickLabel',...
    round( [linspace(-12,12,5)/stckRateHz]*10 )/10,...
    'XColor',[1,1,1],'YColor',[1,1,1],'FontUnits','normalized','FontSize',0.07);
hcb=colorbar;
hcb.Color = [1,1,1];
% hcb.Box = 'off';
ylabel(hcb, '\DeltaF/F');

ylabel('Cell #');
xlabel('Time to Shock (s)');
title('Shock Triggered Average','color',[1,1,1]);
% box off

axes(handles.shockHeat);
imagesc(dffz);
hold on, trigLocation = trigLocations(1).ShockTrigger; %from L04 Shock01
trigLocation_S = round(trigLocation/stckRateHz );
% for k=1:length(trigLocation)
%     plot([trigLocation(k),trigLocation(k)],get(gca,'YLim'),'w-','LineWidth',0.05); 
% end
set(gca,'XTick',trigLocation,'XTickLabel',trigLocation_S,'TickDir','out','TickLength',[0.01, 0.025],'XColor',[1,1,1],'YColor',[1,1,1]);
caxis([0 1]);
set(gca,'FontUnits','normalized','FontSize',0.07);
xtickangle(45)
% colormap hot
xlabel('Times of Shock Stimuli (s)');
ylabel('Cell #');
title('\Delta F/F','color',[1,1,1]);
box off

load ventral_root_quant_data.mat
createfigure_ventralroot(vxr,vxf,handles.ventral_root_processing);
%%%%%%
handles.ls_videoReader = VideoReader([handles.posterHome,fs,'Data',fs,'cAlx_Spots_Segmentation.avi']);
sdat = readFrame(handles.ls_videoReader); %using read will load the WHOLE
% movie into memory. 
image(sdat,'Parent',handles.shockSpinner);
handles.shockSpinner.Visible = 'off';
%%%%%%%%%%%%%%%%%

axes(handles.slicePosMap);
% Update handles structure
guidata(hObject, handles);


% This sets up the initial plot - only do when we are invisible
% so window can get raised using DynamicPoster_DC_SFN2017.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

function handles=loadDataset(handles,eventdata)
%switches data set
fs = filesep;
dataFolder = handles.dataFolderArray{get(handles.Dataset_Selection,'Value')};

handles.recruitMapData = load([handles.posterHome,fs,'Data',fs,dataFolder,fs,'recruitmentMap.mat']);
load([handles.posterHome,fs,'Data',fs,dataFolder,fs,'grm_eae_dF0.mat']);
handles.dF0 = dF0;

handles.stackX = [1:size(handles.dF0{1},2)]/handles.stackAcqFreq;

if ~exist([handles.posterHome,fs,'Data',fs,dataFolder,fs,'costesRecruitmentThreshold.mat'],'file');
    
    amplitudeFits = [handles.recruitMapData.cellFit.Amplitude];
    integAmplitudeSlopes = amplitudeFits(:,2:2:end);
    corrValue = handles.recruitMapData.cellCorrelation.maxCorr;
    
    for k=1:numel(handles.dF0)
        [T1,T2,cellReject] = costesThreshold(integAmplitudeSlopes,corrValue,k,0.05);
        costeThreshold(k,:) = [T1,T2]; %amplitude, frequency
        cellRecruitment(:,k) = ~ cellReject;
    end
    save([handles.posterHome,fs,'Data',fs,dataFolder,fs,'costesRecruitmentThreshold.mat'],'costeThreshold','cellRecruitment');
else
    load([handles.posterHome,fs,'Data',fs,dataFolder,fs,'costesRecruitmentThreshold.mat']);
end

handles.costeThreshold = costeThreshold;
cellRecruitment(1:end) = 1; %let's just display all cells for now. Don't know how to deal...
handles.cellRecruitment = cellRecruitment;

handles.spPos = handles.recruitMapData.posReference.pos;
handles.alx = handles.recruitMapData.posReference.alx;
handles.alxPos = handles.spPos(handles.alx,:);
handles.alxAge = handles.recruitMapData.posReference.alxAge(handles.alx);
handles.inSlice = zeros(size(handles.alx,1),1,'logical');


%for faster loading, generate a compressed data set that has only above
%threshold dFFs
%also, add dots representing the dFF/motor metric combo in the raw traces
%so that people can see where the data comes from. 

handles.trial2use = 1;

% axes(handles.slicePosMap); cla;
% hold on
% plot(handles.spPos(:,1),handles.spPos(:,2),'k.');
% plot(handles.alxPos(:,1),handles.alxPos(:,2),'y.');
% axis equal
% axis tight

% 
% handles.recruitFunction = 'Amplitude';
% heatBins = handles.plotAmplitudeBins;
% handles.heatBins = [-inf heatBins inf];
% handles.heatColor = jet(10);

amplitudeFits = [handles.recruitMapData.cellFit.Amplitude];
frequencyFits = [handles.recruitMapData.cellFit.Frequency];

handles.recruitmentMap.Amplitude = amplitudeFits(:,2:2:end);
handles.recruitmentMap.Frequency = frequencyFits(:,2:2:end);
handles.recruitmentMap.Correlation = handles.recruitMapData.cellCorrelation.maxCorr;
handles.recruitmentMap.xCorrLag = handles.recruitMapData.cellCorrelation.maxCorr_lag;

handles.plotFrequencyBins = linspace(0, 0.5*max(handles.recruitmentMap.Frequency(1:end)),9);
handles.plotAmplitudeBins = linspace(0, 0.5*max(handles.recruitmentMap.Amplitude(1:end)), 9);
handles.plotCorrelationBins = linspace(min(handles.recruitmentMap.Correlation(1:end)), max(handles.recruitmentMap.Correlation(1:end)), 9);
handles.plotAgeBins = 0.5;

handles.linearRange = floor( min(handles.alxPos(:,1)) ):handles.sliceThickness: ceil( max(handles.alxPos(:,1)) );
handles.linearRange = [handles.linearRange inf];
handles.currentSlice = ceil( length(handles.linearRange)/2 );

set(handles.sliceSelector,'Max',length(handles.linearRange)-1);
set(handles.sliceSelector,'Value',handles.currentSlice);
set(handles.sliceSelector,'SliderStep',[ 1/(length(handles.linearRange)-1), 0.1]);
% handles = scatterSlice(handles);
% handles = plotSlice(handles);
handles.cellSelect = 1;

%plots everything!
handles = recruitFunction_Callback(handles.recruitFunction,eventdata,handles);

try
    %%%%light sheet movie
    load( [handles.posterHome,fs,'Data',fs,dataFolder,fs,'aviData',fs,'sliceAvi.mat']);
    handles.currentFrame = 1;
    handles.lsMovieSlice = 2;
    global CURRENT_MOVIE_SLICE 
    CURRENT_MOVIE_SLICE = handles.lsMovieSlice;

    handles.videoData = videoData;
    handles.frameRate = 15;
    imshow(handles.videoData(:,:,1,2),'Parent',handles.moviePlayer);
    handles.moviePlayer.Visible = 'off';
    colormap(handles.moviePlayer,'jet');
catch
    warning('Could not load sliceAvi.mat for this Larva');
end
%%%%%%%%%%%%%%%%%%%%


% UIWAIT makes DynamicPoster_DC_SFN2017 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function handles = scatterSlice(handles)
%plots the slice specific scatter plot and sets the current cell based on
%cellSelector slide

heatColor = handles.heatColor;
% heatBins = linspace(0,120,9);
% heatBins = linspace(0,6,9);
% heatBins = [-inf heatBins inf];
% handles.recruitMapData.f0DelayStacks =2;
xy = handles.recruitMapData.trialXYData;
h_xy = handles.recruitMapData.trialHistogramData;
% updateScatter=false;

mRecruitLvl = calculateRecruitment(handles);

cellRecruitment = handles.cellRecruitment(:,handles.trial2use);
alxPos = handles.spPos(handles.recruitMapData.alx,:);

sK = handles.currentSlice;
RCSegmentation = handles.linearRange;

% inSlice = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & cellRecruitment;
inSlice = findInSlice(handles,cellRecruitment);
nCells = sum(inSlice);
if nCells > 0
    
    set(handles.cellSelector,'Max',nCells);
    set(handles.cellSelector,'Min',1);
    set(handles.cellSelector,'SliderStep',[1/nCells, 5*(1/nCells)]);
    current_selector_value =  get(handles.cellSelector,'Value');
    if current_selector_value > nCells
        set(handles.cellSelector,'Value', round(nCells/2) );
        handles.cellSelect = round(nCells/2);
    end

else
    axes(handles.ScatterPlotAx);
    cla;
    return;
end
% handles.cellSelect = round(nCells/2);

trialColors = hsv(numel(handles.dF0 ));
f_inSlice = find(inSlice);
% randtoSelect = randperm( sum(inSlice) );
cellSelect = f_inSlice( handles.cellSelect);

%updates circled point on coronal slice plot
axes(handles.sliceAx);
if isfield(handles,'circlePoint')
       delete(handles.circlePoint );
end
handles.circlePoint = plot (alxPos(cellSelect,2),alxPos(cellSelect,3),'ro','linewidth',2);

%updates dFF plot with current cell
axes(handles.dffPlot), cla %this bit of code will generate a warning if the cell has changed but trial has not. not important. 
% yyaxis left
hold on
for k=1:numel(handles.dF0)
    if k==handles.trial2use
        plot(handles.stackX( xy(k).validBaseline > 0), handles.dF0{k}(cellSelect,xy(k).validBaseline > 0),'color',trialColors(k,:),'linewidth',2);
%     else
%         plot(handles.dF0{k}(cellSelect,logical( xy(k).validBaseline )),'k--', 'linewidth',0.5);
    end
end
% xlim([0 size(handles.dF0{1},2)]);
xlim([0 handles.stackX(end)]);
ylim([-0.1 1]);
set(gca,'XColor',[1,1,1],'YColor',[1,1,1]);
xlabel('Time (s)');
ylabel('\Delta F/F');
title('Fluorescence Time Series', 'color', [1 1 1]);
set(handles.dffPlot,'FontUnits','normalized','FontSize',0.06);
handles.dffPlot = gca;
% legend(h,'\Delta F/F');

% axes(handles.coloc_axes), cla,
% hold on
% % plot(handles.recruitmentMap.Amplitude(:,handles.trial2use),handles.recruitmentMap.Frequency(:,handles.trial2use),'k.');
% % plot(handles.recruitmentMap.Amplitude(cellSelect,handles.trial2use),handles.recruitmentMap.Frequency(cellSelect,handles.trial2use),'r.');
% plot(handles.recruitmentMap.Amplitude(:,handles.trial2use),handles.recruitmentMap.Correlation(:,handles.trial2use),'k.');
% plot(handles.recruitmentMap.Amplitude(cellSelect,handles.trial2use),handles.recruitmentMap.Correlation(cellSelect,handles.trial2use),'r.');
% xlabel('Amplitude Gain');
% % ylabel('Frequency Gain');
% ylabel('Correlation');
% currentXLim = get(gca,'Xlim');
% currentYLim = get(gca,'YLim');
% TT = handles.costeThreshold(handles.trial2use,:);
% 
% plot([TT(1) TT(1)],currentYLim,'g--');
% plot(currentXLim,[TT(2) TT(2)],'g--');
% axis([-0.01 0.05 0 1])
% set(gca,'XColor',[1 1 1], 'YColor', [1,1,1]);

%updates gain plot
axes(handles.ScatterPlotAx), cla
jxOff = zeros(1,numel(handles.dF0));
% handles.ScatterPlotAx.Legend = [];
legend('off');
hold on
for k=handles.trial2use
    
%     dFFindexDat = xy(k).xDataStacks;
%     xDat = xy(k).amplitude(dFFindexDat-handles.recruitMapData.f0DelayStacks);
    dFFData = handles.dF0{k}(cellSelect,:);
    
%     nBins = length(h_xy(k).amplitudeBins) - 1;
    if strcmp(handles.recruitFunctionX,'Amplitude') || strcmp(handles.recruitFunctionX,'Age') 
        
        nBins = length(h_xy(1).amplitudeBins) - 1;
%         jxOff = linspace(-0.1,0.1,numel(handles.dF0) );
        
        aIdx = h_xy(k).amplitudeIdx;
        xD = h_xy(k).amplitudeBins(2:end-1) + diff(h_xy(k).amplitudeBins(2:3) );
        xD = [xD xD(end)+xD(1)];
        xlabel('Integrated Amlitude (20 ms)');
        cF = handles.recruitMapData.cellFit(k).Amplitude(cellSelect,:);
%         xD = [xD xD(end)+xD(1)];
        xDat = xy(k).amplitude;
%     
    elseif strcmp(handles.recruitFunctionX,'Frequency')
        
        nBins = length(h_xy(1).frequencyBins) - 1;
%         jxOff = linspace(-0.1,2,numel(handles.dF0) );
        
        aIdx = h_xy(k).frequencyIdx';
        xD = h_xy(k).frequencyBins(1:end-1) + h_xy(k).frequencyBins(2)/2;
        cF = handles.recruitMapData.cellFit(k).Frequency(cellSelect,:);
%         xD = [xD xD(end)+xD(1)];
        xlabel('Frequency (Hz)');
        xDat = xy(k).frequency;
        
    elseif strcmp(handles.recruitFunctionX,'Correlation')
        title('Cross-Correlation of \DeltaF/F to Ventral Root');
        [crossCross,lags] = xcorr(dFFData(xy(k).validBaseline>0),xy(k).amplitude(xy(k).validBaseline>0),'coeff');
        plot(lags,crossCross,'color',trialColors(k,:),'linewidth',2);
        xlabel('Lag (Stacks)');
        ylabel('Correlation');
        text(lags(20),0.9,sprintf('Max Correlation %0.2f at lag %2.0f',max(crossCross), lags(crossCross == max(crossCross))),'FontSize',6);
        return
    
%     elseif strcmp(handles.recruitFunctionX,'Age')
%         return
    end
    
            
    for jj=1:nBins
        
%         if strcmp(handles.recruitFunctionX,'Amplitude')
            yIdx = find( [aIdx == jj] & xy(k).validBaseline );
%         else
%             fqIdx = find( [aIdx == jj] & xy(k).validBaseline );
%             yIdx = xy(k).xDataStacks(fqIdx);
%         end
        try
        dffBind(jj) = mean([dFFData(yIdx) + dFFData(yIdx-1)]./2 ); %average two bins as in base analysis
        dffBinError(jj) = std([dFFData(yIdx) + dFFData(yIdx-1)]./2 );
        catch
            disp('Got you!');
        end
    end
    
    plot(xDat(xy(k).validBaseline>0), dFFData(xy(k).validBaseline>0), 'k.');
    plot(xD+jxOff(k),dffBind,'ko','markerfacecolor',trialColors(k,:))
    hold on, 
    errorbar(xD+jxOff(k),dffBind,dffBinError,'.','color',trialColors(k,:));
    
    plot(xD, cF(1) + cF(2)*xD, 'color', trialColors(k,:));
    
    
    
    text(0.5,.9,sprintf('Gain of %3.4f per unit',mRecruitLvl(cellSelect)  ),'FontUnits','normalized','FontSize',0.05 );
    axis tight
    ylim([-0.10 1.05]);
    
end
title('Neuronal Gain Function','color',[1,1,1]);
set(gca,'XColor',[1,1,1],'YColor',[1,1,1]);
% xlabel('Integrated Amlitude (2 ms)');
ylabel('\Delta F/F');
set(gca,'FontUnits','normalized','FontSize',0.05);
handles.ScatterPlotAx = gca;


function mRecruitLevel = calculateRecruitment(handles)
%Calculates the recruitment for each cell based on recruitmentFunction
trial2use = handles.trial2use;
if strcmp(handles.recruitFunctionX,'Amplitude')
%    mRecruitLevel = nanmean(handles.recruitmentMap.Amplitude(:,trial2use),2);
    mRecruitLevel = handles.recruitmentMap.Amplitude(:,trial2use);
   
elseif strcmp(handles.recruitFunctionX,'Frequency')
%     mRecruitLevel = nanmean(handles.recruitmentMap.Frequency(:,trial2use),2);
    mRecruitLevel = handles.recruitmentMap.Frequency(:,trial2use);
    
elseif strcmp(handles.recruitFunctionX,'Correlation')
    mRecruitLevel = handles.recruitmentMap.Correlation(:,trial2use);
%      mRecruitLevel = handles.recruitmentMap.xCorrLag(:,trial2use);

elseif strcmp(handles.recruitFunctionX,'Age')
    mRecruitLevel = handles.alxAge;
    
end


function handles=plotSlice(handles)
%plots a recruitment map onto sliceAx given the current slice selection
%NOTE: CORRECT SLICE IS K:K+1 AND NUMBER OF SLICES IS
%LENGTH(RCSEGMENTATION)- 1

heatColor = handles.heatColor;
heatBins = handles.heatBins;

handles.cellClicker = handles.sliceAx.ButtonDownFcn;

mRecruitLvl = calculateRecruitment(handles);

cellRecruitment = handles.cellRecruitment(:,handles.trial2use);

alxPos = handles.spPos(handles.alx,:);

axes(handles.sliceAx), cla
sK = handles.currentSlice;
RCSegmentation = handles.linearRange;
all_inSlice = handles.spPos(:,1) >  RCSegmentation(sK-1) & handles.spPos(:,1) < RCSegmentation(sK) & ~handles.recruitMapData.alx;
allDots = handles.spPos(all_inSlice,:);

% inSlice = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & cellRecruitment;
inSlice = findInSlice(handles,cellRecruitment & ~isnan(mRecruitLvl));
% inSlice_notRecruited = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & ~cellRecruitment;
inSlice_notRecruited = findInSlice(handles,~cellRecruitment);

[histo,Edge,ActBin] = histcounts(mRecruitLvl(inSlice),heatBins);

handles.inSlice = inSlice;
posRecruit = alxPos(inSlice,:);
posNeutral = alxPos(inSlice_notRecruited,:);
    
plot(allDots(:,2),allDots(:,3),'.','color',[0.2 0.2 0.2],'hittest','off');
hold on
scatter(posRecruit(:,2),posRecruit(:,3),32,heatColor(ActBin,:),'o','filled','hittest','off');
colormap(jet);
axis equal;% axis tight
set(gca,'TickLength',[0,0],'XTick',[],'YTick',[],'FontUnits','normalized','FontSize',0.05);
% 
% cb = colorbar;
% cb.Color = [1,1,1];
% ylabel(cb,'\DeltaF/F');
% 
% % cbTicks = get(cb,'Ticks');
% set(cb,'Ticks',[0:0.2:1]);
% tickLabelsVector = linspace(heatBins(2),heatBins(end-1),length(get(cb,'Ticks'))) ;
% for kjk=1:length(tickLabelsVector)
%    tickLabelsS{kjk,1} = sprintf('%2.4f',tickLabelsVector(kjk)); 
% end
% set(cb,'TickLabels',tickLabelsS);

% if get(handles.draw_AlxNonRecruit,'Value')
%     plot(posNeutral(:,2),posNeutral(:,3),'kx');
% end
title(sprintf('Slice %2.0f',sK-1),'color',[1,1,1]);

handles.sliceAx.ButtonDownFcn = handles.cellClicker;
%%%%%
axes(handles.ageHistogramAx), cla
alxRecruit = mRecruitLvl(inSlice);
alxAge = handles.alxAge(inSlice);
hold on
h1=histogram(alxRecruit(~alxAge), 20 );
h2=histogram(alxRecruit(alxAge),20);
h1.EdgeColor = 'none'; h2.EdgeColor = 'none';
legend([h1,h2],{'Younger Chx10 Neurons', 'Older Chx10 Neurons'} );
ylabel('Cell Count');
xlabel(sprintf('%s Metric',handles.recruitFunctionX) );
axis tight

%%%%%
axes(handles.slicePosMap), cla
hold on

plot(handles.spPos(:,1),handles.spPos(:,2),'k.');
plot(handles.alxPos(:,1),handles.alxPos(:,2),'y.');
axis equal
axis tight

currentLim = get(gca,'YLim');

if isfield(handles,'SliceMap_gLine')
   for gl = 1:numel(handles.SliceMap_gLine)
       delete(handles.SliceMap_gLine{gl} );
   end
end
axis equal
axis tight

handles.SliceMap_gLine{1} = plot([RCSegmentation(sK),RCSegmentation(sK)], currentLim, 'g');
handles.SliceMap_gLine{2} = plot([RCSegmentation(sK-1),RCSegmentation(sK-1)], currentLim, 'g');


function inSlice = findInSlice(handles,cellRecruit)
%consisent framework for identifying which cells are in slice
alxPos = handles.spPos(handles.alx,:);
sK = handles.currentSlice;
RCSegmentation = handles.linearRange;
inSlice = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & cellRecruit;

% --- Outputs from this function are returned to the command line.
function varargout = DynamicPoster_DC_SFN2017_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on slider movement.
function sliceSelector_Callback(hObject, eventdata, handles)
% hObject    handle to sliceSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.currentSlice = round(  get(hObject,'Value') ) +1;
handles = plotSlice(handles);
handles = scatterSlice(handles); 
guidata(hObject,handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliceSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function handles = cellSelector_Callback(hObject, eventdata, handles)
% hObject    handle to cellSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.cellSelect = round( get(hObject,'Value') );
handles = scatterSlice(handles);
guidata(hObject,handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function cellSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cellSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% 
% % --- Executes on button press in draw_AlxNonRecruit.
% function draw_AlxNonRecruit_Callback(hObject, eventdata, handles)
% % hObject    handle to draw_AlxNonRecruit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles=guidata(hObject);
% handles = plotSlice(handles);
% guidata(hObject,handles);
% % Hint: get(hObject,'Value') returns toggle state of draw_AlxNonRecruit


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over draw_AlxNonRecruit.
function draw_AlxNonRecruit_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to draw_AlxNonRecruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in recruitFunction.
function handles = recruitFunction_Callback(hObject, eventdata, handles)
% hObject    handle to recruitFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
string=get(hObject,'String');
value = get(hObject,'Value');
trialColors = hsv(numel(handles.dF0 ));
% handles.heatColor = jet(10);

% axes(handles.motorDriveAx); 
axes(handles.dffPlot); yyaxis right
hold on
switch string{value}
    case 'Amplitude'
        cla;
        handles.recruitFunctionX = 'Amplitude';
        heatBins = handles.plotAmplitudeBins;
        handles.heatBins = [-inf heatBins inf];
        handles.heatColor = jet(10);
        
        xy = handles.recruitMapData.trialXYData;
        hold on
        for k=1:numel(handles.dF0)
            s = smooth(xy(k).amplitude);
            if k==handles.trial2use
               h = plot(handles.stackX,s,'k','linewidth',1);
%             else
%                 plot(s,'k--','linewidth',0.5);
            end
        end
        xlabel('Time (s)');
        ylabel('Max. Integrated Amplitude (20 ms)');
        title('Ventral Root Signal');
        axis tight

    case 'Frequency'
        cla;
        handles.recruitFunctionX = 'Frequency';
        heatBins = handles.plotFrequencyBins;
        handles.heatBins = [-inf heatBins inf];
        handles.heatColor = jet(10);
        
        xy = handles.recruitMapData.trialXYData;
         
        hold on
        for k=1:numel(handles.dF0)
            
            s = smooth( smooth(xy(k).frequency) );
            if k==handles.trial2use
                h = plot(handles.stackX,s,'color','k','linewidth',1);
%             else
%                 plot(s,'k--','linewidth',0.5);
            end
        end
        xlabel('Time (s)');
        ylabel('Max Frequency (Hz)');
        title('Ventral Root Signal');
        axis tight
        
    case 'Correlation'
        cla;
        handles.recruitFunctionX = 'Correlation';
        heatBins = handles.plotCorrelationBins;
        handles.heatBins = [-inf heatBins inf];
        handles.heatColor = jet(10);
        
        xy = handles.recruitMapData.trialXYData;
        hold on
        for k=1:numel(handles.dF0)
            s = smooth(handles.stackX,xy(k).amplitude);
            if k==handles.trial2use
                h = plot(handles.stackX,s,'k','linewidth',1);
%             else
%                 plot(s,'k--','linewidth',0.5);
            end
        end
        xlabel('Time (s)');
        ylabel('Max. Integrated Amplitude (20 ms)');
        title('Ventral Root Signal');
        axis tight
        
    case 'Age'
        handles.recruitFunctionX = 'Age';
        heatBins = handles.plotAgeBins;
        handles.heatBins = [-inf heatBins inf];
        handles.heatColor = [0 1 0; 1 0 0];
end
handles.dffPlot = gca;
axes(handles.dffPlot);
set(gca,'ycolor',[1,1,1]);
yyaxis left
handles = plotSlice(handles);
handles = scatterSlice(handles);

if handles.update_fts_legend
    axes(handles.dffPlot);
    % legend(h,'Ventral Root');
    [ab,L] = legend({'\DeltaF/F','Ventral Root'},'Location','NorthEast','Box','off');%,'Position', [0.725   0.111  0.0623    0.0346]);
    hL=findobj(L,'type','line');  % get the lines, not text
    currentLineLength  = get(hL,'XData');
%     pushforward = diff(currentLineLength{1})/2;
    set(hL([1,3],1),'linewidth',3); %%%'XData',[currentLineLength{1}(1)+pushforward currentLineLength{1}(2)]) ;
    cBoxPos = get(ab,'Position');
    set(ab,'FontSize',5,'Position',[cBoxPos(1)+0.00 cBoxPos(2)+0.03 cBoxPos(3) cBoxPos(4)]);
    handles.update_fts_legend = false;
end
% currbox = get(ab,'Position');
% set(ab,'Position', [currbox(1)+.03, currbox(2)+.03, currbox(3), currbox(4)]);


guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns recruitFunction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from recruitFunction


% --- Executes during object creation, after setting all properties.
function recruitFunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recruitFunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Dataset_Selection.
function Dataset_Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Dataset_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles=loadDataset(handles,eventdata);
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns Dataset_Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Dataset_Selection


% --- Executes during object creation, after setting all properties.
function Dataset_Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dataset_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in colorbar_increase.
function colorbar_increase_Callback(hObject, eventdata, handles)
% hObject    handle to colorbar_increase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

string =get(handles.recruitFunction,'String');
value = get(handles.recruitFunction,'Value');

% axes(handles.motorDriveAx); 
switch string{value}
    case 'Frequency'
    freqMax = handles.plotFrequencyBins(end);
    handles.plotFrequencyBins = linspace(0, freqMax*0.9,9);
    case 'Amplitude'
    ampMax = handles.plotAmplitudeBins(end);
    handles.plotAmplitudeBins = linspace(0, ampMax*0.9, 9);
    case 'Correlation'
	corrMax = handles.plotCorrelationBins(end);
    corrMin = handles.plotCorrelationBins(1);
    handles.plotCorrelationBins = linspace(corrMin, corrMax*0.9,9);        
end

handles = recruitFunction_Callback(handles.recruitFunction,eventdata,handles);

guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of colorbar_increase


% --- Executes on button press in colorbar_decrease.
function colorbar_decrease_Callback(hObject, eventdata, handles)
% hObject    handle to colorbar_decrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

string =get(handles.recruitFunction,'String');
value = get(handles.recruitFunction,'Value');

% axes(handles.motorDriveAx); 
switch string{value}
    case 'Frequency'
    freqMax = handles.plotFrequencyBins(end);
    handles.plotFrequencyBins = linspace(0, freqMax*1.1,9);
    case 'Amplitude'
    ampMax = handles.plotAmplitudeBins(end);
    handles.plotAmplitudeBins = linspace(0, ampMax*1.1, 9);
    case 'Correlation'
	corrMax = handles.plotCorrelationBins(end);
    corrMin = handles.plotCorrelationBins(1);
    handles.plotCorrelationBins = linspace(corrMin, corrMax*1.1,9);        
end

handles = recruitFunction_Callback(handles.recruitFunction,eventdata,handles);

guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of colorbar_decrease


% --- Executes on slider movement.
function trialSelector_Callback(hObject, eventdata, handles)
% hObject    handle to trialSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
handles.trial2use = mod(handles.trial2use, numel(handles.dF0) ) + 1;
handles.update_fts_legend = true;

handles = recruitFunction_Callback(handles.recruitFunction,eventdata,handles);

% axes(handles.triggAvgHeatMap); colormap hot;

guidata(hObject,handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function trialSelector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in playMovie.
function playMovie_Callback(hObject, eventdata, handles)
% hObject    handle to playMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CURRENT_MOVIE_SLICE 
handles = guidata(hObject);
% playBinary = get(hObject,'Value');

while get(hObject,'Value')
% tic
    vidFrame = handles.currentFrame;
%     zSlice = handles.lsMovieSlice;
    imshow(handles.videoData(:,:,vidFrame,CURRENT_MOVIE_SLICE),'Parent',handles.moviePlayer);
    colormap(handles.moviePlayer,'jet');
%     handles.moviePlayer.Visible = 'off';
%     drawnow;
    handles.currentFrame = mod(handles.currentFrame,870) + 1;
% toc;
    pause(1/handles.frameRate);
end
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of playMovie


% --- Executes on button press in startSpinAnimation.
function startSpinAnimation_Callback(hObject, eventdata, handles)
% hObject    handle to startSpinAnimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of startSpinAnimation
handles=guidata(hObject);

while get(hObject,'Value')
% tic
    if ~hasFrame(handles.ls_videoReader)
        handles.ls_videoReader.CurrentTime=0;
    end
    vidFrame = readFrame(handles.ls_videoReader);
    imshow(vidFrame,'Parent',handles.shockSpinner);
% toc;
    pause(1/handles.ls_videoReader.FrameRate);
end

guidata(hObject,handles);


% --- Executes on button press in toggleSlice.
function toggleSlice_Callback(hObject, eventdata, handles)
% hObject    handle to toggleSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CURRENT_MOVIE_SLICE 
handles = guidata(hObject);
CurrZSlice = handles.lsMovieSlice;
handles.lsMovieSlice = mod(CurrZSlice,size(handles.videoData,4))+1;
CURRENT_MOVIE_SLICE = handles.lsMovieSlice;
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of toggleSlice

% --- Executes on mouse press over axes background.
function sliceAx_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to sliceAx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

position_list = handles.alxPos(handles.inSlice,2:3);
clickPosition = eventdata.IntersectionPoint(1:2);
distance_to_list = pdist2(clickPosition,position_list);

[mD,cellSelected] = min(distance_to_list);
handles.cellSelector.Value = cellSelected;

handles = cellSelector_Callback(handles.cellSelector,eventdata,handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

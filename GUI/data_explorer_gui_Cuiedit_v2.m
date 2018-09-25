function varargout = data_explorer_gui_Cuiedit_v2(varargin)
% DATA_EXPLORER_GUI_CUIEDIT_V2 MATLAB code for data_explorer_gui_Cuiedit_v2.fig
%      DATA_EXPLORER_GUI_CUIEDIT_V2, by itself, creates a new DATA_EXPLORER_GUI_CUIEDIT_V2 or raises the existing
%      singleton*.
%
%      H = DATA_EXPLORER_GUI_CUIEDIT_V2 returns the handle to a new DATA_EXPLORER_GUI_CUIEDIT_V2 or the handle to
%      the existing singleton*.
%
%      DATA_EXPLORER_GUI_CUIEDIT_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATA_EXPLORER_GUI_CUIEDIT_V2.M with the given input arguments.
%
%      DATA_EXPLORER_GUI_CUIEDIT_V2('Property','Value',...) creates a new DATA_EXPLORER_GUI_CUIEDIT_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before data_explorer_gui_Cuiedit_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to data_explorer_gui_Cduiedit_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help data_explorer_gui_Cuiedit_v2

% Last Modified by GUIDE v2.5 24-Sep-2018 15:51:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @data_explorer_gui_Cuiedit_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @data_explorer_gui_Cuiedit_v2_OutputFcn, ...
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

% --- Executes just before data_explorer_gui_Cuiedit_v2 is made visible.
function data_explorer_gui_Cuiedit_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to data_explorer_gui_Cuiedit_v2 (see VARARGIN)

% Choose default command line output for data_explorer_gui_Cuiedit_v2
handles.output = hObject;

handles.sliceThickness = 20;

% homepath = getenv('HOME');
% fs = filesep;
% handles.guiHome_functions = [homepath,fs,'Dropbox',fs,'CellAnalysisUtilities',fs,'LightSheet' fs 'light_sheet_data_explorer_gui', fs, 'functions'];
% addpath(handles.guiHome_functions);
handles.stackAcqFreq = 1; %initialized with a dummy value

handles.cellSelect = 1; %for plotting on slice axis
handles.cellSelect_idx = 1; %for keeping track of absolute cell identity

handles.microns_per_z = 5;

handles.roi = struct;
handles.currView = [];
handles.zStack = [];
handles.calling_function = 'None';

handles.sPMap_Ax_roi = [];
handles.cell_roi_list = [];
handles.posmap_img = [];

handles.im=[];

handles.timeListener=addlistener(handles.sliceSelector,'ContinuousValueChange',@sliceSelector_Callback);
handles.zListener=addlistener(handles.cellSelector,'ContinuousValueChange',@cellSelector_Callback);

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using data_explorer_gui_Cuiedit_v2.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

function handles=loadDataset(filepath,handles,~)
%switches data set
load(filepath);
assignin('base','filepath',filepath);
handles.fts = fluorescence_time_series;
assignin('base','fluorescence_time_series',fluorescence_time_series);
assignin('base','spPos',spPos);
handles.spPos = spPos;
handles.Sc = Sc;

if exist('spRadiiXYZ','var')
    handles.spRadiiXYZ = spRadiiXYZ;
else
    handles.spRadiiXYZ = NaN;
end

global AllPosition %Defined a global here!!!!!
AllPosition=spPos;


if exist('dFF','var')
    handles.dFF = dFF;
    assignin('base','dFF',dFF);
else
    handles.dFF = NaN;
end

if exist('roi','var')
    handles.roi=roi;
    handles = display_roiListbox(handles);
else
    handles.roi=[];
    handles.roi.name = 'default';
    handles.roi.members = ones(size(spPos,1),1,'logical');
    handles.roi.ROIcellNo = [1:size(spPos,1)]';
    handles = display_roiListbox(handles);
    update_roiMasterList(handles);
end



% handles.roi = zeros(size(spPos,1),1,'logical');

% handles.linearRange = floor( min(handles.spPos(:,1)) ):handles.sliceThickness: ceil( max(handles.spPos(:,1)) );
% handles.linearRange = [handles.linearRange inf];
handles = get_linear_range(handles);
handles.currentSlice = ceil(length(handles.linearRange)/2);

set(handles.sliceSelector,'Max',length(handles.linearRange)-1);
set(handles.sliceSelector,'Value',handles.currentSlice);
set(handles.sliceSelector,'SliderStep',[ 1/(length(handles.linearRange)-1), 0.1]);

handles.labelVector = zeros(size(handles.spPos,1),1,'logical');
if size(handles.spPos,1)>10000
    handles.labelVector(1:10:end) = 1;
else
    handles.labelVector(1:end) = 1;
end

handles.ColorOfDots= max(handles.fts');%set color of dots to each cells' highest fluo value.
handles.caxis0.String = num2str( min(handles.ColorOfDots), '%4.0f' );
handles.caxis1.String = num2str( max (handles.ColorOfDots), '%4.0f' );

handles = plot_slice_maps(handles);
handles = plot_pos_maps(handles);
handles = plot_fts_in_slice(handles);


function handles = get_linear_range(handles)
% handles.linearRange = floor( min(handles.spPos(:,1))
% ):handles.sliceThickness: ceil( max(handles.spPos(:,1)) ); %Cui changed
% this line
% handles.linearRange =[0:handles.sliceThickness:ceil(max(handles.spPos(:,1))/handles.sliceThickness)*handles.sliceThickness];%Cui
% Changed this line at July 4th 2018
handles.linearRange =[min(handles.spPos(:,1)):handles.sliceThickness:min(handles.spPos(:,1))+ceil((max(handles.spPos(:,1))-min(handles.spPos(:,1)))/handles.sliceThickness)*handles.sliceThickness];
% assignin('base','linearRange',handles.linearRange);
handles.linearRange = [handles.linearRange inf];

function handles = set_label_density(handles, density)
%density is a fraction between 0 and 1
spacing = round(density^-1);
handles.labelVector = zeros(size(handles.spPos,1),1,'logical');
handles.labelVector(1:spacing:end) = 1;

function set_slice_thickness_Callback(hObject, ~, handles)
% hObject    handle to set_slice_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_slice_thickness as text
%        str2double(get(hObject,'String')) returns contents of set_slice_thickness as a double
handles.sliceThickness = str2num(get(hObject,'String') );
handles = get_linear_range(handles);

handles.currentSlice = ceil(length(handles.linearRange)/2);

set(handles.sliceSelector,'Max',length(handles.linearRange)-1);
set(handles.sliceSelector,'Value',handles.currentSlice);
set(handles.sliceSelector,'SliderStep',[ 1/(length(handles.linearRange)-1), 0.1]);

handles = plot_slice_maps(handles);
handles = plot_pos_maps(handles);
handles = plot_fts_in_slice(handles);
guidata(hObject, handles);
% handles = plot_third_dim(handles);


%%%%%%%%%%%%%%%%%%%%

% function handles = plot_third_dim(handles)
% %updates the 3D plot with spots specified in handles.labelVector. By default  this
% %will be a 1/100 sample of the data on loading. 
% posB = handles.spPos(handles.labelVector & ~handles.inSlice,:);
% posR = handles.spPos(handles.labelVector & handles.inSlice,:);
% axes(handles.plot_3Dim); cla
% hold on
% plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
% plot3(posR(:,1),posR(:,2),posR(:,3),'r.');
% grid on
% axis equal

% UIWAIT makes data_explorer_gui_Cuiedit_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function handles = plot_fts_in_slice(handles)
%plots the slice specific scatter plot and sets the current cell based on
%cellSelector slide

% inSlice = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & cellRecruitment;
if handles.PlotSelect.Value < 3
    handles.inSlice = findInSlice(handles);
    inSlice = handles.inSlice;
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

    f_inSlice = find(inSlice);
    % randtoSelect = randperm( sum(inSlice) );
    cellSelect = f_inSlice( handles.cellSelect);
    handles.cellSelect_idx = cellSelect;
    %updates circled point on coronal slice plot
    axes(handles.sliceAx);
    if isfield(handles,'CirclePointinSlice')
           delete(handles.CirclePointinSlice );
    end
    hold on;

    handles.CellNoToFindinSlice=cellSelect;
    handles=DrawCirclePointinSlice(handles);
    handles=plotfts(handles,cellSelect);
else
    
%     axes(handles.dffPlot); hold on
    if isfield(handles,'timeLine') && ishandle(handles.timeLine)
%         set(handles.timeLine, 'XData', [handles.sliceSelector.Value, handles.sliceSelector.Value]);
          handles.timeLine.XData = [handles.sliceSelector.Value, handles.sliceSelector.Value];
          handles.timeLine.Visible = 'on';
    else
        ylim = get(handles.dffPlot, 'YLim');
        handles.timeLine = plot([handles.sliceSelector.Value, handles.sliceSelector.Value], ylim, 'r-', 'Parent', handles.dffPlot);
        handles.timeLine.Visible = 'off';
%         disp('Setting new time line');
        pause(0.01);
    end
%     hold off
end


function handles=DrawCirclePointinSlice(handles)
% axes(handles.sliceAx);

if isfield(handles,'CirclePointinSlice')
    delete(handles.CirclePointinSlice);
end

handles.CirclePointinSlice=plot3(handles.spPos(handles.CellNoToFindinSlice,1), handles.spPos(handles.CellNoToFindinSlice,2),handles.spPos(handles.CellNoToFindinSlice,3),'bo','linewidth',3,'Markersize',10,'Parent',handles.sliceAx);
% guidata(hObject,handles);
% assignin('base','circinslice',handles.CirclePointinSlice);

function handles=DrawCirclePoint(handles)
% axes(handles.PosMap);

if isfield(handles,'CirclePoint')
    delete(handles.CirclePoint);
end
handles.CirclePoint=plot3(handles.spPos(handles.CellNoToFind,1), handles.spPos(handles.CellNoToFind,2),handles.spPos(handles.CellNoToFind,3),'wo','linewidth',3,'Markersize',10);

function handles=plotfts(handles, cellSelect)
%plots fluorescent time series on the time series plot
axes(handles.dffPlot), cla %this bit of code will generate a warning if the cell has changed but trial has not. not important. 
% yyaxis left
hold on

% plot(handles.fts(cellSelect,:),'color',[0 1 0],'linewidth',2);
if get(handles.raw_trace_selector,'Value')
    plotTraces_GUI(handles.fts(cellSelect,:), 200, 2, handles.dffPlot);
    colorbar('off');
    ylabel('Intensity - Mean');
elseif get(handles.dFF_traces,'Value')
    plotTraces_GUI(handles.dFF(cellSelect,:), 1, 3, handles.dffPlot);
    colorbar('off');
    ylabel('\Delta F/F');
elseif get(handles.heatmap_selector,'Value')
    toplot = handles.fts(cellSelect,:);
    toplot = bsxfun(@minus, toplot, median(toplot,2));
    imagesc(handles.dffPlot, toplot );
    axis tight
    colorbar;  %colorbar('on');
    colormap jet;
elseif get(handles.dFF_heatmap,'Value')
    toplot = handles.dFF(cellSelect,:);
    imagesc(handles.dffPlot, toplot );
    axis tight
    colorbar;  %colorbar('on');
    colormap jet;
end

xlabel('Frames');
title('Fluorescence Time Series', 'color', [1 1 1]);
set(handles.dffPlot,'FontUnits','normalized','FontSize',0.06);
% handles.dffPlot = gca;

function handles = plot_pos_maps(handles)
%plots a recruitment map onto sliceAx given the current slice selection
%NOTE: CORRECT SLICE IS K:K+1 AND NUMBER OF SLICES IS
%LENGTH(RCSEGMENTATION)- 1

%%%%%
%%%%%

if handles.PlotSelect.Value==1
    axes(handles.slicePosMap), cla
    hold on
    % posB = handles.spPos(handles.labelVector & ~handles.inSlice,:);
    % posR = handles.spPos(handles.labelVector & handles.inSlice,:);
    % plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
    % plot3(posR(:,1),posR(:,2),posR(:,3),'r.');
    %----Cui Edit----
    posAllCell=handles.spPos;
    scatter3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),10,handles.ColorOfDots,'.','hittest','off', 'Parent', handles.slicePosMap);
    colormap jet;
    caxis([str2double(handles.caxis0.String),str2double(handles.caxis1.String)])
    CB1=colorbar;
    CB1.Color='w';
    
    % CB1.Label.String='Max Fluorescence Intensity of Cells (Dots)';
    posR_C = handles.spPos( handles.inSlice,:);
    plot3(posR_C(:,1),posR_C(:,2),posR_C(:,3),'k.','Parent', handles.slicePosMap);
    grid on;
    xlabel('X','Color','w');
    ylabel('Y','Color','w');
    zlabel('Z','Color','w');
    
    % set(gca,'color',[0 0 0]);
    %----Cui Edit End----
    grid on
    axis equal
    if get(handles.CellFinderStay,'Value')==1
        handles=DrawCirclePoint_Special(hObject,handles);
    end
elseif handles.PlotSelect.Value==2
    ShowAllROI=get(handles.ShowAllROICheck,'Value');
        current_roi = handles.roiMaster.Value;    
        axes(handles.slicePosMap), cla
        hold on
        posB = handles.spPos(handles.labelVector,:);
        if islogical(handles.roi(current_roi).members)==0
            handles.roi(current_roi).members=logical(handles.roi(current_roi).members);
        end
        
        if ShowAllROI==0
            CellNoToDisp=find(handles.roi(current_roi).members==1);
        elseif ShowAllROI==1
            CellNoToDisp=[];
            for i=1:size(handles.roi,2)
                CellNoToDisp=union(CellNoToDisp,find(handles.roi(i).members));
            end
        end
        assignin('base','CellNoToDisp',CellNoToDisp);
        assignin('base','spPos',handles.spPos);
        posR = handles.spPos(CellNoToDisp,:);
        if isfield(handles,'DispColor')==0
            ColorR='r';
        else
            ColorR =handles.DispColor(CellNoToDisp,:);
            assignin('base','ColorR',ColorR);
            assignin('base','posR',posR);
        end
        plot3(posB(:,1),posB(:,2),posB(:,3),'k.','Parent', handles.slicePosMap);
%         plot3(posR(:,1),posR(:,2),posR(:,3),'ro', 'markerfacecolor','r','markersize',3);
%         ColorBarlim = caxis;
        caxis([0 100]);
        CB2=colorbar;
        CB2.Color='w';
%         ColorList=(ColorR./100).*(ColorBarlim(2)-ColorBarlim(1))+ColorBarlim(1)
        scatter3(posR(:,1),posR(:,2),posR(:,3),500,ColorR,'.', 'Parent', handles.slicePosMap);
        grid on
        axis equal

        currentLim = get(gca,'YLim');
        zLim = get(gca,'ZLim');
        RCSegmentation = handles.linearRange;

        if isfield(handles,'SliceMap_gLine')
           for gl = 1:numel(handles.SliceMap_gLine)
               delete(handles.SliceMap_gLine{gl} );
           end
        end
        axis equal
        axis tight

        sK = handles.currentSlice;
        handles.SliceMap_gLine{1} = plot([RCSegmentation(sK),RCSegmentation(sK)], currentLim, 'g', 'Parent', handles.slicePosMap);
        handles.SliceMap_gLine{2} = plot([RCSegmentation(sK+1),RCSegmentation(sK+1)], currentLim, 'g', 'Parent', handles.slicePosMap);
        handles.SliceMap_gLine{3} = plot3([RCSegmentation(sK),RCSegmentation(sK)], currentLim, [zLim(2),zLim(2)], 'g', 'Parent', handles.slicePosMap);
        handles.SliceMap_gLine{4} = plot3([RCSegmentation(sK+1),RCSegmentation(sK+1)], currentLim, [zLim(2), zLim(2)], 'g', 'Parent', handles.slicePosMap);
        if get(handles.CellFinderStay,'Value')==1
            handles=DrawCirclePoint_Special(hObject,handles);
        end
elseif handles.PlotSelect.Value==3 && ~isempty(handles.currView)
    currentT =  round( get(handles.sliceSelector,'Value') );
    climz = [str2double(handles.caxis0.String), str2double(handles.caxis1.String)];
    
    if isnan(climz(1)) 
       frame1 = handles.currView(:,:,1);
       climz = [prctile(frame1(:), 5) prctile(frame1(:), 99)];
       handles.caxis0.String = num2str(climz(1), '%4.0f');
       handles.caxis1.String = num2str(climz(2), '%4.0f');   
    end

    cla(handles.slicePosMap);
    handles.posmap_img = imshow(handles.currView(:,:,currentT), climz, 'Parent', handles.slicePosMap);

    if isempty(handles.cell_roi_list) && sum(handles.roi.members) > 0
        roi_members_cell_no = find(handles.roi.members);
        handles.cell_roi_list = roi_members_cell_no(handles.roiListbox.Value);
    end
    
    if numel(handles.cell_roi_list) > 0 
        cxy = handles.spPos(handles.cell_roi_list,1:2)/handles.Sc(1,1);
        rxy = handles.spRadiiXYZ(handles.cell_roi_list,1)/handles.Sc(1,1);
        handles.sPMap_Ax_roi = viscircles(handles.slicePosMap, cxy, rxy, 'LineWidth', 0.5, 'EnhanceVisibility', 0);
    end
    
    title(handles.slicePosMap,sprintf('Stack %04.0f', currentT));
    drawnow;
    pause(0.005);
%     colormap(gca, 'jet');
%     caxis('auto');
end
set(gca,'TickLength',[0,0],'XTick',[],'YTick',[],'ZTick',[],'FontUnits','normalized','FontSize',0.05);

function handles = plot_slice_maps(handles)

if handles.PlotSelect.Value < 3
    set(handles.load_z_stack, 'Visible', 'off');
    set(handles.sliceAx,'Visible','on');
    handles.cellClicker = handles.sliceAx.ButtonDownFcn;
    axes(handles.sliceAx), cla
    sK = handles.currentSlice;
    % RCSegmentation = handles.linearRange;
    % all_inSlice = handles.spPos(:,1) >  RCSegmentation(sK-1) & handles.spPos(:,1) < RCSegmentation(sK);
    all_inSlice = findInSlice(handles);
    allDots = handles.spPos(all_inSlice,:);
    handles.inSlice = all_inSlice; 

    current_roi = get(handles.roiMaster,'Value');

    if get(handles.ShowAllROICheck,'Value')==0
        roi_in_slice = handles.roi(current_roi).members(all_inSlice);
    elseif get(handles.ShowAllROICheck,'Value')==1
        current_roi_cells=handles.roi(1).members;
        for i=1:size(handles.roi,2)
            current_roi_cells=current_roi_cells+handles.roi(i).members;
        end
        current_roi_cells(current_roi_cells>0)=1;
        roi_in_slice = current_roi_cells(all_inSlice);
    end
    
    if islogical(roi_in_slice)==0
        roi_in_slice =logical(roi_in_slice);
    end

    if isfield(handles,'DispColor')==0
        ColorR='r';
    else
        allColorinSlice=handles.DispColor(all_inSlice,:);
        ColorR =allColorinSlice(roi_in_slice,:);
        assignin('base','DispColorMatrix',handles.DispColor);
        assignin('base','roi_in_slice',roi_in_slice);
        assignin('base','ColorRinSlice',ColorR);
    end
    
    hold(handles.sliceAx,'on');
    plot3(allDots(roi_in_slice,1),allDots(roi_in_slice,2),allDots(roi_in_slice,3),'.','color',[1 0 0],'hittest','off','Parent',handles.sliceAx);
    plot3(allDots(~roi_in_slice,1),allDots(~roi_in_slice,2),allDots(~roi_in_slice,3),'.','color',[0.2 0.2 0.2],'hittest','off','Parent',handles.sliceAx);
    hold(handles.sliceAx,'off');
%     caxis([0 100]);
    

    set(gca,'TickLength',[0,0],'XTick',[],'YTick',[],'ZTick',[],'FontUnits','normalized','FontSize',0.05);
    grid on
    axis equal
    view(90,0);
    title(sprintf('Slice %2.0f',sK),'color',[1,1,1]);

    handles.sliceAx.ButtonDownFcn = handles.cellClicker;
    
elseif ~strcmp(handles.calling_function, 'sliceSelector')
    if ~isempty(handles.zStack)
       zLvl = round( get(handles.cellSelector, 'Value') );
       climz = [str2double(handles.caxis0.String), str2double(handles.caxis1.String)];
    if isnan(climz(1))
       imshow(handles.zStack(:,:,zLvl), [50 1000], 'Parent', handles.sliceAx); 
       caxis(handles.sliceAx, 'auto');
       climz = get(handles.sliceAx, 'CLim');
       handles.caxis0.String = num2str(climz(1), '%4.0f');
       handles.caxis1.String = num2str(climz(2), '%4.0f');       
    else
       cla(handles.sliceAx);
       imshow(handles.zStack(:,:,zLvl), climz, 'Parent', handles.sliceAx); 
    end
       title(handles.sliceAx,sprintf('Z=%3.0f',zLvl));
       inZ = round( handles.spPos(:,3)/handles.microns_per_z ) == zLvl;
       roi_in_slice = handles.roi.members & inZ;
       if sum(roi_in_slice) > 0
           cxy = handles.spPos(roi_in_slice,1:2)/handles.Sc(1,1);
           rxy = handles.spRadiiXYZ(roi_in_slice,1)/handles.Sc(1,1);
           handles.sAx_roi = viscircles(handles.sliceAx, cxy, rxy, 'LineWidth', 0.5, 'EnhanceVisibility', 0);
       end
       drawnow;
       pause(0.005);
%        handles.sliceAx = gca;
%        axis equal
    end
end

function inSlice = findInSlice(handles)
%consisent framework for identifying which cells are in slice
sK = handles.currentSlice;
RCSegmentation = handles.linearRange;
inSlice = handles.spPos(:,1) >  RCSegmentation(sK) & handles.spPos(:,1) < RCSegmentation(sK)+(RCSegmentation(2)-RCSegmentation(1));

% --- Outputs from this function are returned to the command line.
function varargout = data_explorer_gui_Cuiedit_v2_OutputFcn(hObject, eventdata, handles)
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
handles.calling_function = 'sliceSelector';
handles.currentSlice = round( get(hObject,'Value') );
% round( get(hObject,'Value') );
% handles.currentSlice 
% min=get(hObject,'Min')
% max=get(hObject,'Max')
handles = plot_slice_maps(handles);
handles = plot_pos_maps(handles);
handles = plot_fts_in_slice(handles); 
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
handles.calling_function = 'cellSelector';

if handles.PlotSelect.Value == 3
    plot_slice_maps(handles);
else
    handles.cellSelect = round( get(hObject,'Value') );
    handles = plot_fts_in_slice(handles);
end
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


% --- Executes on mouse press over axes background.
function sliceAx_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to sliceAx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

position_list = single( handles.spPos(handles.inSlice,:) );
clickPosition = eventdata.IntersectionPoint;
distance_to_list = pdist2(single( clickPosition ),position_list);

[mD,CellSelectedinSlice] = min(distance_to_list);
handles.cellSelector.Value = CellSelectedinSlice;

handles = cellSelector_Callback(handles.cellSelector,eventdata,handles);

set(handles.ClickCellNo,'String',num2str(handles.cellSelect_idx));

if size(handles.roi,2)>1
    for i=1:size(handles.roi,2)
        if handles.roi(i).members(handles.cellSelect_idx)==1
            set(handles.CateNo,'String',num2str(i));
        end
    end
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function data_set_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to data_set_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_data_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.Currentpath='0';
[currentfile,handles.Currentpath] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(handles.Currentpath);
handles.currentfile = [handles.Currentpath,currentfile];
handles = loadDataset(handles.currentfile,handles,eventdata);
set(handles.DataFileName, 'String', [handles.Currentpath,currentfile]);
guidata(hObject,handles);


% --- Executes on selection change in roiListbox.
function roiListbox_Callback(hObject, eventdata, handles)
% hObject    handle to roiListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiListbox
% CurrentString=get(handles.roiListbox,'String');
handles.SelectedROIListNo=get(hObject,'Value');
roi_members_cell_no = find(handles.roi(handles.roiMaster.Value(1)).members);
% assignin('base','SelectedROIListNo',handles.SelectedROIListNo);
% handles.SelectedROI=str2num(CurrentString{get(hObject,'Value')});
handles.SelectedROI = roi_members_cell_no( get(hObject, 'Value') );

% roiIdx = strcmp({handles.roi.name}, roiname);
% find(handles.roi(roiIdx).members))

% assignin('base','selectedroi',handles.SelectedROI);
% roiIdx = strcmp({handles.roi.name}, roiname);
%-----------Display Corr List------------------
if handles.AlltoAllCheck==1
    [CorrVal handles.CorrCellNo]=sort(handles.Allcor(:,handles.SelectedROI),'descend');
    assignin('base','CorrCellNo',handles.CorrCellNo);
    DistList=DistOfVectorsToVectors(handles.spPos(handles.SelectedROI,:), handles.spPos(handles.CorrCellNo,:));
    StringToDisplay=strcat(num2str(handles.CorrCellNo),' -',num2str(CorrVal),' - ',num2str(DistList));
    set(handles.CorrCellNoList,'String',StringToDisplay);
elseif handles.ROItoROICheck==1
    SelectedROIRowNo=get(hObject,'Value');
	[CorrVal handles.CorrCellNo]=sort(handles.ROIcor(:,SelectedROIRowNo),'descend');
%     assignin('base','CorrCellNo',CorrCellNo);
    DistList=DistOfVectorsToVectors(handles.spPos(handles.SelectedROI,:), handles.spPos(handles.CorrCellNo,:));
%     StringToDisplay=strcat(num2str(C orrCellNo),' -',num2str(CorrVal));

    StringToDisplay=strcat(num2str(handles.CurrentRoiList(handles.CorrCellNo)),'  -',num2str(CorrVal),' - ',num2str(DistList));
    set(handles.CorrCellNoList,'String',StringToDisplay);
    assignin('base','CorrCellNOList',num2str(handles.CurrentRoiList(handles.CorrCellNo)));
end
guidata(hObject,handles);



function [Dist]=DistOfVectorsToVectors(A,B)
% A_{N*d}, B_{M*d},Dist_{N*M}
[N d]=size(A);
[M d]=size(B);
Dist = sqrt(sum(A.*A, 2)*ones(1, M) + ones(N, 1) * sum(B.*B, 2)' - 2*A*B');
Dist=Dist';

% --- Executes during object creation, after setting all properties.
function roiListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function roiMaster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--- Executes on key press with focus on roiListbox and none of its controls.
function roiListbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to roiListbox (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
current_roi = handles.roiMaster.Value;

%note: keep everything inside the if/elif statements to prevent weird
%selection behavior. 
if strcmp(eventdata.Character,'a')
    %     handles.roi = unique( [handles.roi, handles.cellSelect_idx] );
    handles.roi(current_roi).members(handles.cellSelect_idx) = 1;
    handles = display_roiListbox(handles);
    set(handles.ROICellNo,'String',num2str(sum(handles.roi(current_roi).members)));
    guidata(hObject,handles);
elseif strcmp(eventdata.Key,'delete')
    %     handles.roi(handles.roiListbox.Value) = [];
    current_roi = handles.roiMaster.Value;
    roi = find(handles.roi(current_roi).members);
    toDeleteCell = roi(get(handles.roiListbox,'Value') );
%     toDeleteCell = handles.roiListbox.String( handles.roiListbox.Value);
%     toDelete = str2num(toDeleteCell{1});
    handles.roi(current_roi).members(toDeleteCell) = 0;
    handles = display_roiListbox(handles);
    set(handles.ROICellNo,'String',num2str(sum(handles.roi(current_roi).members)));
    guidata(hObject,handles);
    
elseif strcmp(eventdata.Character,'p')||strcmp(eventdata.Character,'q')
    current_roi = handles.roiMaster.Value;
    roi = find(handles.roi(current_roi).members);
    selected_cells = roi(get(handles.roiListbox,'Value') );
    handles=plotfts(handles,selected_cells);
    guidata(hObject,handles);
elseif strcmp(eventdata.Character,'f')    
    current_roi = handles.roiMaster.Value;
    roi = find(handles.roi(current_roi).members);
    selected_cells = roi(get(handles.roiListbox,'Value') );
    handles.CellNoToFind_Special=selected_cells;
    axes(handles.slicePosMap);
%     handles=DrawCirclePoint_Special(hObject,handles);
    
%--------Plot fts, plot slice in slice PosMap, plot SliceAx, draw circle in both-------------
    handles=plotfts(handles, selected_cells);

    handles.currentSlice = ceil( (handles.spPos(selected_cells,1)-min(handles.spPos(:,1)))/handles.sliceThickness);
    handles=plot_slice_maps(handles);
    handles = plot_pos_maps(handles);

    handles.CellNoToFindinSlice=selected_cells;
    axes(handles.sliceAx);
    handles=DrawCirclePointinSlice(handles);

    handles.CellNoToFind=selected_cells;
    axes(handles.slicePosMap);
    handles=DrawCirclePoint(handles);
    
    guidata(hObject,handles);
end


% --- Executes on button press in ExportROI.
function ExportROI_Callback(hObject, eventdata, handles)
% hObject    handle to ExportROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = handles.roi;
[f,path] = uiputfile('roi.mat');
save([path,f],'roi');
%hint: roi struct has fields members (logical), name (string) and ROIcellNo
%(int)

% --- Executes on button press in LoadROI.
function LoadROI_Callback(hObject, eventdata, handles)
% hObject    handle to LoadROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,path] = uigetfile('*.mat');
load([path,f]);
if ~exist('roi','var')
    warning('ROIs must be saved in the variable roi')
end
% handles.roi=[];
handles.roi = [handles.roi roi]; %concatenate ROIs together
assignin('base','roi',handles.roi)
% set(handles.roiListbox,'String',num2cell(find(handles.roi(end).members)));
update_roiMasterList(handles);
handles.roiListbox.Max = sum(handles.roi(end).members);
handles.roiMaster.Value = numel(handles.roi);
handles = display_roiListbox(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function ImarisMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ImarisMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ConnectImaris_Callback(hObject, eventdata, handles)
% hObject    handle to ConnectImaris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.im = GetImaris;
guidata(hObject,handles);

% --------------------------------------------------------------------
function toggleVisibility_im_Callback(hObject, eventdata, handles)
% hObject    handle to toggleVisibility_im (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentVisibility = handles.im.GetVisible;
handles.im.SetVisible(~currentVisibility);

% --------------------------------------------------------------------
function exportSpots_Callback(hObject, eventdata, handles)
% hObject    handle to exportSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pxyz = handles.spPos(find(handles.roi.members),:);
pxyz(:,4) = 0;
sH = MakeImarisSpots(pxyz,[1,0,0,0],'ROI',handles.im);
celldiameters = repmat([5,5,7],size(pxyz,1),1);
sH.SetRadiiXYZ(celldiameters./2);


% --------------------------------------------------------------------
function im_importSpots_Callback(hObject, eventdata, handles)
% hObject    handle to im_importSpots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
objNames = CheckObjects(handles.im,'Spots');
[indx,tf] = listdlg('ListString', objNames(:,1),'SelectionMode','single','PromptString','Select Imaris Spots to Import to ROI');
objSelected = objNames{indx,2};
objPos = objSelected.GetPositionsXYZ;
distance_matrix =  pdist2(handles.spPos, objPos);
exactMatches = distance_matrix == 0;
labelCells = logical(sum(exactMatches,2));
if sum(labelCells) == 0
    error('Could not find exact position matches!!')
end

roi.members = labelCells;
roi.name = objNames{indx,1};
% handles.roi = labelCells;
handles = update_roi(handles, roi,'add');
% handles = display_roiListbox(handles, roi.name);

guidata(hObject,handles);

function handles = display_roiListbox(handles)
%displays current members of roi listed by roiname
% disp( roiname );
% roiIdx = strcmp({handles.roi.name}, roiname);
roiIdx =get(handles.roiMaster,'Value');
% assignin('base','currentmembers',handles.roi(roiIdx).members);
set(handles.roiListbox,'String',num2cell(find(handles.roi(roiIdx).members)))
handles.roiListbox.Max = sum(handles.roi(roiIdx).members);
% handles.roiListbox.Value = 1;
if handles.roiListbox.Value > handles.roiListbox.Max
    handles.roiListbox.Value = handles.roiListbox.Max;
end

function update_roiMasterList(handles)
%updates the roiMaster listbox with current information
set(handles.roiMaster,'String', {handles.roi.name});
handles.roiMaster.Max = numel(handles.roi);
if isempty(handles.roiMaster.Value) || handles.roiMaster.Value > handles.roiMaster.Max
    handles.roiMaster.Value = 1;
end

function handles = update_roi(handles, roi, command)
%function for manipulating the roi textboxes
%handles.roi is a struct with #elements equal to number of roi groups with
%fields members (logical #cells), name (string), and others?
roistruct = handles.roi;

if strcmp(command,'add')
    roistruct(end+1) = roi;
end
handles.roi = roistruct;
update_roiMasterList(handles);

% --- Executes during object creation, after setting all properties.
function set_slice_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_slice_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in roiMaster.
function roiMaster_Callback(hObject, eventdata, handles)
% hObject    handle to roiMaster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiMaster contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiMaster

% set(handles.CorrCellNoList,'Max',9999999999999,'Min',0);
handles = guidata(hObject);

handles = display_roiListbox(handles);
handles = plot_pos_maps(handles);
handles = plot_slice_maps(handles);
handles = plot_fts_in_slice(handles);

%WTF IS THIS? TALK TO DAWNIS BEFORE IMPLEMENTING THIS KIND OF MYSTERIOUS
%CRAP. 09/25/2018
% axes(handles.dffPlot), cla %this bit of code will generate a warning if the cell has changed but trial has not. not important. 
% % yyaxis left
% hold on
% 
% CellNoToPlot=find(handles.roi(get(hObject,'Value')).members);
% if get(handles.raw_trace_selector,'Value')
%     ftsToPlot=mean(handles.fts(CellNoToPlot,:),1);
%     plotTraces_GUI(ftsToPlot, 100, 4, handles.dffPlot);
%     colorbar('off');
%     ylabel('Intensity - Mean');
% elseif get(handles.dFF_traces,'Value')
%     dFFToPlot=mean(handles.dFF(CellNoToPlot,:),1);
%     plotTraces_GUI(dFFToPlot, 1, 3, handles.dffPlot);
%     colorbar('off');
%     ylabel('\Delta F/F');
% elseif get(handles.heatmap_selector,'Value')
%     toplot=mean(handles.fts(CellNoToPlot,:),1);
% %     toplot = handles.fts(cellSelect,:);
%     toplot = bsxfun(@minus, toplot, median(toplot,2));
%     imagesc(handles.dffPlot, toplot );
%     axis tight
%     CB1=colorbar;  %colorbar('on');
%     CB1.Color='w';
%     colormap jet;
% elseif get(handles.dFF_heatmap,'Value')
%     toplot=mean(handles.dFF(CellNoToPlot,:),1);
% %     toplot = handles.dFF(cellSelect,:);
%     imagesc(handles.dffPlot, toplot );
%     axis tight
%     CB1=colorbar;  %colorbar('on');
%     CB1.Color='w';
%     colormap jet;
% end
% 
% xlabel('Frames');
% title('Fluorescence Time Series', 'color', [1 1 1]);
% set(handles.dffPlot,'FontUnits','normalized','FontSize',0.06);
% handles.dffPlot = gca;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uipanel12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in roi_plot_fts_button.
function roi_plot_fts_button_Callback(hObject, eventdata, handles)
% hObject    handle to roi_plot_fts_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_roi = handles.roiMaster.Value;

roi = find(handles.roi(current_roi).members);
selected_cells = roi(get(handles.roiListbox,'Value') );
handles=plotfts(handles,selected_cells);
guidata(hObject,handles);


% --- Executes on button press in dffButton.
function dffButton_Callback(hObject, eventdata, handles)
% hObject    handle to dffButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dFF = GUI_dFF_Awesome(handles.fts);
guidata(hObject,handles);


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Corr_AlltoAll_Callback(hObject, eventdata, handles)
% hObject    handle to Corr_AlltoAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer1 = inputdlg('What data do you want to calculate correlation with? Type dff for dF/F, fluo for fluorescence.','dFF or Fluo');
    if strcmp(answer1{1,1},'dff')||strcmp(answer1{1,1},'dFF')
        if isnan(handles.dFF)
            warndlg('dFF of this dataset is not calculated yet, please calculate dFF first! Or use Fluo for correlation instead.','No dFF calculated')
        else
            print1='Calculating...';
            disp(print1);
            if strcmp(handles.CurrentCorrType,'Pearson')
                [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.dFF);
            elseif strcmp(handles.CurrentCorrType,'Spearman')
                [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Spearman(handles.dFF);                
            else
                warndlg('Wrong function input! Please select Spearmn or Pearson');
            end
        end
        
    elseif strcmp(answer1{1,1},'fluo')||strcmp(answer1{1,1},'Fluo')
        print1='Calculating...';
        disp(print1);
        if strcmp(handles.CurrentCorrType,'Pearson')
            [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.fts);            
        elseif strcmp(handles.CurrentCorrType,'Spearman')
            [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.fts);            
        else
            warndlg('Wrong function input! Please select Spearmn or Pearson');
        end
        
    else
        warndlg('Wrong data type input! Please type dff or fluo, case unsensitive');
end

    
% assignin('base','AllMaxcorno',handles.AllMaxcorno);
% assignin('base','AllMaxcorval',handles.AllMaxcorval);
% assignin('base','Allcor',handles.Allcor);
% msgbox({'All Correlation Calculated!' 'Stored in handles.AllMaxcorval and handles.AllMaxcorno' 'You can also check AllMaxcorval & AllMaxcorno in the base workspace'});
msgbox({'All Correlation Calculated!'});
guidata(hObject,handles);

% --------------------------------------------------------------------
function Corr_ROItoROI_Callback(hObject, eventdata, handles)
% hObject    handle to Corr_ROItoROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer1 = inputdlg('What data do you want to calculate correlation with? Type dff for dF/F, fluo for fluorescence.','dFF or Fluo');

CurrentROIString=get(handles.roiListbox,'String');
% assignin('base','CurrentROIString',CurrentROIString);
handles.CurrentRoiList=find(handles.roi.members);
% assignin('base','CurrentRoiList',handles.CurrentRoiList);


if strcmp(answer1{1,1},'dff')||strcmp(answer1{1,1},'dFF')
    if isnan(handles.dFF)
        warndlg('dFF of this dataset is not calculated yet, please calculate dFF first! Or use Fluo for correlation instead.','No dFF calculated')
    else
        print1='Calculating...';
        disp(print1);
        if strcmp(handles.CurrentCorrType,'Pearson')
            [handles.ROIcor, handles.ROIPVal]=TopCorrCellInGroup_Pearson(handles.dFF(handles.CurrentRoiList,:));
        elseif strcmp(handles.CurrentCorrType,'Spearman')
            [handles.ROIcor, handles.ROIPVal]=TopCorrCellInGroup_Spearman(handles.dFF(handles.CurrentRoiList,:));
        else
            warndlg('Wrong function input! Please select Spearmn or Pearson');
        end
    end

elseif strcmp(answer1{1,1},'fluo')||strcmp(answer1{1,1},'Fluo')
    print1='Calculating...';
    disp(print1);
    if strcmp(handles.CurrentCorrType,'Pearson')
        [handles.ROIcor, handles.ROIPVal]=TopCorrCellInGroup_Pearson(handles.fts(handles.CurrentRoiList,:));
    elseif strcmp(handles.CurrentCorrType,'Spearman')
        [handles.ROIcor, handles.ROIPVal]=TopCorrCellInGroup_Spearman(handles.fts(handles.CurrentRoiList,:));
    else
        warndlg('Wrong function input! Please select Spearmn or Pearson');
    end

else
    warndlg('Wrong data type input! Please type dff or fluo, case unsensitive');
end

% assignin('base','ROIMaxcorno',handles.ROIMaxcorno);
% assignin('base','ROIMaxcorval',handles.ROIMaxcorval);
% assignin('base','ROI',handles.CurrentRoiList);
% assignin('base','ROIcor',handles.ROIcor);
msgbox({'All Correlation Calculated!'});
guidata(hObject,handles);


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Corr_SeriesToAll_Callback(hObject, eventdata, handles)
% hObject    handle to Corr_SeriesToAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%-------------Load Target Files-------------
[SeriesFile,SeriesFilePath] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(SeriesFilePath);
SeriesFileDirectory= [SeriesFilePath,SeriesFile];
SeriesFile=load(SeriesFileDirectory);

answer1 = inputdlg('What data do you want to calculate correlation with? Type dff for dF/F, fluo for fluorescence.','dFF or Fluo');
answer2 = inputdlg('Do you want to only consider certain range of frames for correlation? If yes please type 1, otherwise type whatever else.','Time Range?');
if strcmp(answer2{1,1},'1')
    answer3=inputdlg({'Lower Time Range','Upper Time Range'});
    LowerTimeRange=str2num(answer3{1,1});
    UpperTimeRange=str2num(answer3{2,1});
else
    LowerTimeRange=1;
    UpperTimeRange=size(handles.fts,2);
end

SeriesFile_Cell=struct2cell(SeriesFile);
TargetSeries=SeriesFile_Cell{1,1};
handles.TargetSeries=TargetSeries;

answer4= inputdlg('Threshold of dFF or Fluo? If no threshold, type 0');
threshold=str2num(answer4{1,1});

%----------Save some variables for export------------
handles.CurrentTimeRange=[LowerTimeRange,UpperTimeRange];
handles.DataTypeUsedForCorr=answer1{1,1};
guidata(hObject,handles);

%------------Calculating part----------------
k=size(handles.fts,1);%Here I set the program to display the correlation of all cells to input series. Change here if you only want to see top K ones.

TargetSeriesNo=size(TargetSeries,1);%Change here if your target series are more than one.

% plot(1:length(SeriesFile.StimTimeSeries),SeriesFile.StimTimeSeries,'Parent',handles.StimTSPlot);
plot(1:length(TargetSeries),TargetSeries,'Parent',handles.StimTSPlot);
xLimits = get(handles.dffPlot,'XLim');
set(handles.StimTSPlot,'XLim',xLimits);
set(handles.StimTSPlot,'XColor',[1 1 1]);
        if strcmp(answer1{1,1},'dff')||strcmp(answer1{1,1},'dFF')
            if isnan(handles.dFF)
                warndlg('dFF of this dataset is not calculated yet, please calculate dFF first! Or use Fluo for correlation instead.','No dFF calculated')
            else
                print1='Calculating...';
                disp(print1);
                if strcmp(handles.CurrentCorrType,'Pearson')
                    assignin('base','inputfts',handles.dFF(:,LowerTimeRange:UpperTimeRange));
                    assignin('base','inputseries',TargetSeries(:,LowerTimeRange:UpperTimeRange));
                    [handles.CorrValToSeries, handles.CorrCellNoToSeries, handles.AllCorrToSeries, handles.AllPValtoSeries,LowdFFList]=TopkCorr_Pearson(handles.dFF(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k,threshold);                    
                elseif strcmp(handles.CurrentCorrType,'Spearman')
                    [handles.CorrValToSeries, handles.CorrCellNoToSeries, handles.AllCorrToSeries, handles.AllPValtoSeries,LowdFFList]=TopkCorr_Spearman(handles.dFF(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k,threshold);
                else
                    warndlg('Wrong function input! Please select Spearmn or Pearson');
                end
            end
%             handles.CorrValToSeries=handles.CorrValToSeries';handles.CorrCellNoToSeries=handles.CorrCellNoToSeries'; 
%             handles.AllCorrToSeries=handles.AllCorrToSeries';handles.AllPValtoSeries=handles.AllPValtoSeries';

            StringToDisplay=strcat(num2str(handles.CorrCellNoToSeries),' -',num2str(handles.CorrValToSeries),' -',num2str(handles.AllPValtoSeries));
            set(handles.CorrCellNoList,'String',StringToDisplay);
        elseif strcmp(answer1{1,1},'fluo')||strcmp(answer1{1,1},'Fluo')
            print1='Calculating...';
            disp(print1);
            if strcmp(handles.CurrentCorrType,'Pearson')
                [handles.CorrValToSeries, handles.CorrCellNoToSeries, handles.AllCorrToSeries, handles.AllPValtoSeries,LowdFFList]=TopkCorr_Pearson(handles.fts(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k,threshold);
            elseif strcmp(handles.CurrentCorrType,'Spearman');
                [handles.CorrValToSeries, handles.CorrCellNoToSeries, handles.AllCorrToSeries, handles.AllPValtoSeries,LowdFFList]=TopkCorr_Spearman(handles.fts(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k,threshold);
            else
                warndlg('Wrong function input! Please select Spearmn or Pearson');
            end
%             handles.CorrValToSeries=handles.CorrValToSeries';handles.CorrCellNoToSeries=handles.CorrCellNoToSeries';
%             handles.AllCorrToSeries=handles.AllCorrToSeries';handles.AllPValtoSeries=handles.AllPValtoSeries';

            StringToDisplay=strcat(num2str(handles.CorrCellNoToSeries),' -',num2str(handles.CorrValToSeries),' -',num2str(handles.AllPValtoSeries));
            set(handles.CorrCellNoList,'String',StringToDisplay);

        else
            warndlg('Wrong data type input! Please type dff or fluo, case unsensitive');
        end
    assignin('base','FilteredOutCellList',LowdFFList);
    msgbox({strcat('Calculation Done! ',num2str(size(LowdFFList,1)),' cells were excluded for correlation rank (for detail check FilteredOutCellList in base workspace). Result displayed in Correlation Listbox')});
%----------Save some variables for export.2------------
    handles.ThresholdUsedForCorrRank=threshold;
    handles.CellList_BelowThreshold=LowdFFList;
 guidata(hObject,handles);



% --- Executes on selection change in CorrCellNoList.
function CorrCellNoList_Callback(hObject, eventdata, handles)
% hObject    handle to CorrCellNoList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrCellNoList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrCellNoList

% assignin('base','CorrCellNoList',get(handles.CorrCellNoList,'String'));
% handles.SelectedCorCellListNo=get(hObject,'Value');
% CurrentString=str2num(get(handles.CorrCellNoList,'String'));

set(handles.CorrCellNoList,'Max',99999999999,'Min',0);
set(handles.CurrentCorrCellSeqNo,'String',num2str(get(hObject,'Value')));

if handles.SeriestoAllCheck==1
    handles.SelectedCorrCellNo=handles.CorrCellNoToSeries(get(hObject,'Value'));
    selectedcellNo=length(handles.SelectedCorrCellNo)
else
    handles.SelectedCorrCellNo=handles.CorrCellNo(get(hObject,'Value'));
end

guidata(hObject,handles);
% plot(1:size(ndles.fts,2),handles.fts(handles.SelectedCorrCellNo),'Parent',handles.dffPlot);

assignin('base','SelectedCorrCellNo',handles.SelectedCorrCellNo);
% if handles.SeriesToAllCalculated==1

%--------Plot fts, plot slice in slice PosMap, plot SliceAx, draw circle in both-------------
handles=plotfts(handles, handles.SelectedCorrCellNo);

handles.currentSlice = ceil( (handles.spPos(handles.SelectedCorrCellNo,1)-min(handles.spPos(:,1)))/handles.sliceThickness);%ceil( handles.spPos(handles.SelectedCorrCellNo,1)/handles.sliceThickness);
handles=plot_slice_maps(handles);
handles = plot_pos_maps(handles);

handles.CellNoToFindinSlice=handles.SelectedCorrCellNo;
axes(handles.sliceAx);
handles=DrawCirclePointinSlice(handles);

handles.CellNoToFind=handles.SelectedCorrCellNo;
axes(handles.slicePosMap);
handles=DrawCirclePoint(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function CorrCellNoList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrCellNoList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in PlotCellPos.
function PlotCellPos_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCellPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ExportROIMean.
function ExportROIMean_Callback(hObject, eventdata, handles)
% hObject    handle to ExportROIMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_roi = handles.roiMaster.Value;
roi = find(handles.roi(current_roi).members);
selected_cells = roi(get(handles.roiListbox,'Value') );


% CurrentRoi=find(handles.roi.members);
% assignin('base','CurrentRoiList',CurrentRoi);
ROIMeanFluo=mean(handles.fts(selected_cells,:),1);
if isfield(handles,'dFF')
    ROIMeandFF=mean(handles.dFF(selected_cells,:),1);
end
ROIList=selected_cells;
[f2,path2] = uiputfile('ROIMeanFluo&ROIList.mat');
save([path2,f2],'ROIMeanFluo','ROIMeandFF','ROIList');


Fig1=figure;
plot((1:length(ROIMeanFluo)),ROIMeanFluo);
title('Mean Fluo of Selected Cells');
savefig(Fig1,strcat(path2,f2,'-Fluo.fig'));

if isfield(handles,'dFF')
    Fig2=figure;
    plot((1:length(ROIMeandFF)),ROIMeandFF);
    title('Mean dFF of Selected Cells');
    savefig(Fig2,strcat(path2,f2,'-dFF.fig'));
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
if get(hObject,'Value')==1
    handles.SeriestoAllCheck=1;
    handles.AlltoAllCheck=0;
    handles.ROItoROICheck=0;
elseif get(hObject,'Value')==2
    handles.SeriestoAllCheck=0;
    handles.AlltoAllCheck=1;
    handles.ROItoROICheck=0;
elseif get(hObject,'Value')==3
    handles.SeriestoAllCheck=0;
    handles.AlltoAllCheck=0;
    handles.ROItoROICheck=1;
end
guidata(hObject,handles);
    

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    handles.SeriestoAllCheck=1;
    handles.AlltoAllCheck=0;
    handles.ROItoROICheck=0;
guidata(hObject,handles);



function CellNoToFind_Callback(hObject, eventdata, handles)
% hObject    handle to CellNoToFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNoToFind as text
%        str2double(get(hObject,'String')) returns contents of CellNoToFind as a double
handles.CellNoToFind_Special=str2double(get(hObject,'String'));
%--------Plot fts, plot slice in slice PosMap, plot SliceAx, draw circle in both-------------
handles=plotfts(handles, handles.CellNoToFind_Special);

handles.currentSlice = ceil((handles.spPos(handles.CellNoToFind_Special,1)-min(handles.spPos(:,1)))/handles.sliceThickness);
handles = plot_slice_maps(handles);
handles = plot_pos_maps(handles);

handles.CellNoToFindinSlice=handles.CellNoToFind_Special;
axes(handles.sliceAx);
handles=DrawCirclePointinSlice(handles);

handles.CellNoToFind=handles.CellNoToFind_Special;
axes(handles.slicePosMap);
handles=DrawCirclePoint(handles);

set(handles.XVal, 'String', num2str(handles.spPos(handles.CellNoToFind_Special,1)));
set(handles.YVal, 'String', num2str(handles.spPos(handles.CellNoToFind_Special,2)));
set(handles.ZVal, 'String', num2str(handles.spPos(handles.CellNoToFind_Special,3)));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function CellNoToFind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellNoToFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FTSViewer.
function FTSViewer_Callback(hObject, eventdata, handles)
% hObject    handle to FTSViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=CXYsDataExplorer_Copied2(handles);
guidata(hObject,handles);

% --- Executes on button press in CompressCellNo.
function CompressCellNo_Callback(hObject, eventdata, handles)
% hObject    handle to CompressCellNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CorrMapGUI.
function CorrMapGUI_Callback(hObject, eventdata, handles)
% hObject    handle to CorrMapGUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TopCorrViewer_v2;


function CellNoToAdd_Callback(hObject, eventdata, handles)
% hObject    handle to CellNoToAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNoToAdd as text
%        str2double(get(hObject,'String')) returns contents of CellNoToAdd as a double

current_roi = handles.roiMaster.Value;
handles.CellNoToAdd=str2double(get(hObject,'String'));
handles.roi(current_roi).members(handles.CellNoToAdd) = 1;
handles = display_roiListbox(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function CellNoToAdd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellNoToAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CropDataset.
function CropDataset_Callback(hObject, eventdata, handles)
% hObject    handle to CropDataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.CroppedDataPointList = GUI_DataCropParameterInput(handles.spPos,handles.fts,handles.CurrentRoiList);

DataCropper;
guidata(hObject,handles);



% --- Executes on selection change in PlotSelect.
function PlotSelect_Callback(hObject, eventdata, handles)
% hObject    handle to PlotSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotSelect
handles=guidata(hObject);
%rescale sliceSelector and cellSelector to accomdate differing functions. 
handles = reset_slider_handles(handles);

if get(hObject, 'Value') < 3
    axes(handles.sliceAx); view(3);
    set(handles.load_z_stack, 'Visible', 'off');
    set(handles.throw_roi_button, 'Visible', 'off');
    set(handles.slicePosMap, 'Visible', 'on');
    colormap(handles.slicePosMap, 'jet')
else
    axes(handles.sliceAx); view(2); cla;
    axes(handles.slicePosMap); cla;
    set(handles.load_z_stack, 'Visible', 'on');
    set(handles.throw_roi_button, 'Visible', 'on');
    handles.caxis0.String = '';
    handles.caxis1.String = '';
    handles.sPMap_Ax_roi = [];
    handles.cell_roi_list = [];
    handles.posmap_img = [];
    
end

handles = plot_pos_maps(handles);
plot_slice_maps(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PlotSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dFF_traces.
function dFF_traces_Callback(hObject, eventdata, handles)
% hObject    handle to dFF_traces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dFF_traces
% handles=plotfts(handles,handles.CellSelectedinSlice);
% handles = plot_fts_in_slice(handles);
handles=plotfts(handles,handles.CellNoToFindinSlice);

% --- Executes on button press in raw_trace_selector.
function raw_trace_selector_Callback(hObject, eventdata, handles)
% hObject    handle to raw_trace_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw_trace_selector
% handles=plotfts(handles,handles.CellSelectedinSlice);
% handles = plot_fts_in_slice(handles);
handles=plotfts(handles,handles.CellNoToFindinSlice);


% --- Executes on button press in heatmap_selector.
function heatmap_selector_Callback(hObject, eventdata, handles)
% hObject    handle to heatmap_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of heatmap_selector
% handles=plotfts(handles,handles.CellSelectedinSlice);
% handles = plot_fts_in_slice(handles);
handles=plotfts(handles,handles.CellNoToFindinSlice);

% --- Executes on button press in dFF_heatmap.
function dFF_heatmap_Callback(hObject, eventdata, handles)
% hObject    handle to dFF_heatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dFF_heatmap

% handles=plotfts(handles,handles.CellSelectedinSlice);
% handles = plot_fts_in_slice(handles);
handles=plotfts(handles,handles.CellNoToFindinSlice);

% --- Executes on button press in CheckROIInRawData.
function CheckROIInRawData_Callback(hObject, eventdata, handles)
% hObject    handle to CheckROIInRawData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

current_roi = handles.roiMaster.Value;
roi = find(handles.roi(current_roi).members);
assignin('base','roi',roi);
genFIJI_ROIs_CuiEdited(handles.spPos(roi,:), roi, handles.spRadiiXYZ(roi,:));


% --- Executes on selection change in CorrTypeSelect.
function CorrTypeSelect_Callback(hObject, eventdata, handles)
% hObject    handle to CorrTypeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrTypeSelect contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrTypeSelect
if get(hObject,'Value') ==1
    handles.CurrentCorrType='Spearman';
    print1='Corr type set to Spearman';
    disp(print1);
elseif get(hObject,'Value')==2
    handles.CurrentCorrType='Pearson';
    print1='Corr type set to Pearson';
    disp(print1);
end
guidata(hObject,handles);
    

% --- Executes during object creation, after setting all properties.
function CorrTypeSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrTypeSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.CurrentCorrType='Spearman';
guidata(hObject,handles);


% --- Executes on button press in CellFinderStay.
function CellFinderStay_Callback(hObject, eventdata, handles)
% hObject    handle to CellFinderStay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CellFinderStay


% --- Executes on key press with focus on CorrCellNoList and none of its controls.
function CorrCellNoList_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to CorrCellNoList (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Character,'a')
    %     handles.roi = unique( [handles.roi, handles.cellSelect_idx] );
    current_roi = handles.roiMaster.Value;
    handles.roi(current_roi).members(handles.SelectedCorrCellNo) = 1;
    AddCellNo=length(handles.SelectedCorrCellNo);
    handles = display_roiListbox(handles);
    guidata(hObject,handles);    
    set(handles.ROICellNo,'String',num2str(sum(handles.roi(current_roi).members)));
end


% --------------------------------------------------------------------
function CrossCorr_SeriesToAll_Callback(hObject, eventdata, handles)
% hObject    handle to CrossCorr_SeriesToAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[SeriesFile,SeriesFilePath] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(SeriesFilePath);
SeriesFileDirectory= [SeriesFilePath,SeriesFile];
SeriesFile=load(SeriesFileDirectory);

answer1 = inputdlg('What data do you want to calculate correlation with? Type dff for dF/F, fluo for fluorescence.','dFF or Fluo');
answer2 = inputdlg('Do you want to only consider certain range of frames for correlation? If yes please type 1','Time Range?');
if strcmp(answer2{1,1},'1')
    answer3=inputdlg({'Lower Time Range','Upper Time Range'});
    LowerTimeRange=str2num(answer3{1,1});
    UpperTimeRange=str2num(answer3{2,1});
else
    LowerTimeRange=1;
    UpperTimeRange=size(handles.fts,2);
end

answer4=inputdlg('LagRange?');
LagRange=str2num(answer4{1,1});

k=size(handles.fts,1);%Here I set the program to display the correlation of all cells to input series. Change here if you only want to see top K ones.

% TargetSeries=SeriesFile.ROIMeanFluo;
TargetSeries=SeriesFile.StimTimeSeries; %Change this back to the line above if you want to use ROI mean as target fts!

TargetSeriesNo=size(TargetSeries,1);%Change here if your target series are more than one.

plot(1:length(SeriesFile.StimTimeSeries),SeriesFile.StimTimeSeries,'Parent',handles.StimTSPlot);
xLimits = get(handles.dffPlot,'XLim');
set(handles.StimTSPlot,'XLim',xLimits);
set(handles.StimTSPlot,'XColor',[1 1 1]);
        if strcmp(answer1{1,1},'dff')||strcmp(answer1{1,1},'dFF')
            if isnan(handles.dFF)
                warndlg('dFF of this dataset is not calculated yet, please calculate dFF first! Or use Fluo for correlation instead.','No dFF calculated')
            else
                print1='Calculating...';
                disp(print1);

                [handles.CrsCorrValToSeries, handles.CrsCorrCellNoToSeries, handles.AllCrsCorrToSeries, handles.AllLagTimeToSeries]=TopkCrossCorr(handles.dFF(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k, LagRange);
                
            end
%             handles.CorrValToSeries=handles.CorrValToSeries';handles.CorrCellNoToSeries=handles.CorrCellNoToSeries'; 
%             handles.AllCorrToSeries=handles.AllCorrToSeries';handles.AllPValtoSeries=handles.AllPValtoSeries';

            StringToDisplay=strcat(num2str(handles.CorrCellNoToSeries),' -',num2str(handles.CorrValToSeries),' -',num2str(handles.AllPValtoSeries));
            set(handles.CorrCellNoList,'String',StringToDisplay);
        elseif strcmp(answer1{1,1},'fluo')||strcmp(answer1{1,1},'Fluo')
            print1='Calculating...';
            disp(print1);
                
            [handles.CrsCorrValToSeries, handles.CrsCorrCellNoToSeries, handles.AllCrsCorrToSeries, handles.AllLagTimeToSeries]=TopkCrossCorr(handles.fts(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k, LagRange);
                
%             handles.CorrValToSeries=handles.CorrValToSeries';handles.CorrCellNoToSeries=handles.CorrCellNoToSeries';
%             handles.AllCorrToSeries=handles.AllCorrToSeries';handles.AllPValtoSeries=handles.AllPValtoSeries';

            StringToDisplay=strcat(num2str(handles.CorrCellNoToSeries),' -',num2str(handles.CorrValToSeries),' -',num2str(handles.AllPValtoSeries));
            set(handles.CorrCellNoList,'String',StringToDisplay);

        else
            warndlg('Wrong data type input! Please type dff or fluo, case unsensitive');
    end
    msgbox({'Calculation Done! Result displayed in Correlation Listbox'});
    

guidata(hObject,handles);



% --- Executes on button press in ExportCorr.
function ExportCorr_Callback(hObject, eventdata, handles)
% hObject    handle to ExportCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CorrValToSeries=handles.CorrValToSeries;
CorrCellNoToSeries=handles.CorrCellNoToSeries;
AllCorrToSeries=handles.AllCorrToSeries;
AllPValtoSeries=handles.AllPValtoSeries;
CurrentCorrType=handles.CurrentCorrType;
CurrentTimeRange=handles.CurrentTimeRange;
DataTypeUsedForCorr=handles.DataTypeUsedForCorr;

ThresholdUsedForCorrRank=handles.ThresholdUsedForCorrRank;
CellList_BelowThreshold=handles.CellList_BelowThreshold;


[f,path] = uiputfile('CorrResult.mat');
save([path,f],'CorrValToSeries','CorrCellNoToSeries','AllCorrToSeries','AllPValtoSeries','CurrentCorrType','CurrentTimeRange','DataTypeUsedForCorr','ThresholdUsedForCorrRank','CellList_BelowThreshold');
guidata(hObject,handles);


% --- Executes on button press in LoadCorr.
function LoadCorr_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[f,path] = uigetfile('*.mat');
load([path,f]);

handles.CorrValToSeries=CorrValToSeries;
handles.CorrCellNoToSeries=CorrCellNoToSeries;
handles.AllCorrToSeries=AllCorrToSeries;
handles.AllPValtoSeries=AllPValtoSeries;

handles.CurrentCorrType=CurrentCorrType;
if strcmp(handles.CurrentCorrType,'Spearman')==1
%     set(handles.CorrTypeSelect,'Value')=1;
    warndlg('Corr type should be set to Spearman!');
elseif strcmp(handles.CurrentCorrType,'Pearson')==1
%     set(handles.CorrTypeSelect,'Value')=2;
    warndlg('Corr type should be set to Pearson!');
end

handles.CurrentTimeRange=CurrentTimeRange;
handles.DataTypeUsedForCorr=DataTypeUsedForCorr;
handles.ThresholdUsedForCorrRank=ThresholdUsedForCorrRank;
handles.CellList_BelowThreshold=CellList_BelowThreshold;

StringToDisplay=strcat(num2str(handles.CorrCellNoToSeries),' -',num2str(handles.CorrValToSeries),' -',num2str(handles.AllPValtoSeries));
set(handles.CorrCellNoList,'String',StringToDisplay);

guidata(hObject,handles);


% --- Executes on button press in ExportFig.
function ExportFig_Callback(hObject, eventdata, handles)
% hObject    handle to ExportFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(gca, handles.slicePosMap);

Fig2 = figure;
newAxes=copyobj(handles.slicePosMap, Fig2);
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(Fig2,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
    xlabel('X','Color','k');
    ylabel('Y','Color','k');
    zlabel('Z','Color','k');
    colormap jet;
    grid on
    axis equal
%     CB1=colorbar
%     caxis([-0.5 1]);
%     CB1.Label.String='Correlation Coef of Cells (Dots)';
%     tit2= strcat('Total Correlation Map of',answer3{1,1});
    view(0,90);
%     title(tit2);
%     savefig(Fig3,strcat(path,'\',tit2));
[f,path] = uiputfile('Figure.fig');
savefig(Fig2,[path,f]);


% --- Executes on button press in CorrMovie.
function CorrMovie_Callback(hObject, eventdata, handles)
% hObject    handle to CorrMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=inputdlg({'Printing Range is from how many';'to how many in the list?';'spacing interval?'});
PrintingRange=[str2num(answer{1,1}),str2num(answer{2,1})];
SpacingInt=str2num(answer{3,1});

path= uigetdir();

answer2=inputdlg('Is your target series a cell fluo fts? 1 for yes 0 for no');
if strcmp(answer2{1,1},'1')
    answer3=inputdlg('What is its cell No?');
    InputCellNo=str2num(answer3{1,1});
    answer4=inputdlg('Is it a correlation based on dFF or Fluo? Type dff or fluo');
    Fig2=figure;
    if strcmp(answer4{1,1},'dff')||strcmp(answer4{1,1},'dFF')
        TargetSeries=handles.dFF(str2num(answer3{1,1}),:);
        plot(1:size(handles.dFF,2),TargetSeries);
        title(strcat('dFF of cell ',answer3{1,1}));
        f2=strcat('dFF of cell ',answer3{1,1});
    elseif strcmp(answer4{1,1},'fluo')||strcmp(answer4{1,1},'Fluo')
        TargetSeries=handles.fts(str2num(answer3{1,1}),:);
        plot(1:size(handles.fts,2),TargetSeries);
        title(strcat('Fluo of cell ',answer3{1,1}));
        f2=strcat('Fluo of cell ',answer3{1,1});        
    end
    savefig(Fig2,strcat(path,'\',f2));
else
    answer3=inputdlg('Input Target Series Name');
    Fig2=figure;
    plot(1:length(handles.TargetSeries),handles.TargetSeries);
    f2=strcat('Stim Time Series of ',answer3{1,1});        
    title(f2);
    savefig(Fig2,strcat(path,'\',f2,'.fig'));
end

answer5=inputdlg('Plot first and last top k correlated cells concurrently? 1 for yes; 0 for no; 3 for both');

loopNo=abs(ceil((PrintingRange(2)-PrintingRange(1))/SpacingInt));


    Fig3 = figure;
    ColorBasedOnCorr=handles.CorrValToSeries;
    NewPos=handles.spPos(handles.CorrCellNoToSeries,:);
    scatter3(NewPos(:,1),NewPos(:,2),NewPos(:,3),20,ColorBasedOnCorr,'.','hittest','off');
    hold on
    if strcmp(answer2{1,1},'1')
        plot3(handles.spPos(InputCellNo,1),handles.spPos(InputCellNo,2),handles.spPos(InputCellNo,3),'ko','linewidth',3,'Markersize',10);
    end
    colormap jet;
    grid on
    axis equal
    xlabel('X','Color','k');
    ylabel('Y','Color','k');
    zlabel('Z','Color','k');
    CB1=colorbar;
%     caxis([-0.5 1]);
    CB1.Label.String='Correlation Coef of Cells (Dots)';
    tit2= strcat('Total Correlation Map of',answer3{1,1});
    view(0,90);
    title(tit2);
    savefig(Fig3,strcat(path,'\',tit2,'.fig'));
    
    
for i=1:loopNo
    if PrintingRange(2)-PrintingRange(1)>0
        selectedCell=handles.CorrCellNoToSeries(PrintingRange(1):(PrintingRange(1)+i*SpacingInt));
        CurrentCorrVal=handles.CorrValToSeries(PrintingRange(1)+i*SpacingInt);
        if strcmp(answer5{1,1},'1')||strcmp(answer5{1,1},'3')
            AllCellNo=size(handles.CorrCellNoToSeries,1);
            anti_selectedCell=handles.CorrCellNoToSeries((AllCellNo-i*SpacingInt):AllCellNo);
            anti_CurrentCorrVal=handles.CorrValToSeries(AllCellNo-i*SpacingInt);
        end
    else
        selectedCell=handles.CorrCellNoToSeries((PrintingRange(1)-i*SpacingInt):PrintingRange(1));
%         assignin('base','selectedCell',selectedCell);
        CurrentCorrVal=handles.CorrValToSeries(PrintingRange(1)-i*SpacingInt);
        if strcmp(answer5{1,1},'1')||strcmp(answer5{1,1},'3')
            anti_selectedCell=handles.CorrCellNoToSeries(1:(1+i*SpacingInt));
            anti_CurrentCorrVal=handles.CorrValToSeries(1+i*SpacingInt);
        end
    end

    Fig3 = figure;
    posB = handles.spPos(handles.labelVector,:);
    posR = handles.spPos(selectedCell,:);
    plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
    hold on;
    plot3(posR(:,1),posR(:,2),posR(:,3),'ro', 'markerfacecolor','r','markersize',3);
    if strcmp(answer5{1,1},'1')||strcmp(answer5{1,1},'3')
        posG = handles.spPos(anti_selectedCell,:);
        plot3(posG(:,1),posG(:,2),posG(:,3),'go', 'markerfacecolor','g','markersize',3);
    end
    if strcmp(answer2{1,1},'1')
        plot3(handles.spPos(InputCellNo,1),handles.spPos(InputCellNo,2),handles.spPos(InputCellNo,3),'bo','linewidth',3,'Markersize',10);
    end
    grid on
    axis equal
    xlabel('X','Color','k');
    ylabel('Y','Color','k');
    zlabel('Z','Color','k');
    
    if strcmp(answer5{1,1},'1')||strcmp(answer5{1,1},'3')
        tit=strcat('Top ',num2str(PrintingRange(1)+i*SpacingInt),'(Corr>=',num2str(CurrentCorrVal),')and last',num2str(PrintingRange(1)+i*SpacingInt),'(Corr<= ',num2str(anti_CurrentCorrVal),')Correlated Cells of ',answer3{1,1});
        f=strcat('Top ',num2str(PrintingRange(1)+i*SpacingInt),'and last',num2str(PrintingRange(1)+i*SpacingInt),'Correlated Cells of ',answer3{1,1});
    else
        if PrintingRange(2)-PrintingRange(1)>0
            tit=strcat('Top ',num2str(PrintingRange(1)),' to ',num2str(PrintingRange(1)+i*SpacingInt),' Correlated Cells of ',answer3{1,1},', Corr>=',num2str(CurrentCorrVal));
            f=strcat('Top ',num2str(PrintingRange(1)),' to ',num2str(PrintingRange(1)+i*SpacingInt),' Correlated Cells of ',answer3{1,1});
        else
            tit=strcat('Top ',num2str(PrintingRange(1)-(loopNo-i)*SpacingInt),' to ',num2str(PrintingRange(1)),' Correlated Cells of ',answer3{1,1},', Corr<=',num2str(CurrentCorrVal));
            f=strcat('Top ',num2str(PrintingRange(1)-(loopNo-i)*SpacingInt),' to ',num2str(PrintingRange(1)),' Correlated Cells of ',answer3{1,1});
        end
    end
    view(0,90);
    title(tit);
    savefig(Fig3,strcat(path,'\',f,'.fig'));
    
    if strcmp(answer5{1,1},'3')
        Fig3 = figure;
        posB = handles.spPos(handles.labelVector,:);
        posR = handles.spPos(selectedCell,:);
        plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
        hold on;
        plot3(posR(:,1),posR(:,2),posR(:,3),'ro', 'markerfacecolor','r','markersize',3);
        if strcmp(answer2{1,1},'1')
            plot3(handles.spPos(InputCellNo,1),handles.spPos(InputCellNo,2),handles.spPos(InputCellNo,3),'bo','linewidth',3,'Markersize',10);
        end
        grid on
        axis equal
        xlabel('X','Color','k');
        ylabel('Y','Color','k');
        zlabel('Z','Color','k');
        if PrintingRange(2)-PrintingRange(1)>0
            tit=strcat('Top ',num2str(PrintingRange(1)),' to ',num2str(PrintingRange(1)+i*SpacingInt),' Correlated Cells of ',answer3{1,1},', Corr>=',num2str(CurrentCorrVal));
            f=strcat('Top ',num2str(PrintingRange(1)),' to ',num2str(PrintingRange(1)+i*SpacingInt),' Correlated Cells of ',answer3{1,1});
        else
            tit=strcat('Top ',num2str(PrintingRange(1)-(loopNo-i)*SpacingInt),' to ',num2str(PrintingRange(1)),' Correlated Cells of ',answer3{1,1},', Corr<=',num2str(CurrentCorrVal));
            f=strcat('Top ',num2str(PrintingRange(1)-(loopNo-i)*SpacingInt),' to ',num2str(PrintingRange(1)),' Correlated Cells of ',answer3{1,1});
        end
        view(0,90);
        title(tit);
        savefig(Fig3,strcat(path,'\',f,'.fig'));
    end
    
end


% current_roi = handles.roiMaster.Value;
% handles.roi(current_roi).members(handles.SelectedCorrCellNo) = 1;
% AddCellNo=length(handles.SelectedCorrCellNo)
% handles = display_roiListbox(handles, handles.roi(current_roi).name);
% guidata(hObject,handles);   
% set(handles.ROICellNo,'String',num2str(sum(handles.roi(current_roi).members)));


% --- Executes on button press in ClearColor.
function ClearColor_Callback(hObject, eventdata, handles)
% hObject    handle to ClearColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allRoiCells=handles.roi(1).members;
for i=1:size(handles.roi,2)
    allRoiCells=allRoiCells+handles.roi(i).members;
end
allRoiCells(allRoiCells>0)=1;
allRoiCells=logical(allRoiCells);
handles.DispColor(allRoiCells)=30;
handles=plot_pos_maps(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ShowAllROICheck.
function ShowAllROICheck_Callback(hObject, eventdata, handles)
% hObject    handle to ShowAllROICheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ShowAllROICheck
handles = plot_pos_maps(handles);
handles = plot_slice_maps(handles);
handles = plot_fts_in_slice(handles);

% --- Executes on key press with focus on roiMaster and none of its controls.
function roiMaster_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to roiMaster (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
current_roi = handles.roiMaster.Value;
if isfield(handles,'DispColor')==0
    handles.DispColor=ones(size(handles.fts,1),1);
end
    
%note: keep everything inside the if/elif statements to prevent weird
%selection behavior. 
if strcmp(eventdata.Character,'a')%add ROI Master
    answer0=inputdlg('What is the name of the ROI Category you want?');
    newROIname=answer0{1,1};
    newROISize=size(handles.roi,2)+1;
    handles.roi(newROISize).name= newROIname;
    handles.roi(newROISize).members=zeros(size(handles.spPos,1),1);
    handles.roi(newROISize).members=logical(handles.roi(newROISize).members);
    set(handles.roiMaster,'String',{handles.roi.name});
    assignin('base','roi',handles.roi);

elseif strcmp(eventdata.Key, 'delete')
    handles.roi(current_roi) = [];
    update_roiMasterList(handles);
    handles = display_roiListbox(handles);
    
elseif strcmp(eventdata.Character,'v')
    answer1=inputdlg('What color do you want this ROI group to be plotted? please input 0 to 100 digit.');

    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=str2num(answer1{1,1});
    assignin('base','DispColor',handles.DispColor);
%     
%     namestr = cellstr(get(hObject, 'String'));
%     validx = get(hObject, 'Value');
%     newstr = regexprep(namestr{validx}, '"red"','"green"');
%     namestr{validx} = newstr;
%     set(hObject, 'String', namestr);

    guidata(hObject,handles);
    
elseif strcmp(eventdata.Character,'p')%Rythmic Neurons-spike like
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=76;
elseif strcmp(eventdata.Character,'m')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=85;
elseif strcmp(eventdata.Character,'b')%Rythmic Neurons-wave like
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=60;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'y')%Responding Cells
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=18;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'o')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=10;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'g')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=35;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'r')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=0;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'n')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=50;
    assignin('base','DispColor',handles.DispColor);
    %---------------------------
elseif strcmp(eventdata.Character,'1')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=10;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'2')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=20;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'3')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=30;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'4')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=40;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'5')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=50;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'6')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=60;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'7')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=70;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'8')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=80;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'9')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=90;
    assignin('base','DispColor',handles.DispColor);
elseif strcmp(eventdata.Character,'0')
    roi = find(handles.roi(current_roi).members);
    handles.DispColor(roi)=0;
    assignin('base','DispColor',handles.DispColor);
end
    guidata(hObject,handles);

function [htmlname]=regexprep(Entry)
htmlname = sprintf('<HTML><BODY bgcolor="%s">%s', 'red', Entry);


% --- Executes on button press in LoadSequence.
function LoadSequence_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[SeriesFile,SeriesFilePath] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(SeriesFilePath);
SeriesFileDirectory= [SeriesFilePath,SeriesFile];
SeriesFile=load(SeriesFileDirectory);

plot(1:length(SeriesFile.StimTimeSeries),SeriesFile.StimTimeSeries,'Parent',handles.StimTSPlot);
xLimits = get(handles.dffPlot,'XLim');
set(handles.StimTSPlot,'XLim',xLimits);
set(handles.StimTSPlot,'XColor',[1 1 1]);
guidata(hObject,handles);


% --- Executes on button press in ColorSelection.
function ColorSelection_Callback(hObject, eventdata, handles)
% hObject    handle to ColorSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = plot_pos_maps(handles);
handles = plot_slice_maps(handles);
handles = plot_fts_in_slice(handles);


% --- Executes on button press in ExportROIData.
function ExportROIData_Callback(hObject, eventdata, handles)
% hObject    handle to ExportROIData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_roi = handles.roiMaster.Value;
roi = find(handles.roi(current_roi).members);

% CurrentRoi=find(handles.roi.members);
% assignin('base','CurrentRoiList',CurrentRoi);
spPos=handles.spPos(roi,:);
fluorescence_time_series=handles.fts(roi,:);
spRadiiXYZ=handles.spRadiiXYZ(roi,:);
% cellSegmentation=handles.cellSegmentation(:,roi);
Sc=handles.Sc;
ROIList=roi;

if isfield(handles,'dFF')
    dFF=handles.dFF(roi,:);
end
[f2,path2] = uiputfile('ROIData.mat');
if isfield(handles,'dFF')==0
    save([path2,f2],'ROIList','spPos','fluorescence_time_series','spRadiiXYZ','Sc');
else
    save([path2,f2],'ROIList','spPos','fluorescence_time_series','dFF','spRadiiXYZ','Sc');
end


% --- Executes on button press in RandColor.
function RandColor_Callback(hObject, eventdata, handles)
% hObject    handle to RandColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

allRoiCells=handles.roi(1).members;
for i=1:size(handles.roi,2)
    allRoiCells=allRoiCells+handles.roi(i).members;
end
allRoiCells(allRoiCells>0)=1;
allRoiCells=logical(allRoiCells);

assignin('base','allroi',allRoiCells);
assignin('base','DispColor',handles.DispColor);
assignin('base','spPos',handles.spPos);
handles.DispColor(allRoiCells)=zeros(size(handles.spPos,1),1);

ROIMasterNo=size(handles.roi,2);
StepNo=(100/ROIMasterNo);
tempvector=[0:StepNo:100];
tempvector=tempvector';

assignin('base','tempvector',tempvector);
for i1=1:ROIMasterNo
    currentRoiCell=handles.roi(i1).members;
    handles.DispColor(currentRoiCell)=tempvector(i1);
end
handles=plot_pos_maps(handles);
guidata(hObject,handles);


% --- Executes on button press in gen_time_slice_movies.
function gen_time_slice_movies_Callback(hObject, eventdata, handles)
% hObject    handle to gen_time_slice_movies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

% set(handles.status_bar,'String','Loading...');
% set(handles.slider2,'enable','off');
% set(handles.slider1,'enable','off');
% drawnow;

currentDirectory = pwd;

stackDirectory = uigetdir(pwd, 'Select .klb directory to generate time slice movies');
cd(stackDirectory);

mkdir('tsView');
lsStackFiles= dir ('*.klb');
stackdata = readImage(lsStackFiles(1).name);

zSize = size(stackdata,3);
tSize = numel(stackDirectory);

% fs=filesep;
cd('tsView');

for z=1:zSize
    stackName = sprintf('ts%03d',z);
    tiffName{z} = [stackName,'.tif'];

    if exist(tiffName{z},'file');
       delete(tiffName{z});
    end
    
%     tiffObj{z} = Tiff(tiffName,'a');

end

for t=1:tSize
    fullStack = readImage(['..\',lsStackFiles(t).name]);
    
    if t==1 %detect whether uint8 or uint16
        if isa(fullStack,'uint8')
            BitsPerSample = 8;
        elseif isa(fullStack,'uint16')
            BitsPerSample = 16; 
        end
    end
    
    for m=1:zSize
%         objt = tiffObj{m};

        tiffObj = Tiff(tiffName{m},'a');
        
        tiffObj.setTag('Photometric',Tiff.Photometric.LinearRaw);
        tiffObj.setTag('BitsPerSample',BitsPerSample);
        tiffObj.setTag('ImageWidth',size(fullStack,2));
        tiffObj.setTag('ImageLength',size(fullStack,1));
        tiffObj.setTag('SamplesPerPixel',1);
        tiffObj.setTag('Compression',Tiff.Compression.PackBits);
        tiffObj.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    
        tiffObj.write( fullStack(:,:,m) );        
    
    end
    
    if t==tSize
        writeDirectory(tiffObj);
        tiffObj.close();
    end
    
    disp(sprintf('Processed frame %03d',t));
end    

currZ = floor(zSize/2);
stackName = sprintf('ts%03d',currZ);

cd (currentDirectory);
guidata(hObject,handles);


% --- Executes on button press in load_slice_movie.
function load_slice_movie_Callback(hObject, eventdata, handles)
% hObject    handle to load_slice_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
[file,path] = uigetfile('*.tif', 'Select tsView file');

fs = filesep;
disp('Reading file...');
tic; handles.currView = readImage([path,fs,file]); toc;
handles = reset_slider_handles(handles);

guidata(hObject, handles);

function handles = reset_slider_handles(handles)
    if handles.PlotSelect.Value == 3
        if ~isempty(handles.zStack)
        zMx = size(handles.zStack,3);
        currScale = handles.cellSelector.Value/handles.cellSelector.Max;
        set(handles.cellSelector,'Max',zMx);
        set(handles.cellSelector,'Min',1);
        set(handles.cellSelector,'SliderStep',[1/zMx, 5/zMx]); 
        
        newValue = round(currScale * zMx);
        if newValue < 1 || newValue > handles.cellSelector.Max
            newValue = round(handles.cellSelector.Max/2);
        end
        handles.cellSelector.Value = newValue;
        
        end

        if ~isempty(handles.currView)
        tMx = size(handles.currView,3);
        currScale = handles.sliceSelector.Value/handles.sliceSelector.Max;
        set(handles.sliceSelector,'Max', tMx);
        set(handles.sliceSelector,'SliderStep',[ 1/tMx, 5/tMx]);
        newValue = round(currScale * tMx);
        if newValue < 1 || newValue > handles.sliceSelector.Max
            newValue = round(handles.sliceSelector.Max/2);
        end
        handles.sliceSelector.Value = newValue;
        end
    else
        currScale = handles.sliceSelector.Value/handles.sliceSelector.Max;
        set(handles.sliceSelector,'Max',length(handles.linearRange)-1);
        set(handles.sliceSelector,'SliderStep',[ 1/(length(handles.linearRange)-1), 0.1]);
        newValue = round(currScale * length(handles.linearRange)-1);
        if newValue < 1 || newValue > handles.sliceSelector.Max
            newValue = round(handles.sliceSelector.Max/2);
        end
        handles.sliceSelector.Value = newValue;
        handles.currentSlice = newValue;
        
        handles.inSlice = findInSlice(handles);
        inSlice = handles.inSlice;
        nCells = sum(inSlice);
        currScale = handles.cellSelector.Value/handles.cellSelector.Max;
        set(handles.cellSelector, 'Max', nCells);
        set(handles.cellSelector,'SliderStep',[1/nCells, 5*(1/nCells)]);
        
        newValue = round(currScale * nCells);
        if newValue < 1 || newValue > handles.cellSelector.Max
            newValue = round(handles.cellSelector.Max/2);
        end
        handles.cellSelector.Value = newValue;
        
    end


% --- Executes on button press in load_z_stack.
function load_z_stack_Callback(hObject, eventdata, handles)
% hObject    handle to load_z_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[file,path] = uigetfile({'*.klb','*.tif'}, 'Select reference z-stack to load');
tic; handles.zStack = readImage([path,filesep,file]); toc; 
handles = reset_slider_handles(handles);
handles.calling_function = 'Load Z';
plot_slice_maps(handles);
guidata(hObject, handles);



function caxis0_Callback(hObject, eventdata, handles)
% hObject    handle to caxis0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caxis0 as text
%        str2double(get(hObject,'String')) returns contents of caxis0 as a double


% --- Executes during object creation, after setting all properties.
function caxis0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caxis0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function caxis1_Callback(hObject, eventdata, handles)
% hObject    handle to caxis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caxis1 as text
%        str2double(get(hObject,'String')) returns contents of caxis1 as a double


% --- Executes during object creation, after setting all properties.
function caxis1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caxis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in throw_roi_button.
function throw_roi_button_Callback(hObject, eventdata, handles)
% hObject    handle to throw_roi_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
zLvl = round( get(handles.cellSelector, 'Value') );
roi_z = handles.spPos(handles.roi.members,3);
inZ = round( roi_z/handles.microns_per_z ) == zLvl;
handles.roiListbox.Value = find(inZ);
handles.cell_roi_list = [];
guidata(hObject, handles);

function varargout = data_explorer_gui(varargin)
% DATA_EXPLORER_GUI MATLAB code for data_explorer_gui.fig
%      DATA_EXPLORER_GUI, by itself, creates a new DATA_EXPLORER_GUI or raises the existing
%      singleton*.
%
%      H = DATA_EXPLORER_GUI returns the handle to a new DATA_EXPLORER_GUI or the handle to
%      the existing singleton*.
%
%      DATA_EXPLORER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATA_EXPLORER_GUI.M with the given input arguments.
%
%      DATA_EXPLORER_GUI('Property','Value',...) creates a new DATA_EXPLORER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before data_explorer_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to data_explorer_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help data_explorer_gui

% Last Modified by GUIDE v2.5 24-Apr-2018 12:15:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @data_explorer_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @data_explorer_gui_OutputFcn, ...
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

% --- Executes just before data_explorer_gui is made visible.
function data_explorer_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to data_explorer_gui (see VARARGIN)

% Choose default command line output for data_explorer_gui
handles.output = hObject;

handles.sliceThickness = 20;

homepath = getenv('HOME');
fs = filesep;
handles.guiHome_functions = [homepath,fs,'Dropbox',fs,'CellAnalysisUtilities',fs,'LightSheet' fs 'light_sheet_data_explorer_gui', fs, 'functions'];
addpath(handles.guiHome_functions);
handles.stackAcqFreq = 1; %initialized with a dummy value

handles.cellSelect = 1; %for plotting on slice axis
handles.cellSelect_idx = 1; %for keeping track of absolute cell identity

handles.roi = struct;

handles.im=[]
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using data_explorer_gui.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

function handles=loadDataset(filepath,handles,eventdata)
%switches data set
load(filepath);
handles.fts = fluorescence_time_series;
handles.spPos = spPos;
handles.Sc = Sc;

if exist('dFF','var')
    handles.dFF = dFF;
else
    handles.dFF = NaN;
end

if exist('roi','var')
    handles.roi=roi;
    handles = display_roiListbox(handles, roi(1).name);
else
    handles.roi.name = 'default';
    handles.roi.members = zeros(size(spPos,1),1,'logical');
    handles = display_roiListbox(handles, 'default');
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
handles.labelVector(1:10:end) = 1;

handles = plot_pos_maps(handles);
handles = plot_fts_in_slice(handles);


function handles = get_linear_range(handles)
handles.linearRange = floor( min(handles.spPos(:,1)) ):handles.sliceThickness: ceil( max(handles.spPos(:,1)) );
handles.linearRange = [handles.linearRange inf];

function handles = set_label_density(handles, density)
%density is a fraction between 0 and 1
spacing = round(density^-1);
handles.labelVector = zeros(size(handles.spPos,1),1,'logical');
handles.labelVector(1:spacing:end) = 1;

function set_slice_thickness_Callback(hObject, eventdata, handles)
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

% UIWAIT makes data_explorer_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function handles = plot_fts_in_slice(handles)
%plots the slice specific scatter plot and sets the current cell based on
%cellSelector slide

% inSlice = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & cellRecruitment;
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
hold on;
if isfield(handles,'circlePoint')
       delete(handles.circlePoint );
end
handles.circlePoint = plot3 (handles.spPos(cellSelect,1), handles.spPos(cellSelect,2),handles.spPos(cellSelect,3),'ro','linewidth',2);
handles=plotfts(handles,cellSelect);

%updates dFF plot with current cell
% axes(handles.dffPlot), cla %this bit of code will generate a warning if the cell has changed but trial has not. not important. 
% % yyaxis left
% hold on
% 
% % plot(handles.fts(cellSelect,:),'color',[0 1 0],'linewidth',2);
% plotTraces_GUI(handles.fts(cellSelect,:), 100, 10, handles.dffPlot);
% 
% xlabel('Time (s)');
% ylabel('\Delta F/F');
% title('Fluorescence Time Series', 'color', [1 1 1]);
% set(handles.dffPlot,'FontUnits','normalized','FontSize',0.06);
% handles.dffPlot = gca;

function handles=plotfts(handles, cellSelect)
%plots fluorescent time series on the time series plot
axes(handles.dffPlot), cla %this bit of code will generate a warning if the cell has changed but trial has not. not important. 
% yyaxis left
hold on

% plot(handles.fts(cellSelect,:),'color',[0 1 0],'linewidth',2);
if get(handles.raw_trace_selector,'Value')
    plotTraces_GUI(handles.fts(cellSelect,:), 100, 10, handles.dffPlot);
    colorbar('off');
    ylabel('Intensity - Mean');
elseif get(handles.dFF_traces,'Value')
    plotTraces_GUI(handles.dFF(cellSelect,:), 1, 10, handles.dffPlot);
    colorbar('off');
    ylabel('\Delta F/F');
elseif get(handles.heatmap_selector,'Value')
    toplot = handles.fts(cellSelect,:);
    toplot = bsxfun(@minus, toplot, median(toplot,2));
    imagesc(handles.dffPlot, toplot );
    axis tight
    colorbar;  %colorbar('on');
    colormap(jet);
elseif get(handles.dFF_heatmap,'Value')
    toplot = handles.dFF(cellSelect,:);
    imagesc(handles.dffPlot, toplot );
    axis tight
    colorbar;  %colorbar('on');
    colormap(jet);
end

xlabel('Frames');
title('Fluorescence Time Series', 'color', [1 1 1]);
set(handles.dffPlot,'FontUnits','normalized','FontSize',0.06);
handles.dffPlot = gca;

function handles=plot_pos_maps(handles)
%plots a recruitment map onto sliceAx given the current slice selection
%NOTE: CORRECT SLICE IS K:K+1 AND NUMBER OF SLICES IS
%LENGTH(RCSEGMENTATION)- 1

handles.cellClicker = handles.sliceAx.ButtonDownFcn;

axes(handles.sliceAx), cla
sK = handles.currentSlice;
% RCSegmentation = handles.linearRange;
% all_inSlice = handles.spPos(:,1) >  RCSegmentation(sK-1) & handles.spPos(:,1) < RCSegmentation(sK);
all_inSlice = findInSlice(handles);
allDots = handles.spPos(all_inSlice,:);
handles.inSlice = all_inSlice; 
current_roi = get(handles.roiMaster,'Value');
roi_in_slice = handles.roi(current_roi).members(all_inSlice);
plot3(allDots(~roi_in_slice,1),allDots(~roi_in_slice,2),allDots(~roi_in_slice,3),'.','color',[0.2 0.2 0.2],'hittest','off');
plot3(allDots(roi_in_slice,1),allDots(roi_in_slice,2),allDots(roi_in_slice,3),'.','color','r','hittest','off');
set(gca,'TickLength',[0,0],'XTick',[],'YTick',[],'ZTick',[],'FontUnits','normalized','FontSize',0.05);
grid on
axis equal
view(90,0);
title(sprintf('Slice %2.0f',sK-1),'color',[1,1,1]);

handles.sliceAx.ButtonDownFcn = handles.cellClicker;
%%%%%
%%%%%

if handles.radiobutton_Slice.Value
axes(handles.slicePosMap), cla
hold on
posB = handles.spPos(handles.labelVector & ~handles.inSlice,:);
posR = handles.spPos(handles.labelVector & handles.inSlice,:);
plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
plot3(posR(:,1),posR(:,2),posR(:,3),'r.');
grid on
axis equal

elseif handles.radiobutton_roi.Value
current_roi = handles.roiMaster.Value;    
axes(handles.slicePosMap), cla
hold on
posB = handles.spPos(handles.labelVector,:);
posR = handles.spPos(handles.roi(current_roi).members,:);
plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
plot3(posR(:,1),posR(:,2),posR(:,3),'ro', 'markerfacecolor','r','markersize',5);
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

handles.SliceMap_gLine{1} = plot([RCSegmentation(sK),RCSegmentation(sK)], currentLim, 'g');
handles.SliceMap_gLine{2} = plot([RCSegmentation(sK-1),RCSegmentation(sK-1)], currentLim, 'g');
handles.SliceMap_gLine{3} = plot3([RCSegmentation(sK),RCSegmentation(sK)], currentLim, [zLim(2),zLim(2)], 'g');
handles.SliceMap_gLine{4} = plot3([RCSegmentation(sK-1),RCSegmentation(sK-1)], currentLim, [zLim(2), zLim(2)], 'g');
end
set(gca,'TickLength',[0,0],'XTick',[],'YTick',[],'ZTick',[],'FontUnits','normalized','FontSize',0.05);
% plot(handles.spPos(:,1),handles.spPos(:,2),'k.');
% axis equal
% axis tight
% 
% currentLim = get(gca,'YLim');
% 
% if isfield(handles,'SliceMap_gLine')
%    for gl = 1:numel(handles.SliceMap_gLine)
%        delete(handles.SliceMap_gLine{gl} );
%    end
% end
% axis equal
% axis tight
% 
% handles.SliceMap_gLine{1} = plot([RCSegmentation(sK),RCSegmentation(sK)], currentLim, 'g');
% handles.SliceMap_gLine{2} = plot([RCSegmentation(sK-1),RCSegmentation(sK-1)], currentLim, 'g');

function inSlice = findInSlice(handles)
%consisent framework for identifying which cells are in slice
sK = handles.currentSlice;
RCSegmentation = handles.linearRange;
inSlice = handles.spPos(:,1) >  RCSegmentation(sK-1) & handles.spPos(:,1) < RCSegmentation(sK);

% --- Outputs from this function are returned to the command line.
function varargout = data_explorer_gui_OutputFcn(hObject, eventdata, handles)
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
handles.cellSelect = round( get(hObject,'Value') );
handles = plot_fts_in_slice(handles);
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

[mD,cellSelected] = min(distance_to_list);
handles.cellSelector.Value = cellSelected;

handles = cellSelector_Callback(handles.cellSelector,eventdata,handles);
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
[currentfile,path] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(path);
handles.currentfile = [path,currentfile];
handles = loadDataset(handles.currentfile,handles,eventdata);
guidata(hObject,handles);


% --- Executes on selection change in roiListbox.
function roiListbox_Callback(hObject, eventdata, handles)
% hObject    handle to roiListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiListbox


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
    handles = display_roiListbox(handles, handles.roi(current_roi).name);
    guidata(hObject,handles);

elseif strcmp(eventdata.Character,'d')
    %     handles.roi(handles.roiListbox.Value) = [];
    toDeleteCell = handles.roiListbox.String( handles.roiListbox.Value);
    toDelete = str2num(toDeleteCell{1});
    handles.roi(current_roi).members(toDelete) = 0;
    handles = display_roiListbox(handles, handles.roi(current_roi).name);
    guidata(hObject,handles);

end

% --- Executes on button press in radiobutton_Slice.
function radiobutton_Slice_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_pos_maps(handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton_Slice


% --- Executes on button press in ExportROI.
function ExportROI_Callback(hObject, eventdata, handles)
% hObject    handle to ExportROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi = handles.roi;
[f,path] = uiputfile('roi.mat');
save([path,f],'roi');

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
handles.roi = roi;
% set(handles.roiListbox,'String',num2cell(handles.roi))
% handles.roiListbox.Max = length(handles.roi);
% handles.roiListbox.Value = 1;
% guidata(hObject,handles);

set(handles.roiListbox,'String',num2cell(find(handles.roi)))
handles.roiListbox.Max = sum(handles.roi);
handles.roiListbox.Value = 1;
guidata(hObject,handles);


% --- Executes on button press in radiobutton_roi.
function radiobutton_roi_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_pos_maps(handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton_roi


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
pxyz = handles.spPos(handles.roi,:);
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

function handles = display_roiListbox(handles, roiname)
%displays current members of roi listed by roiname
roiIdx = strcmp({handles.roi.name}, roiname);
set(handles.roiListbox,'String',num2cell(find(handles.roi(roiIdx).members)))
handles.roiListbox.Max = sum(handles.roi(roiIdx).members);
handles.roiListbox.Value = 1;

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
handles = display_roiListbox(handles, handles.roi( hObject.Value).name);
handles = plot_pos_maps(handles);
handles = plot_fts_in_slice(handles);

% --- Executes during object creation, after setting all properties.
function roiMaster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiMaster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

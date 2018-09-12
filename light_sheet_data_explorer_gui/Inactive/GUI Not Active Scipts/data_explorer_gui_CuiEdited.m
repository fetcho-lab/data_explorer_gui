function varargout = data_explorer_gui_CuiEdited(varargin)
% DATA_EXPLORER_GUI_CUIEDITED MATLAB code for data_explorer_gui_CuiEdited.fig
%      DATA_EXPLORER_GUI_CUIEDITED, by itself, creates a new DATA_EXPLORER_GUI_CUIEDITED or raises the existing
%      singleton*.
%
%      H = DATA_EXPLORER_GUI_CUIEDITED returns the handle to a new DATA_EXPLORER_GUI_CUIEDITED or the handle to
%      the existing singleton*.
%
%      DATA_EXPLORER_GUI_CUIEDITED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATA_EXPLORER_GUI_CUIEDITED.M with the given input arguments.
%
%      DATA_EXPLORER_GUI_CUIEDITED('Property','Value',...) creates a new DATA_EXPLORER_GUI_CUIEDITED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before data_explorer_gui_CuiEdited_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to data_explorer_gui_CuiEdited_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help data_explorer_gui_CuiEdited

% Last Modified by GUIDE v2.5 30-Apr-2018 17:20:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @data_explorer_gui_CuiEdited_OpeningFcn, ...
                   'gui_OutputFcn',  @data_explorer_gui_CuiEdited_OutputFcn, ...
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

% --- Executes just before data_explorer_gui_CuiEdited is made visible.
function data_explorer_gui_CuiEdited_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to data_explorer_gui_CuiEdited (see VARARGIN)

% Choose default command line output for data_explorer_gui_CuiEdited
handles.output = hObject;

handles.sliceThickness = 20;

homepath = getenv('HOME');
fs = filesep;
handles.guiHome_functions = [homepath,fs,'Dropbox',fs,'CellAnalysisUtilities',fs,'LightSheet' fs 'light_sheet_data_explorer_gui', fs, 'functions'];
addpath(handles.guiHome_functions);
handles.stackAcqFreq = 1; %initialized with a dummy value

handles.cellSelect = 1; %for plotting on slice axis
handles.cellSelect_idx = 1; %for keeping track of absolute cell identity

handles.roi = [];

handles.im=[]
% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using data_explorer_gui_CuiEdited.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
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
assignin('base', 'handles',handles);
guidata(hObject,handles);


function handles=loadDataset(filepath,handles,eventdata)
%switches data set
load(filepath);
handles.fts = fluorescence_time_series;
handles.spPos = spPos;
handles.Sc = Sc;

handles.linearRange = floor( min(handles.spPos(:,1)) ):handles.sliceThickness: ceil( max(handles.spPos(:,1)) );
handles.linearRange = [handles.linearRange inf];

handles.currentSlice = ceil(length(handles.linearRange)/2);

set(handles.sliceSelector,'Max',length(handles.linearRange)-1);
set(handles.sliceSelector,'Value',handles.currentSlice);
set(handles.sliceSelector,'SliderStep',[ 1/(length(handles.linearRange)-1), 0.1]);

handles.labelVector = zeros(size(handles.spPos,1),1,'logical');
handles.labelVector(1:10:end) = 1;

handles = plot_pos_maps(handles);
handles = plot_fts(handles);
handels = plot_HeatMap(handles);


function handles = plot_HeatMap(handles)
axes(handles.BrainHeatMap);
clims=[min(min(handles.fts)), max(max(handles.fts))];%Set colorbar scale of heatmap here
imagesc(handles.fts(:,:),clims);
colorbar
colormap jet
ax = gca;
ax.XAxis.Color = 'white';
ax.YAxis.Color = 'white';
CB2=colorbar
CB2.Color='w';
% CB2.Label.String='Max Fluorescence Intensity of Cells (Dots)';


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

% UIWAIT makes data_explorer_gui_CuiEdited wait for user response (see UIRESUME)
% uiwait(handles.figure1);
function handles = plot_fts(handles)
%plots the slice specific scatter plot and sets the current cell based on
%cellSelector slide

% inSlice = alxPos(:,1) >  RCSegmentation(sK-1) & alxPos(:,1) < RCSegmentation(sK) & cellRecruitment;
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

handles.f_inSlice = find(inSlice);
% randtoSelect = randperm( sum(inSlice) );
cellSelect = handles.f_inSlice( handles.cellSelect);
handles.cellSelect_idx = cellSelect;
%updates circled point on coronal slice plot
axes(handles.sliceAx);
hold on;
if isfield(handles,'circlePoint')
       delete(handles.circlePoint );
end
handles.circlePoint = plot3 (handles.spPos(cellSelect,1), handles.spPos(cellSelect,2),handles.spPos(cellSelect,3),'ro','linewidth',2);

%updates dFF plot with current cell
axes(handles.dffPlot), cla %this bit of code will generate a warning if the cell has changed but trial has not. not important. 
% yyaxis left
hold on

plot(handles.fts(cellSelect,:),'color',[0 1 0],'linewidth',2);
%----Cui modified-----
xlabel('Time (s)','color','w');
% ylabel('\Delta F/F','color','w'); %Cui: Our data ploted here is not dF/F!
%it's actually just fluorescence intensity values.
ax = gca;
ax.XAxis.Color = 'white';
ax.YAxis.Color = 'white';
%----Cui modify end---
title('Fluorescence Time Series', 'color', [1 1 1]);
set(handles.dffPlot,'FontUnits','normalized','FontSize',0.06);
handles.dffPlot = gca;


function handles=plot_pos_maps(handles)
%plots a recruitment map onto sliceAx given the current slice selection
%NOTE: CORRECT SLICE IS K:K+1 AND NUMBER OF SLICES IS
%LENGTH(RCSEGMENTATION)- 1

handles.cellClicker = handles.sliceAx.ButtonDownFcn;
assignin('base','cellclicker2',handles.cellClicker);

axes(handles.sliceAx), cla
sK = handles.currentSlice;
% RCSegmentation = handles.linearRange;
% all_inSlice = handles.spPos(:,1) >  RCSegmentation(sK-1) & handles.spPos(:,1) < RCSegmentation(sK);
all_inSlice = findInSlice(handles);
allDots = handles.spPos(all_inSlice,:);
handles.inSlice = all_inSlice; 
plot3(allDots(:,1),allDots(:,2),allDots(:,3),'.','color',[0.2 0.2 0.2],'hittest','off');
set(gca,'TickLength',[0,0],'XTick',[],'YTick',[],'ZTick',[],'FontUnits','normalized','FontSize',0.05);
grid on
axis equal
view(90,0);
title(sprintf('Slice %2.0f',sK-1),'color',[1,1,1]);

handles.sliceAx.ButtonDownFcn = handles.cellClicker;

assignin('base','cellclicker22',handles.cellClicker);

%%%%%
%%%%%

if handles.radiobutton_Slice.Value
axes(handles.slicePosMap), cla
hold on
%----Part modified by Cui----
% posB = handles.spPos(handles.labelVector & ~handles.inSlice,:);
% posR = handles.spPos(handles.labelVector & handles.inSlice,:);
% plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
% plot3(posR(:,1),posR(:,2),posR(:,3),'r.');
%----Part modify by Cui end---------
%----Cui has modified code above to:-----
posAllCell=handles.spPos;
handles.ColorOfDots=max(handles.fts');%set color of dots to each cells' highest fluo value.
scatter3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),10,handles.ColorOfDots,'.','hittest','off');
colormap jet;
CB1=colorbar
CB1.Color='w';
CB1.Label.String='Max Fluorescence Intensity of Cells (Dots)';
posR_C = handles.spPos( handles.inSlice,:);
plot3(posR_C(:,1),posR_C(:,2),posR_C(:,3),'k.');
%---Cui's modify end---------
% assignin('base','posR',posR);

grid on
axis equal

elseif handles.radiobutton_roi.Value
    
axes(handles.slicePosMap), cla
hold on
posB = handles.spPos(handles.labelVector,:);
posR = handles.spPos(handles.roi,:);
plot3(posB(:,1),posB(:,2),posB(:,3),'k.');
plot3(posR(:,1),posR(:,2),posR(:,3),'ro','linewidth',2);
grid on
axis equal
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
function varargout = data_explorer_gui_CuiEdited_OutputFcn(hObject, eventdata, handles)
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
handles = plot_fts(handles); 
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
handles = plot_fts(handles);
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
assignin('base','clickpos',clickPosition);
assignin('base','position_list',position_list);
distance_to_list = pdist2(single( clickPosition ),position_list);

[mD,cellSelected] = min(distance_to_list);
assignin('base','cellselected',cellSelected);
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
handles = guidata(hObject);

if strcmp(eventdata.Character,'a')
    handles.roi = unique( [handles.roi, handles.cellSelect_idx] );
elseif strcmp(eventdata.Character,'d')
    handles.roi(handles.roiListbox.Value) = [];
end
set(handles.roiListbox,'String',num2cell(handles.roi))
handles.roiListbox.Max = length(handles.roi);
handles.roiListbox.Value = 1;
guidata(hObject,handles);


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
set(handles.roiListbox,'String',num2cell(handles.roi))
handles.roiListbox.Max = length(handles.roi);
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
objPos = objSelected.GetX


% --- Executes on mouse press over axes background.
function slicePosMap_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to slicePosMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);

position_list = single( handles.spPos );
clickPosition = eventdata.IntersectionPoint;
assignin('base','clickpos',clickPosition);
assignin('base','position_list',position_list);
distance_to_list = pdist2(single( clickPosition ),position_list);

[mD,cellSelected] = min(distance_to_list);
assignin('base','cellselected',cellSelected);
handles.cellSelector.Value = cellSelected;
assignin('base','clickpos',clickPosition);
handles = cellSelector_Callback(handles.cellSelector,eventdata,handles);
PlotSelectedCellOnPosMap(hObject, eventdata, handles);
guidata(hObject,handles);

function PlotSelectedCellOnPosMap(hObject, eventdata, handles)
hold on;
axes(handles.slicePosMap);
if isfield(handles,'circlePoint')
       delete(handles.circlePoint );
end
handles.circlePoint = plot3 (handles.spPos(handles.cellSelector.Value,1), handles.spPos(handles.cellSelector.Value,2),handles.spPos(handles.cellSelector.Value,3),'ro','linewidth',2);



% --- Executes on button press in AddAllCell.
function AddAllCell_Callback(hObject, eventdata, handles)
% hObject    handle to AddAllCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.roi=[1:size(handles.fts,1)];
set(handles.roiListbox,'String',num2cell(handles.roi))
handles.roiListbox.Max = length(handles.roi);
handles.roiListbox.Value = 1;
assignin('base','handles',handles);
guidata(hObject, handles);


% --- Executes on button press in ExportROIMeanFluo.
function ExportROIMeanFluo_Callback(hObject, eventdata, handles)
% hObject    handle to ExportROIMeanFluo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ROIMeanFluo=mean(handles.fts(handles.roi,:));
ROIList=handles.roi;
[f2,path2] = uiputfile('ROIMeanFluo&ROIList.mat');
save([path2,f2],'ROIMeanFluo','ROIList');
guidata(hObject, handles);


% --- Executes on button press in Export_Fluo_of_ROI.
function Export_Fluo_of_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to Export_Fluo_of_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ROIFluo=handles.fts(handles.roi,:);
ROIList=handles.roi;
ROIList=ROIList';
% assignin('base','Datatowrite',DataToWrite);
% assignin('base','handles',handles);
[f2,path2] = uiputfile('ROIFluo&ROIList.mat');
save([path2,f2],'ROIFluo','ROIList');
guidata(hObject, handles);



% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.CellNoToFind=str2double(get(hObject,'String'));
axes(handles.slicePosMap);
hold on
if isfield(handles,'CirclePoint');
    delete(handles.CirclePoint);
end
handles.CirclePoint=plot3(handles.spPos(handles.CellNoToFind,1), handles.spPos(handles.CellNoToFind,2),handles.spPos(handles.CellNoToFind,3),'ro','linewidth',2);

set(handles.XVal, 'String', num2str(handles.spPos(handles.CellNoToFind,1)));
set(handles.YVal, 'String', num2str(handles.spPos(handles.CellNoToFind,2)));
set(handles.ZVal, 'String', num2str(handles.spPos(handles.CellNoToFind,3)));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


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


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.roi=[handles.roi;handles.f_inSlice];
assignin('base','f_inSlice',handles.f_inSlice);
assignin('base','Currentroi',handles.roi);
set(handles.roiListbox,'String',num2cell(handles.roi))
handles.roiListbox.Max = length(handles.roi);
handles.roiListbox.Value = 1;
% assignin('base','handles',handles);
guidata(hObject, handles);


% function PlotTopCorr(AllMaxcor, AllMaxcorno, AllPos, AllFluo,handles);
% %------------Adjustable--------------
% handles.Dis_thre=[0,999];
% handles.Cor_thre=0.995;
% handles.WidthRangeForLines=[2,8];
% handles.RandomLineColor=1;
% handles.FluoLineColorMaxVal=max(max(AllFluo));%set FluoLineColorMaxVal here!!!!!!
% %---------------------------------------
% 
% %---------The code following seems should have the same plotting function
% %but the result is werid---------------
% % DispList1=find(AllMaxcor>handles.Cor_thre);
% % DispList2=AllMaxcorno(DispList1);
% % for i=1:length(DispList1)
% %     tempx=[AllPos(DispList1,1),AllPos(DispList2,1)];
% %     tempy=[AllPos(DispList1,2),AllPos(DispList2,2)];
% %     tempz=[AllPos(DispList1,3),AllPos(DispList2,3)];
% %     linewidth=(((AllMaxcor(DispList1(i))-handles.Cor_thre)/(1-handles.Cor_thre))*(handles.WidthRangeForLines(2)-handles.WidthRangeForLines(1)))+handles.WidthRangeForLines(1);
% %     Corrplot=plot3(tempx, tempy, tempz, '-','LineWidth',linewidth);
% % end
% %------------------------------------
% Plotted_Correlation_No=0;
%         for i3=1:size(AllFluo,1)
%             tempx=[AllPos(i3,1),AllPos(AllMaxcorno(i3),1)];
%             tempy=[AllPos(i3,2),AllPos(AllMaxcorno(i3),2)];
%             tempz=[AllPos(i3,3),AllPos(AllMaxcorno(i3),3)];
%             dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
%             if dist>handles.Dis_thre(1)
%                 if dist<handles.Dis_thre(2)
%                     if AllMaxcor(i3)>handles.Cor_thre
%                         if AllMaxcor(i3)<=1
%                         linewidth=(((AllMaxcor(i3)-handles.Cor_thre)/(1-handles.Cor_thre))*(handles.WidthRangeForLines(2)-handles.WidthRangeForLines(1)))+handles.WidthRangeForLines(1);
%                             if handles.RandomLineColor==1
%                                 Corrplot=plot3(tempx, tempy, tempz, '-','LineWidth',linewidth);
%                             else
%                                 colorR=max(AllFluo(i3,:))/handles.FluoLineColorMaxVal;
%                                 if colorR>1
%                                     colorR=1;
%                                 end
%                                 colorG=max(AllFluo(AllMaxcorno(i3),:))/handles.FluoLineColorMaxVal;
%                                 if colorG>1
%                                     colorG=1;
%                                 end
%                                 colorB=0;
%                                 Corrplot=plot3(tempx, tempy, tempz, '-','Color',[colorR, colorG, colorB],'LineWidth',linewidth);
%                                 hold on
%                                 Plotted_Correlation_No=Plotted_Correlation_No+1
%                             end
%                         end
%                     end
%                 end
%             end
%         end 
% guidata(hObject,handles);

% --- Executes on button press in DeleteAllROI.
function DeleteAllROI_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteAllROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.roi=[];
set(handles.roiListbox,'String',num2cell(handles.roi))
handles.roiListbox.Max = length(handles.roi);
handles.roiListbox.Value = 1;
guidata(hObject,handles);


% --------------------------------------------------------------------
function CorrelationMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CorrelationMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FindTopCorrCellForAllCells_Callback(hObject, eventdata, handles)
% hObject    handle to FindTopCorrCellForAllCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.AllMaxcorval,handles.AllMaxcorno,handles.Allcor]=TopCorrCellInGroup_corrcoef(handles.fts)
% axes(handles.slicePosMap);
% PlotTopCorr(handles.AllMaxcor,handles.AllMaxcorno,handles.spPos,handles.fts,handles);
% h=dialog('name','Report')
assignin('base','AllMaxcorno',handles.AllMaxcorno);
assignin('base','AllMaxcorval',handles.AllMaxcorval);
msgbox({'All Correlation Calculated!' 'Stored in handles.AllMaxcorval and handles.AllMaxcorno' 'You can also check AllMaxcorval & AllMaxcorno in the base workspace'});
guidata(hObject,handles);


% --------------------------------------------------------------------
function TopCorrForROIs_Callback(hObject, eventdata, handles)
% hObject    handle to TopCorrForROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.ROIMaxcorval,handles.ROIMaxcorno,handles.ROIcor]=TopCorrCellInGroup_corrcoef(handles.fts(handles.roi,:));
% assignin('base','AllMaxCorno1',handles.AllMaxcorno);
handles.ROIMaxcorno=handles.roi(handles.ROIMaxcorno);
assignin('base','ROIMaxcorno',handles.ROIMaxcorno);
assignin('base','ROIMaxcorval',handles.ROIMaxcorval);
assignin('base','ROI',handles.roi);
% axes(handles.slicePosMap);
% PlotTopCorr(handles.AllMaxcor,handles.AllMaxcorno,handles.spPos,handles.fts,handles);
msgbox({'All Correlation Calculated!' 'Stored in handles.ROIMaxcorval and handles.ROIMaxcorno' ...
    'You can also check ROIMaxcorval1 & ROIMaxcorno & ROI in the base workspace'});
guidata(hObject,handles);


% --------------------------------------------------------------------
function FindCorrCellOfFluoSeries_Callback(hObject, eventdata, handles)
% hObject    handle to FindCorrCellOfFluoSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[SeriesFile,SeriesFilePath] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(SeriesFilePath);
SeriesFileDirectory= [SeriesFilePath,SeriesFile];
SeriesFile=load(SeriesFileDirectory);
answer1=inputdlg({'Top k correlated cells will be stored. K number?'},'Input K'); 
k=str2num(answer1{1,1});
TargetSeries=SeriesFile.ROIMeanFluo;
TargetSeriesNo=size(TargetSeries,1);
[handles.TopKCorrVal, handles.TopKCorrCellNo, handles.AllCorrToSeries]=MaxkCorr_corrcoef(handles.fts, TargetSeries, TargetSeriesNo, k);
assignin('base','TopKCorrVal1',handles.TopKCorrVal);
assignin('base','TopKCorrCellNo1',handles.TopKCorrCellNo);
assignin('base','AllCorrToSeries1',handles.AllCorrToSeries);

msgbox({'Calculation Done!' 'Top k correlated cell no & correlation value are stored in handles.ROIMaxcorval and handles.ROIMaxcorno' ...
    'The matrix of correlation between each cell and input series is stored in handles.AllCorrToSeries' ...
    'You can also check TopKCorrVal1 & TopKCorrCellNo1 & AllCorrToSeries1 in the base workspace'});


% --------------------------------------------------------------------
function FindCorrCellsOfRoiMean_Callback(hObject, eventdata, handles)
% hObject    handle to FindCorrCellsOfRoiMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer=inputdlg({'Top k correlated cells will be stored. K number?'},'Input K'); 
k=str2num(answer{1,1});
% RoiCellNo=handles.roi';
handles.ROIMeanFluo=mean(handles.fts(handles.roi,:));
InputSeries=handles.ROIMeanFluo;
InputSeriesNo=size(InputSeries,1);
[handles.TopKCorrVal, handles.TopKCorrCellNo, handles.AllCorrToSeries]=MaxkCorr_corrcoef(handles.fts, InputSeries, InputSeriesNo, k);
assignin('base','TopKCorrVal2',handles.TopKCorrVal);
assignin('base','TopKCorrCellNo2',handles.TopKCorrCellNo);
assignin('base','AllCorrToSeries2',handles.AllCorrToSeries);
assignin('base','ROIMeanFluo',handles.ROIMeanFluo);
msgbox({'Calculation Done!' 'Top k correlated cell no & correlation value are stored in handles.ROIMaxcorval and handles.ROIMaxcorno' ...
    'The matrix of correlation between each cell and input series is stored in handles.AllCorrToSeries' 'Mean FLuo of ROI is stored in handles.ROIMeanFluo' ...
    'You can also check TopKCorrVal2 & TopKCorrCellNo2 & AllCorrToSeries2 & ROIMeanFluo in the base workspace'});


% --------------------------------------------------------------------
function TopCorrViewer_Callback(hObject, eventdata, handles)
% hObject    handle to TopCorrViewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
TopCorrViewer_v2;
guidata(hObject,handles);

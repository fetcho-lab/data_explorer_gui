function varargout = ClusterViewer(varargin)
% CLUSTERVIEWER MATLAB code for ClusterViewer.fig
%      CLUSTERVIEWER, by itself, creates a new CLUSTERVIEWER or raises the existing
%      singleton*.
%
%      H = CLUSTERVIEWER returns the handle to a new CLUSTERVIEWER or the handle to
%      the existing singleton*.
%
%      CLUSTERVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLUSTERVIEWER.M with the given input arguments.
%
%      CLUSTERVIEWER('Property','Value',...) creates a new CLUSTERVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClusterViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClusterViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClusterViewer

% Last Modified by GUIDE v2.5 09-Jul-2018 15:42:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClusterViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ClusterViewer_OutputFcn, ...
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


% --- Executes just before ClusterViewer is made visible.
function ClusterViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClusterViewer (see VARARGIN)

% Choose default command line output for ClusterViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ClusterViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% 
% varargout{1}=handles.AllMaxcor;
% varargout{2}=handles.AllMaxcorno;
% delete(handles.figure1);

function handles=loadDataset(filepath,handles,eventdata)
%switches data set
load(filepath);
handles.fts = fluorescence_time_series;
handles.spPos = spPos;
handles.spRadiiXYZ=spRadiiXYZ;
handles.cellSegmentation=cellSegmentation;
handles.extractParams=extractParams;
handles.Sc=Sc;

if exist('dFF','var');
    handles.dFF=dFF;
    print1='dFF Detected'
end
if exist('BackGroundFTS','var');
    handles.BackGroundFTS=BackGroundFTS;
    print1='BackGroundFTS Detected'
end
if exist('FluoAfterSubtraction','var');
    handles.FluoAfterSubtraction=FluoAfterSubtraction;
    print1='FluoAfterSubtraction Detected'
end


assignin('base','AllPos',handles.spPos);
global AllPosition %Defined a global here!
AllPosition=spPos;
handles.Sc = Sc;
handles.SizeOfDots=20;%Set default dot size here

% assignin('base', 'handles0',handles);
axes(handles.PosMap), cla
plot_pos_maps(handles.spPos,handles.fts,handles.SizeOfDots);
hold on

function plot_pos_maps(posAllCell,fts,SizeOfDots)
% posAllCell=handles.spPos;
ColorOfDots=max(fts');
PosMap=scatter3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),SizeOfDots,ColorOfDots,'.');%,'PickableParts','none'
xlabel(gca,'X');
ylabel(gca,'Y');
zlabel(gca,'Z');
set(gca,'Xcolor','w');
set(gca,'Ycolor','w');
set(gca,'Zcolor','w');
set(gca,'Gridcolor','k');
% PosMap.XColor='w';
% set(gca,'buttondownfcn',@clicky);
% assignin('base','x',posAllCell(:,1));
% plot3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),'.','color',[0.2 0.2 0.2],'hittest','off');
set(gca,'Visible','on');
colormap jet;
% caxis(lim1);
CB1=colorbar;
CB1.Color='w';
CB1.Label.String='Fluorescence Intensity of Cells (Dots)';
grid on


% --- Outputs from this function are returned to the command line.
function varargout = ClusterViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
[handles.filename,path] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(path);
handles.currentfile = [path,handles.filename];
handles = loadDataset(handles.currentfile,handles,eventdata);
handles.DeleteList=[];
handles.CurrentCorrCalculated=0;
guidata(hObject,handles);


function SizeOfDots_Callback(hObject, eventdata, handles)
% hObject    handle to SizeOfDots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeOfDots as text
%        str2double(get(hObject,'String')) returns contents of SizeOfDots as a double
handles.SizeOfDots=str2double(get(hObject,'String'))

axes(handles.PosMap), cla
plot_pos_maps(handles.spPos,handles.fts,handles.SizeOfDots);
% assignin('base','Sizeofdots',handles.SizeofDots);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SizeOfDots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeOfDots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowOriginalCellPlot.
function ShowOriginalCellPlot_Callback(hObject, eventdata, handles)
% hObject    handle to ShowOriginalCellPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.PosMap), cla
plot_pos_maps(handles.spPos,handles.fts,handles.SizeOfDots);

handles=guidata(hObject);

function CellNoToFind_Callback(hObject, eventdata, handles)
% hObject    handle to CellNoToFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNoToFind as text
%        str2double(get(hObject,'String')) returns contents of CellNoToFind as a double
handles.CellNoToFind=str2double(get(hObject,'String'));
axes(handles.PosMap);
hold on
if isfield(handles,'CirclePoint');
    delete(handles.CirclePoint);
end
assignin('base','spPos',handles.spPos);
assignin('base','CellNoToFind',handles.CellNoToFind);
handles.CirclePoint=plot3(handles.spPos(handles.CellNoToFind,1), handles.spPos(handles.CellNoToFind,2),handles.spPos(handles.CellNoToFind,3),'ro','linewidth',3,'Markersize',10);

set(handles.XVal, 'String', num2str(handles.spPos(handles.CellNoToFind,1)));
set(handles.YVal, 'String', num2str(handles.spPos(handles.CellNoToFind,2)));
set(handles.ZVal, 'String', num2str(handles.spPos(handles.CellNoToFind,3)));
guidata(hObject,handles);


% --- Executes on button press in LoadCluster.
function LoadCluster_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[f,path] = uigetfile('*.mat');
load([path,f]);

handles.CorrValToSeries=CorrValToSeries;
handles.CorrCellNoToSeries=CorrCellNoToSeries;
handles.AllCorrToSeries=AllCorrToSeries;
handles.AllPValtoSeries=AllPValtoSeries;

% --- Executes on button press in ExportCluster.
function ExportCluster_Callback(hObject, eventdata, handles)
% hObject    handle to ExportCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function CorrThreEdit_Callback(hObject, eventdata, handles)
% hObject    handle to CorrThreEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CorrThreEdit as text
%        str2double(get(hObject,'String')) returns contents of CorrThreEdit as a double


% --- Executes during object creation, after setting all properties.
function CorrThreEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrThreEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CorrTypeSelection.
function CorrTypeSelection_Callback(hObject, eventdata, handles)
% hObject    handle to CorrTypeSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrTypeSelection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrTypeSelection


% --- Executes during object creation, after setting all properties.
function CorrTypeSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrTypeSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

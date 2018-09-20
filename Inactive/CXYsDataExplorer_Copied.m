function varargout = CXYsDataExplorer_Copied(varargin)
% CXYSDATAEXPLORER_COPIED MATLAB code for CXYsDataExplorer_Copied.fig
%      CXYSDATAEXPLORER_COPIED, by itself, creates a new CXYSDATAEXPLORER_COPIED or raises the existing
%      singleton*.
%
%      H = CXYSDATAEXPLORER_COPIED returns the handle to a new CXYSDATAEXPLORER_COPIED or the handle to
%      the existing singleton*.
%
%      CXYSDATAEXPLORER_COPIED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CXYSDATAEXPLORER_COPIED.M with the given input arguments.
%
%      CXYSDATAEXPLORER_COPIED('Property','Value',...) creates a new CXYSDATAEXPLORER_COPIED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CXYsDataExplorer_Copied_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CXYsDataExplorer_Copied_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CXYsDataExplorer_Copied

% Last Modified by GUIDE v2.5 09-May-2018 12:03:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CXYsDataExplorer_Copied_OpeningFcn, ...
                   'gui_OutputFcn',  @CXYsDataExplorer_Copied_OutputFcn, ...
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


% --- Executes just before CXYsDataExplorer_Copied is made visible.
function CXYsDataExplorer_Copied_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CXYsDataExplorer_Copied (see VARARGIN)

% Choose default command line output for CXYsDataExplorer_Copied
handles.output = hObject;

assignin('base','varargin',varargin)
handles.Currentpath=varargin{1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CXYsDataExplorer_Copied wait for user response (see UIRESUME)
uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = CXYsDataExplorer_Copied_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.figure1);

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_data_file_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_file_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = guidata(hObject);
handles = guidata(hObject);
[currentfile,handles.Currentpath] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(handles.Currentpath);
handles.currentfile = [handles.Currentpath,currentfile];
handles.CurrentTPNo=1; %Set initial Time point here!!
handles = loadDataset(handles.currentfile,handles,eventdata);
handles.TotalTPNo=size(handles.fts,2);
set(handles.TotalTP,'String',num2str(handles.TotalTPNo));
set(handles.CurrentTP,'String',num2str(handles.CurrentTPNo));
handles.CurrentCorrCalculated=0;
guidata(hObject,handles);

% --- Executes on button press in LoadFromParent.
function LoadFromParent_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFromParent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [currentfile,path] = uigetfile('.mat','Please select fluorescence time series .mat file');
assignin('base','handles0',handles);
cd(handles.Currentpath);
handles.currentfile = [path,handles.Currentpath];
handles.CurrentTPNo=1; %Set initial Time point here!!
handles = loadDataset(handles.currentfile,handles,eventdata);
handles.TotalTPNo=size(handles.fts,2);
set(handles.TotalTP,'String',num2str(handles.TotalTPNo));
set(handles.CurrentTP,'String',num2str(handles.CurrentTPNo));
handles.CurrentCorrCalculated=0;
guidata(hObject,handles);


function handles=loadDataset(filepath,handles,eventdata)
%switches data set
load(filepath);
handles.fts = fluorescence_time_series;
handles.spPos = spPos;
% global AllPosition %Defined a global here!
% AllPosition=spPos;
handles.Sc = Sc;
handles.SizeOfDots=20;%Set default dot size here
set(handles.FileName,'String',filepath);
assignin('base', 'handles0',handles);
handles = PlotCurrentTimepointMap(handles);
handles = PlotHeatMap(handles)

function handles=PlotCurrentTimepointMap(handles)
axes(handles.TPMap), cla
 
posAllCell=handles.spPos;
handles.ColorOfDots=handles.fts(:,handles.CurrentTPNo);

TPMap=scatter3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),handles.SizeOfDots,handles.ColorOfDots,'.');%,'PickableParts','none'
hold on

xlabel(handles.TPMap,'X');
ylabel(handles.TPMap,'Y');
zlabel(handles.TPMap,'Z');
set(gca,'Xcolor','w');
set(gca,'Ycolor','w');
set(gca,'Zcolor','w');
set(gca,'Gridcolor','k');
% TPMap.XColor='w';
% set(gca,'buttondownfcn',@clicky);
% assignin('base','x',posAllCell(:,1));
% plot3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),'.','color',[0.2 0.2 0.2],'hittest','off');
set(handles.TPMap,'Visible','on');
colormap jet;
lim1=[min(min(handles.fts)),max(max(handles.fts))];
caxis(lim1);

CB1=colorbar;
CB1.Color='w';
CB1.Label.String='Max Fluorescence Intensity of Cells (Dots)';

grid on

function handles=PlotHeatMap(handles)
axes(handles.HeatMap),cla;
% Fluomin=min(handles.fts,[],2);
% Fluomax=max(handles.fts,[],2);
MaxFluoOfAllCells=max(max(handles.fts));
MinFluoOfAllCells=min(min(handles.fts));
clims=[MinFluoOfAllCells,MaxFluoOfAllCells];
imagesc(handles.fts(1:end,:),clims);
colorbar
colormap jet
ylim=[0 size(handles.fts,1)];


% --- Executes on slider movement.
function TimeLineSlider_Callback(hObject, eventdata, handles)
% hObject    handle to TimeLineSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
CurrentSliderVal=get(hObject,'Value');
% assignin('base','CurrentTP',handles.CurrentTPNo);
handles.CurrentTPNo=round(CurrentSliderVal*(handles.TotalTPNo-1))+1;
handles = PlotCurrentTimepointMap(handles);
set(handles.CurrentTP,'String',num2str(handles.CurrentTPNo));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function TimeLineSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeLineSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function CurrentTP_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentTP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrentTP as text
%        str2double(get(hObject,'String')) returns contents of CurrentTP as a double

handles.CurrentTPNo=str2double(get(hObject,'String'));
handles = PlotCurrentTimepointMap(handles);
set(handles.TimeLineSlider,'Value',(handles.CurrentTPNo-1)/(handles.TotalTPNo-1));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function CurrentTP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentTP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MoveForward.
function MoveForward_Callback(hObject, eventdata, handles)
% hObject    handle to MoveForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CurrentTPNo=handles.CurrentTPNo+1;
handles = PlotCurrentTimepointMap(handles);
set(handles.TimeLineSlider,'Value',(handles.CurrentTPNo-1)/(handles.TotalTPNo-1));
set(handles.CurrentTP,'String',num2str(handles.CurrentTPNo));
guidata(hObject,handles);

% --- Executes on button press in MoveBackward.
function MoveBackward_Callback(hObject, eventdata, handles)
% hObject    handle to MoveBackward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CurrentTPNo=handles.CurrentTPNo-1;
handles = PlotCurrentTimepointMap(handles);
set(handles.TimeLineSlider,'Value',(handles.CurrentTPNo-1)/(handles.TotalTPNo-1));
set(handles.CurrentTP,'String',num2str(handles.CurrentTPNo));
guidata(hObject,handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.CellNoToFind=str2double(get(hObject,'String'));
axes(handles.TPMap);
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
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end

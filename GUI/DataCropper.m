function varargout = DataCropper(varargin)
% DATACROPPER MATLAB code for DataCropper.fig
%      DATACROPPER, by itself, creates a new DATACROPPER or raises the existing
%      singleton*.
%
%      H = DATACROPPER returns the handle to a new DATACROPPER or the handle to
%      the existing singleton*.
%
%      DATACROPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATACROPPER.M with the given input arguments.
%
%      DATACROPPER('Property','Value',...) creates a new DATACROPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DataCropper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DataCropper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DataCropper

% Last Modified by GUIDE v2.5 02-Jul-2018 17:09:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataCropper_OpeningFcn, ...
                   'gui_OutputFcn',  @DataCropper_OutputFcn, ...
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


% --- Executes just before DataCropper is made visible.
function DataCropper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataCropper (see VARARGIN)

% Choose default command line output for DataCropper
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DataCropper wait for user response (see UIRESUME)
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

% handles.XLowerLimit=0;
% set(handles.XLowerLimitEdit,'String',num2str(0));
% handles.XUpperLimit=max(handles.spPos(:,1));
% set(handles.XUpperLimitEdit,'String',num2str(handles.XUpperLimit));

% handles.YLowerLimit=0;
% set(handles.YLowerLimitEdit,'String',num2str(0));
% handles.YUpperLimit=max(handles.spPos(:,2));
% set(handles.YUpperLimitEdit,'String',num2str(handles.YUpperLimit));

handles.ZLowerLimit=0;
set(handles.ZLowerLimitEdit,'String',num2str(0));
handles.ZUpperLimit=max(handles.spPos(:,3));
set(handles.ZUpperLimitEdit,'String',num2str(handles.ZUpperLimit));


handles.LowerTimeRange=1;
set(handles.TimeRangeLowerLimit,'String',num2str(1));

handles.UpperTimeRange=size(handles.fts,2);
set(handles.TimeRangeUpperLimit,'String',num2str(handles.UpperTimeRange));


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
function varargout = DataCropper_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in ListToDelet.
function ListToDelet_Callback(hObject, eventdata, handles)
% hObject    handle to ListToDelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListToDelet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListToDelet


% --- Executes during object creation, after setting all properties.
function ListToDelet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListToDelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XLowerLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to XLowerLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XLowerLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of XLowerLimitEdit as a double
handles.XLowerLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function XLowerLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XLowerLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function XUpperLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to XUpperLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XUpperLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of XUpperLimitEdit as a double

handles.XUpperLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function XUpperLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XUpperLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YLowerLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to YLowerLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YLowerLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of YLowerLimitEdit as a double
handles.YLowerLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function YLowerLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLowerLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YUpperLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to YUpperLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YUpperLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of YUpperLimitEdit as a double
handles.YUpperLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function YUpperLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YUpperLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZLowerLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ZLowerLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZLowerLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of ZLowerLimitEdit as a double
handles.ZLowerLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ZLowerLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZLowerLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZUpperLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ZUpperLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZUpperLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of ZUpperLimitEdit as a double
handles.ZUpperLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ZUpperLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZUpperLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CropData.
function CropData_Callback(hObject, eventdata, handles)
% hObject    handle to CropData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CellNoListBeforeCrop=[1:size(handles.fts,1)];
xlist=find(handles.spPos(:,1)>=handles.XLowerLimit&handles.spPos(:,1)<=handles.XUpperLimit);
ylist=find(handles.spPos(:,2)>=handles.YLowerLimit&handles.spPos(:,2)<=handles.YUpperLimit);
zlist=find(handles.spPos(:,3)>=handles.ZLowerLimit&handles.spPos(:,3)<=handles.ZUpperLimit);
handles.CroppedCellList=intersect(intersect(xlist,ylist),zlist);
% assignin('base','CroppedCellList',handles.CroppedCellList);


%------------Variables to write-----------
% assignin('base','fts1',handles.fts);
fluorescence_time_series=handles.fts(handles.CroppedCellList,handles.LowerTimeRange:handles.UpperTimeRange);
% assignin('base','fts2',fluorescence_time_series);
% handles.LowerTimeRange
% handles.UpperTimeRange
spPos=handles.spPos(handles.CroppedCellList,:);
spRadiiXYZ=handles.spRadiiXYZ(handles.CroppedCellList,:);
cellSegmentation=handles.cellSegmentation(:,handles.CroppedCellList);

extractParams=handles.extractParams;
Sc=handles.Sc;

if isfield(handles,'dFF')
    dFF=handles.dFF(handles.CroppedCellList,handles.LowerTimeRange:handles.UpperTimeRange);
    print2='dFF detected, prepared to write...'
end
if isfield(handles,'BackGroundFTS')
    BackGroundFTS=handles.BackGroundFTS;
end
if isfield(handles,'FluoAfterSubtraction')
    FluoAfterSubtraction=handles.FluoAfterSubtraction(handles.CroppedCellList,handles.LowerTimeRange:handles.UpperTimeRange);
end
    
    
axes(handles.PosMap), cla
plot_pos_maps(spPos,fluorescence_time_series,handles.SizeOfDots);

XLimit=[handles.XLowerLimit,handles.XUpperLimit];
YLimit=[handles.YLowerLimit,handles.YUpperLimit];
ZLimit=[handles.ZLowerLimit,handles.ZUpperLimit];
CroppedCellNoList=handles.CroppedCellList;
Note='CroppedCellNoList can be used to check the original Cell No in original Data files';

[f1,path1] = uiputfile(strcat(handles.filename,'-Cropped.mat'));
save([path1,f1],'fluorescence_time_series','spPos','spRadiiXYZ','cellSegmentation','extractParams','Sc','XLimit','YLimit','ZLimit','CroppedCellNoList','Note');

if isfield(handles,'dFF')
    save([path1,f1],'dFF','-append');
    print2='dFF wrote!'
end
if isfield(handles,'BackGroundFTS')
    save([path1,f1],'BackGroundFTS','-append');
end
if isfield(handles,'FluoAfterSubtraction')
    save([path1,f1],'FluoAfterSubtraction','-append');
end
msgbox('All Variables Wrote!');
% uiresume(handles.figure1);
guidata(hObject,handles);

% --- Executes on selection change in DeleteCellList.
function DeleteCellList_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteCellList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DeleteCellList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DeleteCellList


% --- Executes during object creation, after setting all properties.
function DeleteCellList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeleteCellList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SizeOfDots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CellNoToAdd_Callback(hObject, eventdata, handles)
% hObject    handle to CellNoToAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNoToAdd as text
%        str2double(get(hObject,'String')) returns contents of CellNoToAdd as a double
handles.CellNoToAdd=str2double(get(hObject,'String'));
handles.DeleteList=[handles.DeleteList;handles.CellNoToAdd];
set(handles.DeleteCellList,'String',handles.DeleteList);
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


% --- Executes on button press in ShowCroppedPlot.
function ShowCroppedPlot_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCroppedPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.PosMap), cla;
plot_pos_maps(handles.spPosCropped,handles.ftsCropped,handles.SizeOfDots);


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function TimeRangeLowerLimit_Callback(hObject, eventdata, handles)
% hObject    handle to TimeRangeLowerLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeRangeLowerLimit as text
%        str2double(get(hObject,'String')) returns contents of TimeRangeLowerLimit as a double

handles.LowerTimeRange=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function TimeRangeLowerLimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeRangeLowerLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TimeRangeUpperLimit_Callback(hObject, eventdata, handles)
% hObject    handle to TimeRangeUpperLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeRangeUpperLimit as text
%        str2double(get(hObject,'String')) returns contents of TimeRangeUpperLimit as a double
handles.UpperTimeRange=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function TimeRangeUpperLimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeRangeUpperLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

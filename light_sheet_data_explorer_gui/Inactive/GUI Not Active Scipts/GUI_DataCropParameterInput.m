function varargout = GUI_DataCropParameterInput(varargin)
% GUI_DATACROPPARAMETERINPUT MATLAB code for GUI_DataCropParameterInput.fig
%      GUI_DATACROPPARAMETERINPUT, by itself, creates a new GUI_DATACROPPARAMETERINPUT or raises the existing
%      singleton*.
%
%      H = GUI_DATACROPPARAMETERINPUT returns the handle to a new GUI_DATACROPPARAMETERINPUT or the handle to
%      the existing singleton*.
%
%      GUI_DATACROPPARAMETERINPUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DATACROPPARAMETERINPUT.M with the given input arguments.
%
%      GUI_DATACROPPARAMETERINPUT('Property','Value',...) creates a new GUI_DATACROPPARAMETERINPUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_DataCropParameterInput_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_DataCropParameterInput_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_DataCropParameterInput

% Last Modified by GUIDE v2.5 15-May-2018 18:52:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_DataCropParameterInput_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_DataCropParameterInput_OutputFcn, ...
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


% --- Executes just before GUI_DataCropParameterInput is made visible.
function GUI_DataCropParameterInput_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_DataCropParameterInput (see VARARGIN)

% Choose default command line output for GUI_DataCropParameterInput
handles.output = hObject;
handles.dFF = NaN;
handles.fts = varargin{1};
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_DataCropParameterInput wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_DataCropParameterInput_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% guidata(hObject);
% varargout{1} =  getappdata(hObject,'dFF');
varargout{1} =  handles.CroppedDataPointList;
delete(handles.figure1);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusbox,'String','Status: Computing...');
% drawnow;
% stack_frequency = str2num( get(handles.XLowerRange,'String') );
% YLowerRange = str2num( get(handles.YLowerRange, 'String') );
% ZLowerRange = str2num(get(handles.ZLowerRange, 'String') ); 
% lowpass_on = get(handles.lowfilt_check,'Value');
% handles.dFF = lsExplorer_Compute_dFF(handles.fts, stack_frequency, YLowerRange, ZLowerRange, lowpass_on);
% % handles.dFF=1;
set(handles.statusbox,'String','Status: Done!');
% setappdata(hObject,'dFF', handles.dFF);


guidata(hObject, handles);


function XLowerRange_Callback(hObject, eventdata, handles)
% hObject    handle to XLowerRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XLowerRange as text
%        str2double(get(hObject,'String')) returns contents of XLowerRange as a double


% --- Executes during object creation, after setting all properties.
function XLowerRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XLowerRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YLowerRange_Callback(hObject, eventdata, handles)
% hObject    handle to YLowerRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YLowerRange as text
%        str2double(get(hObject,'String')) returns contents of YLowerRange as a double


% --- Executes during object creation, after setting all properties.
function YLowerRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLowerRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZLowerRange_Callback(hObject, eventdata, handles)
% hObject    handle to ZLowerRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZLowerRange as text
%        str2double(get(hObject,'String')) returns contents of ZLowerRange as a double


% --- Executes during object creation, after setting all properties.
function ZLowerRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZLowerRange (see GCBO)
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
% delete(hObject);
if isequal(get(hObject, 'waitstatus'), 'waiting')
    uiresume(hObject);
else
    delete(hObject);
end


% --- Executes on button press in load_pushbutton.
function load_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, path] = uigetfile('*.mat', 'Select file containing variable named dFF');
load([path,filename]);
handles.dFF = dFF;
setappdata(hObject,'dFF', handles.dFF);
guidata(hObject,handles);


% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, path] = uiputfile('*.mat','Select file to save (adds, does not overwrite if file exists)');
CroppedDataPointList = handles.CroppedDataPointList;
if isnan(dFF)
    error('dFF has not been computed!');
end
savepath = [path,filename];
if exist(savepath,'file')
    save(savepath, 'CroppedDataPointList', '-append');
else
    save(savepath, 'CroppedDataPointList', '-v7.3');
end



function XUpperRange_Callback(hObject, eventdata, handles)
% hObject    handle to XUpperRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of XUpperRange as text
%        str2double(get(hObject,'String')) returns contents of XUpperRange as a double


% --- Executes during object creation, after setting all properties.
function XUpperRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XUpperRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YUpperRange_Callback(hObject, eventdata, handles)
% hObject    handle to YUpperRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YUpperRange as text
%        str2double(get(hObject,'String')) returns contents of YUpperRange as a double


% --- Executes during object creation, after setting all properties.
function YUpperRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YUpperRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ZUpperRange_Callback(hObject, eventdata, handles)
% hObject    handle to ZUpperRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZUpperRange as text
%        str2double(get(hObject,'String')) returns contents of ZUpperRange as a double


% --- Executes during object creation, after setting all properties.
function ZUpperRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZUpperRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

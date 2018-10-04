function varargout = GUI_dFF_Awesome(varargin)
% GUI_DFF_AWESOME MATLAB code for GUI_dFF_Awesome.fig
%      GUI_DFF_AWESOME, by itself, creates a new GUI_DFF_AWESOME or raises the existing
%      singleton*.
%
%      H = GUI_DFF_AWESOME returns the handle to a new GUI_DFF_AWESOME or the handle to
%      the existing singleton*.
%
%      GUI_DFF_AWESOME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_DFF_AWESOME.M with the given input arguments.
%
%      GUI_DFF_AWESOME('Property','Value',...) creates a new GUI_DFF_AWESOME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_dFF_Awesome_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_dFF_Awesome_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_dFF_Awesome

% Last Modified by GUIDE v2.5 02-May-2018 14:00:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_dFF_Awesome_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_dFF_Awesome_OutputFcn, ...
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


% --- Executes just before GUI_dFF_Awesome is made visible.
function GUI_dFF_Awesome_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_dFF_Awesome (see VARARGIN)

% Choose default command line output for GUI_dFF_Awesome
handles.output = hObject;
handles.dFF = NaN;
handles.fts = varargin{1};
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_dFF_Awesome wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_dFF_Awesome_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% guidata(hObject);
% varargout{1} =  getappdata(hObject,'dFF');
varargout{1} =  handles.dFF;
delete(handles.figure1);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusbox,'String','Status: Computing...');
drawnow;
stack_frequency = str2num( get(handles.frequency,'String') );
window_size = str2num( get(handles.window_size, 'String') );
basepercent = str2num(get(handles.basepercent, 'String') ); 
lowpass_on = get(handles.lowfilt_check,'Value');
handles.dFF = lsExplorer_Compute_dFF(handles.fts, stack_frequency, window_size, basepercent, lowpass_on);
% handles.dFF=1;
set(handles.statusbox,'String','Status: Done!');
setappdata(hObject,'dFF', handles.dFF);
guidata(hObject, handles);


function frequency_Callback(hObject, eventdata, handles)
% hObject    handle to frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frequency as text
%        str2double(get(hObject,'String')) returns contents of frequency as a double


% --- Executes during object creation, after setting all properties.
function frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function window_size_Callback(hObject, eventdata, handles)
% hObject    handle to window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of window_size as text
%        str2double(get(hObject,'String')) returns contents of window_size as a double


% --- Executes during object creation, after setting all properties.
function window_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function basepercent_Callback(hObject, eventdata, handles)
% hObject    handle to basepercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basepercent as text
%        str2double(get(hObject,'String')) returns contents of basepercent as a double


% --- Executes during object creation, after setting all properties.
function basepercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basepercent (see GCBO)
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


% --- Executes on button press in lowfilt_check.
function lowfilt_check_Callback(hObject, eventdata, handles)
% hObject    handle to lowfilt_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lowfilt_check


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
dFF = handles.dFF;
if isnan(dFF)
    error('dFF has not been computed!');
end
savepath = [path,filename];
if exist(savepath,'file')
    save(savepath, 'dFF', '-append');
else
    save(savepath, 'dFF', '-v7.3');
end


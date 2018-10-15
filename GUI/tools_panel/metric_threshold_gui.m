function varargout = metric_threshold_gui(varargin)
% METRIC_THRESHOLD_GUI MATLAB code for metric_threshold_gui.fig
%      METRIC_THRESHOLD_GUI, by itself, creates a new METRIC_THRESHOLD_GUI or raises the existing
%      singleton*.
%
%      H = METRIC_THRESHOLD_GUI returns the handle to a new METRIC_THRESHOLD_GUI or the handle to
%      the existing singleton*.
%
%      METRIC_THRESHOLD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in METRIC_THRESHOLD_GUI.M with the given input arguments.
%
%      METRIC_THRESHOLD_GUI('Property','Value',...) creates a new METRIC_THRESHOLD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before metric_threshold_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to metric_threshold_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help metric_threshold_gui

% Last Modified by GUIDE v2.5 15-Oct-2018 15:22:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @metric_threshold_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @metric_threshold_gui_OutputFcn, ...
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


% --- Executes just before metric_threshold_gui is made visible.
function metric_threshold_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to metric_threshold_gui (see VARARGIN)

% Choose default command line output for metric_threshold_gui
handles.output = hObject;
input_handles = varargin{1};

current_metric = input_handles.metric_listbox.Value;
handles.metric_name = input_handles.metric(current_metric).description;
handles.metric_name_textbox.String = handles.metric_name;
handles.metric_value = input_handles.metric(current_metric).value;

handles.threshold_passed = zeros(size(handles.metric_value,1), 1, 'logical');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes metric_threshold_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = metric_threshold_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.threshold_passed;
varargout{2} = handles.threshold_name;
delete(handles.figure1);



function greater_than_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to greater_than_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of greater_than_threshold as text
%        str2double(get(hObject,'String')) returns contents of greater_than_threshold as a double


% --- Executes during object creation, after setting all properties.
function greater_than_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to greater_than_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function less_than_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to less_than_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of less_than_threshold as text
%        str2double(get(hObject,'String')) returns contents of less_than_threshold as a double


% --- Executes during object creation, after setting all properties.
function less_than_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to less_than_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lower_threshold = str2num(handles.greater_than_threshold.String);
upper_threshold = str2num(handles.less_than_threshold.String);
handles.threshold_name = [handles.metric_name,handles.greater_than_threshold.String,'_<_X','_<_',handles.less_than_threshold.String];
handles.threshold_passed = (handles.metric_value' > lower_threshold) & (handles.metric_value' < upper_threshold);
guidata(hObject, handles);
uiresume(handles.figure1);

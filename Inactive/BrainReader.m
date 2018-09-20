function varargout = BrainReader(varargin)
% BRAINREADER MATLAB code for BrainReader.fig
%      BRAINREADER, by itself, creates a new BRAINREADER or raises the existing
%      singleton*.
%
%      H = BRAINREADER returns the handle to a new BRAINREADER or the handle to
%      the existing singleton*.
%
%      BRAINREADER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BRAINREADER.M with the given input arguments.
%
%      BRAINREADER('Property','Value',...) creates a new BRAINREADER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrainReader_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrainReader_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrainReader

% Last Modified by GUIDE v2.5 12-Apr-2018 16:24:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BrainReader_OpeningFcn, ...
                   'gui_OutputFcn',  @BrainReader_OutputFcn, ...
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


% --- Executes just before BrainReader is made visible.
function BrainReader_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BrainReader (see VARARGIN)



% Choose default command line output for BrainReader
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BrainReader wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BrainReader_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[name,path]=uigetfile('*.mat','Select To Open');
if name==0
return
end
filename=[path name];
handles.CurrentData=load(filename);
PlotCellMap(handles.CurrentData.spPos, 0, handles.CurrentData.fluorescence_time_series, [5,30], 0.2);
assignin('base','currentdata',handles.CurrentData.spPos);
assignin('base','hObject',hObject);
guidata(hObject, handles);




function Transparency_For_Dots_Callback(hObject, eventdata, handles)
% hObject    handle to Transparency_For_Dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Transparency_For_Dots as text
%        str2double(get(hObject,'String')) returns contents of Transparency_For_Dots as a double
handles.TransForDots=str2double(get(hObject,'String'))
assignin('base','hObject2',hObject);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Transparency_For_Dots_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Transparency_For_Dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RefreshPlot.
function RefreshPlot_Callback(hObject, eventdata, handles)
% hObject    handle to RefreshPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PlotCellMap(handles.CurrentData.spPos, 0, handles.CurrentData.fluorescence_time_series, [handles.LowerLimitOfDotSize,handles.UpperLimitOfDotSize], handles.TransForDots);
assignin('base','hObject5',hObject);
guidata(hObject, handles);


function UpperLimitOfDotSize_Callback(hObject, eventdata, handles)
% hObject    handle to UpperLimitOfDotSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UpperLimitOfDotSize as text
%        str2double(get(hObject,'String')) returns contents of UpperLimitOfDotSize as a double
handles.UpperLimitOfDotSize=str2num(get(hObject,'String'))
assignin('base','handle',handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function UpperLimitOfDotSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UpperLimitOfDotSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LowerLimitOfDotSize_Callback(hObject, eventdata, handles)
% hObject    handle to LowerLimitOfDotSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowerLimitOfDotSize as text
%        str2double(get(hObject,'String')) returns contents of LowerLimitOfDotSize as a double
handles.LowerLimitOfDotSize=str2num(get(hObject,'String'))
assignin('base','hObject4',hObject);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LowerLimitOfDotSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowerLimitOfDotSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

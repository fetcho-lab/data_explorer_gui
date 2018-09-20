function varargout = FindCorrCellsOfSeries(varargin)
% FINDCORRCELLSOFSERIES MATLAB code for FindCorrCellsOfSeries.fig
%      FINDCORRCELLSOFSERIES, by itself, creates a new FINDCORRCELLSOFSERIES or raises the existing
%      singleton*.
%
%      H = FINDCORRCELLSOFSERIES returns the handle to a new FINDCORRCELLSOFSERIES or the handle to
%      the existing singleton*.
%
%      FINDCORRCELLSOFSERIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINDCORRCELLSOFSERIES.M with the given input arguments.
%
%      FINDCORRCELLSOFSERIES('Property','Value',...) creates a new FINDCORRCELLSOFSERIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FindCorrCellsOfSeries_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FindCorrCellsOfSeries_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FindCorrCellsOfSeries

% Last Modified by GUIDE v2.5 15-Apr-2018 22:18:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FindCorrCellsOfSeries_OpeningFcn, ...
                   'gui_OutputFcn',  @FindCorrCellsOfSeries_OutputFcn, ...
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


% --- Executes just before FindCorrCellsOfSeries is made visible.
function FindCorrCellsOfSeries_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FindCorrCellsOfSeries (see VARARGIN)

% Choose default command line output for FindCorrCellsOfSeries
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FindCorrCellsOfSeries wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FindCorrCellsOfSeries_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function k_value_Callback(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_value as text
%        str2double(get(hObject,'String')) returns contents of k_value as a double


% --- Executes during object creation, after setting all properties.
function k_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Refresh_Plot.
function Refresh_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadSourceData.
function LoadSourceData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSourceData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadSeries.
function LoadSeries_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportResultList.
function ExportResultList_Callback(hObject, eventdata, handles)
% hObject    handle to ExportResultList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

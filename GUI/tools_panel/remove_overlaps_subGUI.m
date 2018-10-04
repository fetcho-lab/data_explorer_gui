function varargout = remove_overlaps_subGUI(varargin)
% REMOVE_OVERLAPS_SUBGUI MATLAB code for remove_overlaps_subGUI.fig
%      REMOVE_OVERLAPS_SUBGUI by itself, creates a new REMOVE_OVERLAPS_SUBGUI or raises the
%      existing singleton*.
%
%      H = REMOVE_OVERLAPS_SUBGUI returns the handle to a new REMOVE_OVERLAPS_SUBGUI or the handle to
%      the existing singleton*.
%
%      REMOVE_OVERLAPS_SUBGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REMOVE_OVERLAPS_SUBGUI.M with the given input arguments.
%
%      REMOVE_OVERLAPS_SUBGUI('Property','Value',...) creates a new REMOVE_OVERLAPS_SUBGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before remove_overlaps_subGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to remove_overlaps_subGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help remove_overlaps_subGUI

% Last Modified by GUIDE v2.5 03-Oct-2018 09:30:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @remove_overlaps_subGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @remove_overlaps_subGUI_OutputFcn, ...
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

% --- Executes just before remove_overlaps_subGUI is made visible.
function remove_overlaps_subGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to remove_overlaps_subGUI (see VARARGIN)

% Choose default command line output for remove_overlaps_subGUI
handles.output = 'Yes';

handles.Intensity = varargin{1};
handles.Pos = varargin{2};
handles.roi = struct;

% Update handles structure
guidata(hObject, handles);

% % Insert custom Title and Text if specified by the user
% % Hint: when choosing keywords, be sure they are not easily confused 
% % with existing figure properties.  See the output of set(figure) for
% % a list of figure properties.
% if(nargin > 3)
%     for index = 1:2:(nargin-3)
%         if nargin-3==index, break, end
%         switch lower(varargin{index})
%          case 'title'
%           set(hObject, 'Name', varargin{index+1});
%          case 'string'
%           set(handles.text1, 'String', varargin{index+1});
%         end
%     end
% end

% % Determine the position of the dialog - centered on the callback figure
% % if available, else, centered on the screen
% FigPos=get(0,'DefaultFigurePosition');
% OldUnits = get(hObject, 'Units');
% set(hObject, 'Units', 'pixels');
% OldPos = get(hObject,'Position');
% FigWidth = OldPos(3);
% FigHeight = OldPos(4);
% if isempty(gcbf)
%     ScreenUnits=get(0,'Units');
%     set(0,'Units','pixels');
%     ScreenSize=get(0,'ScreenSize');
%     set(0,'Units',ScreenUnits);
% 
%     FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
%     FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
% else
%     GCBFOldUnits = get(gcbf,'Units');
%     set(gcbf,'Units','pixels');
%     GCBFPos = get(gcbf,'Position');
%     set(gcbf,'Units',GCBFOldUnits);
%     FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
%                    (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
% end
% FigPos(3:4)=[FigWidth FigHeight];
% set(hObject, 'Position', FigPos);
% set(hObject, 'Units', OldUnits);
% 
% % Show a question icon from dialogicons.mat - variables questIconData
% % and questIconMap
% % load dialogicons.mat
% 
% % IconData=questIconData;
% % questIconMap(256,:) = get(handles.figure1, 'Color');
% % IconCMap=questIconMap;
% 
% % Img=image(IconData, 'Parent', handles.axes1);
% % set(handles.figure1, 'Colormap', IconCMap);
% % 
% % set(handles.axes1, ...
% %     'Visible', 'off', ...
% %     'YDir'   , 'reverse'       , ...
% %     'XLim'   , get(Img,'XData'), ...
% %     'YLim'   , get(Img,'YData')  ...
% %     );
% 
% % Make the GUI modal
% set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes remove_overlaps_subGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = remove_overlaps_subGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% handles = guidata(hObject);
varargout{1} = handles.roi;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Update handles structure
handles = guidata(hObject);
handles.figure1.Visible = 'off';
drawnow;
threshold = str2num( handles.um_threshold.String );
adjCheck = handles.checkbox1.Value;

[toKeep, toRemove ]= remove_overlaps_fts(handles.Intensity, handles.Pos, threshold, adjCheck);
newRoi(1).name = 'toKeep';
newRoi(1).members = toKeep;
newRoi(2).name = 'toRemove';
newRoi(2).members = toRemove;

handles.roi = newRoi;

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    



function um_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to um_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of um_threshold as text
%        str2double(get(hObject,'String')) returns contents of um_threshold as a double


% --- Executes during object creation, after setting all properties.
function um_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to um_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

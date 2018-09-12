function varargout = TopCorrViewer(varargin)
% TOPCORRVIEWER MATLAB code for TopCorrViewer.fig
%      TOPCORRVIEWER, by itself, creates a new TOPCORRVIEWER or raises the existing
%      singleton*.
%
%      H = TOPCORRVIEWER returns the handle to a new TOPCORRVIEWER or the handle to
%      the existing singleton*.
%
%      TOPCORRVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOPCORRVIEWER.M with the given input arguments.
%
%      TOPCORRVIEWER('Property','Value',...) creates a new TOPCORRVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TopCorrViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TopCorrViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TopCorrViewer

% Last Modified by GUIDE v2.5 27-Apr-2018 17:57:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TopCorrViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @TopCorrViewer_OutputFcn, ...
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


% --- Executes just before TopCorrViewer is made visible.
function TopCorrViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TopCorrViewer (see VARARGIN)

% Choose default command line output for TopCorrViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TopCorrViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function handles=loadDataset(filepath,handles,eventdata)
%switches data set
load(filepath);
handles.fts = fluorescence_time_series;
handles.spPos = spPos;
global AllPosition %Defined a global here!
AllPosition=spPos;
handles.Sc = Sc;
handles.SizeOfDots=20;%Set default dot size here

assignin('base', 'handles0',handles);
handles = plot_pos_maps(handles);

function handles=plot_pos_maps(handles)
axes(handles.PosMap), cla
assignin('base','handles1',handles);
posAllCell=handles.spPos;
handles.ColorOfDots=max(handles.fts');
PosMap=scatter3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),handles.SizeOfDots,handles.ColorOfDots,'.');%,'PickableParts','none'
hold on
xlabel(handles.PosMap,'X');
ylabel(handles.PosMap,'Y');
zlabel(handles.PosMap,'Z');
set(gca,'Xcolor','w');
set(gca,'Ycolor','w');
set(gca,'Zcolor','w');
set(gca,'Gridcolor','k');
% PosMap.XColor='w';
% set(gca,'buttondownfcn',@clicky);
% assignin('base','x',posAllCell(:,1));
% plot3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),'.','color',[0.2 0.2 0.2],'hittest','off');
set(handles.PosMap,'Visible','on');
colormap jet;
% caxis(lim1);
CB1=colorbar;
CB1.Color='w';
CB1.Label.String='Fluorescence Intensity of Cells (Dots)';

grid on


% --- Outputs from this function are returned to the command line.
function varargout = TopCorrViewer_OutputFcn(hObject, eventdata, handles) 
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
[currentfile,path] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(path);
handles.currentfile = [path,currentfile];
handles = loadDataset(handles.currentfile,handles,eventdata);
handles.CurrentCorrCalculated=0;
guidata(hObject,handles);


% --- Executes on button press in PlotTopCorr.
function PlotTopCorr_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTopCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'CurrentCorrCalculated')==0||handles.CurrentCorrCalculated==0;
        [handles.AllMaxcor,handles.AllMaxcorno,handles.Allcor]=TopCorrCellInGroup_corrcoef(handles.fts);
        handles.CurrentCorrCalculated=1;
end
assignin('base','AllMaxcor',handles.AllMaxcor);
assignin('base','AllMaxcorno',handles.AllMaxcorno);
assignin('base','Allcor',handles.Allcor);

axes(handles.PosMap);
PlotTopCorr(handles.AllMaxcor,handles.AllMaxcorno,handles.spPos,handles.fts,handles);
guidata(hObject,handles);

function PlotTopCorr(AllMaxcor, AllMaxcorno, AllPos, AllFluo,handles);
handles.Dis_thre=[handles.DisThreLowerLimit,handles.DisThreUpperLimit];
handles.Cor_thre=handles.CorrThre;
handles.WidthRangeForLines=[handles.LineWidthLowerLimit,handles.LineWidthUpperLimit];
handles.RandomLineColor=handles.RandomLineColorCheck;
%------------Adjustable--------------
% handles.Dis_thre=[0,999];
% handles.Cor_thre=0.995;
% handles.WidthRangeForLines=[2,8];
% handles.RandomLineColor=1;
handles.FluoLineColorMaxVal=max(max(AllFluo));%set FluoLineColorMaxVal here!!!!!!
%---------------------------------------

%---------The code following seems should have the same plotting function
%but the result is werid---------------
DispList1=find(AllMaxcor>handles.Cor_thre);
DispList2=AllMaxcorno(DispList1);
% for i=1:length(DispList1)
%     tempx=[AllPos(DispList1,1),AllPos(DispList2,1)];
%     tempy=[AllPos(DispList1,2),AllPos(DispList2,2)];
%     tempz=[AllPos(DispList1,3),AllPos(DispList2,3)];
%     linewidth=(((AllMaxcor(DispList1(i))-handles.Cor_thre)/(1-handles.Cor_thre))*(handles.WidthRangeForLines(2)-handles.WidthRangeForLines(1)))+handles.WidthRangeForLines(1);
%     Corrplot=plot3(tempx, tempy, tempz, 'k-','LineWidth',linewidth);
% end
%------------------------------------
Plotted_Correlation_No=0;
        for i3=1:size(AllFluo,1)
            tempx=[AllPos(i3,1),AllPos(AllMaxcorno(i3),1)];
            tempy=[AllPos(i3,2),AllPos(AllMaxcorno(i3),2)];
            tempz=[AllPos(i3,3),AllPos(AllMaxcorno(i3),3)];
            dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
            if dist>handles.Dis_thre(1)
                if dist<handles.Dis_thre(2)
                    if AllMaxcor(i3)>handles.Cor_thre
                        if AllMaxcor(i3)<=1
                        linewidth=(((AllMaxcor(i3)-handles.Cor_thre)/(1-handles.Cor_thre))*(handles.WidthRangeForLines(2)-handles.WidthRangeForLines(1)))+handles.WidthRangeForLines(1);
                            if handles.RandomLineColor==1
                                Corrplot=plot3(tempx, tempy, tempz, '-','LineWidth',linewidth);
                            else
                                colorR=max(AllFluo(i3,:))/handles.FluoLineColorMaxVal;
                                if colorR>1
                                    colorR=1;
                                end
                                colorG=max(AllFluo(AllMaxcorno(i3),:))/handles.FluoLineColorMaxVal;
                                if colorG>1
                                    colorG=1;
                                end
                                colorB=0;
                                Corrplot=plot3(tempx, tempy, tempz, '-','Color',[colorR, colorG, colorB],'LineWidth',linewidth);
                                hold on
                                Plotted_Correlation_No=Plotted_Correlation_No+1
                            end
                        end
                    end
                end
            end
        end 
% guidata(hObject,handles);

function SizeOfDots_Callback(hObject, eventdata, handles)
% hObject    handle to SizeOfDots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeOfDots as text
%        str2double(get(hObject,'String')) returns contents of SizeOfDots as a double
handles.SizeOfDots=str2double(get(hObject,'String'))
handles=plot_pos_maps(handles);
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


% --- Executes on button press in RefreshCellPlot.
function RefreshCellPlot_Callback(hObject, eventdata, handles)
% hObject    handle to RefreshCellPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=plot_pos_maps(handles);


% --- Executes on mouse press over axes background.


% hObject    handle to PosMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
% assignin('base','eventdata2',eventdata);
% assignin('base','handles',handles);
% 
% position_list = single( handles.spPos );
% assignin('base','poslist',position_list);
% % clickPosition = eventdata.IntersectionPoint;
% clickpoint=ginput(1)
% clickPosition=clickpoint(1,:);
% distance_to_list = pdist2(single( clickPosition ),position_list);
% 
% [mD,cellSelected] = min(distance_to_list);
% handles.cellSelector.Value = cellSelected;
% assignin('base','handles',handles);
% assignin('base','cellselected',cellSelected);
% guidata(hObject,handles);


% --- Executes on button press in PlotTopCorrInROIs.
function PlotTopCorrInROIs_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTopCorrInROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.RoiMaxCorVal,handles.RoiMaxCorCellNo,handles.AllRoiCor]=TopCorrCellInGroup_corrcoef(handles.FluoOfRoi);
axes(handles.PosMap);
PlotTopCorr(handles.RoiMaxCorVal,handles.RoiMaxCorCellNo,handles.PosOfRoi,handles.FluoOfRoi,handles);

guidata(hObject,handles);


% --- Executes on button press in ExportCorrResult.
function ExportCorrResult_Callback(hObject, eventdata, handles)
% hObject    handle to ExportCorrResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MaxCorrelationValOfAllSelectedCells=handles.AllMaxcor;
MaxCorrelatedCellNoOfAllSelectedCells=handles.AllMaxcorno;
EntireCorrelationMatrix=handles.Allcor;
[f1,path1] = uiputfile('CorrCalculationResult.mat');
save([path1,f1],'MaxCorrelationValOfAllSelectedCells','MaxCorrelatedCellNoOfAllSelectedCells','EntireCorrelationMatrix');

guidata(hObject,handles);

% --- Executes on button press in LoadCorrResult.
function LoadCorrResult_Callback(hObject, eventdata, handles)
% hObject    handle to LoadCorrResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Corrcurrentfile,Corrpath] = uigetfile('.mat','Please select CorrCalculationResult.mat file');
cd(Corrpath);
handles.currentCorrfile = [Corrpath,Corrcurrentfile];
load(handles.currentCorrfile);

handles.AllMaxcor=MaxCorrelationValOfAllSelectedCells;
handles.AllMaxcorno=MaxCorrelatedCellNoOfAllSelectedCells;
handles.Allcor=EntireCorrelationMatrix;

handles.CurrentCorrCalculated=1;
guidata(hObject,handles);


% --- Executes on button press in ExportROICorrResult.
function ExportROICorrResult_Callback(hObject, eventdata, handles)
% hObject    handle to ExportROICorrResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.RoiMaxCorVal,handles.RoiMaxCorCellNo,handles.AllRoiCor
MaxCorrelationValOfAllSelectedCells=handles.RoiMaxCorVal;
MaxCorrelatedCellNoOfAllSelectedCells=handles.RoiMaxCorCellNo;
EntireCorrelationMatrix=handles.AllRoiCor;
ROICellNo=handles.CurrentRoi;
[f2,path2] = uiputfile('CorrCalculationResult_ROI.mat');
save([path2,f2],'MaxCorrelationValOfAllSelectedCells','MaxCorrelatedCellNoOfAllSelectedCells','EntireCorrelationMatrix','ROICellNo');

guidata(hObject,handles);

function FTS1NoInput_Callback(hObject, eventdata, handles)
% hObject    handle to FTS1NoInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FTS1NoInput as text
%        str2double(get(hObject,'String')) returns contents of FTS1NoInput as a double
handles.CellNo1=str2double(get(hObject,'String'));
axes(handles.FTS1);
PlotFTS(handles.fts((handles.CellNo1),:));

set(handles.TopCorrCoef1, 'String', num2str(handles.AllMaxcor(handles.CellNo1)));
set(handles.TopCorrCellNo1, 'String', num2str(handles.AllMaxcorno(handles.CellNo1)));

guidata(hObject, handles);




function FTS2NoInput_Callback(hObject, eventdata, handles)
% hObject    handle to FTS2NoInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FTS2NoInput as text
%        str2double(get(hObject,'String')) returns contents of FTS2NoInput as a double
handles.CellNo2=str2double(get(hObject,'String'));
axes(handles.FTS2);
PlotFTS(handles.fts((handles.CellNo2),:));

set(handles.TopCorrCoef2, 'String', num2str(handles.AllMaxcor(handles.CellNo2)));
set(handles.TopCorrCellNo2, 'String', num2str(handles.AllMaxcorno(handles.CellNo2)));

guidata(hObject, handles);

function PlotFTS(fluoseries)
tempx=[1:length(fluoseries)];
plot(tempx, fluoseries,'g');

% --- Executes during object creation, after setting all properties.
function FTS2NoInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FTS2NoInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PlotCorrCell1.
function PlotCorrCell1_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCorrCell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CellNo1=handles.AllMaxcorno(handles.CellNo2);
set(handles.FTS1NoInput, 'String', num2str(handles.CellNo1));
axes(handles.FTS1);
PlotFTS(handles.fts((handles.CellNo1),:));

set(handles.TopCorrCoef1, 'String', num2str(handles.AllMaxcor(handles.CellNo1)));
set(handles.TopCorrCellNo1, 'String', num2str(handles.AllMaxcorno(handles.CellNo1)));
guidata(hObject, handles);


% --- Executes on button press in PlotCorCell2.
function PlotCorCell2_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCorCell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CellNo2=handles.AllMaxcorno(handles.CellNo1);
set(handles.FTS2NoInput, 'String', num2str(handles.CellNo2));
axes(handles.FTS2);
PlotFTS(handles.fts((handles.CellNo2),:));

set(handles.TopCorrCoef2, 'String', num2str(handles.AllMaxcor(handles.CellNo2)));
set(handles.TopCorrCellNo2, 'String', num2str(handles.AllMaxcorno(handles.CellNo2)));
guidata(hObject, handles);


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
handles.CirclePoint=plot3(handles.spPos(handles.CellNoToFind,1), handles.spPos(handles.CellNoToFind,2),handles.spPos(handles.CellNoToFind,3),'ro','linewidth',2);

set(handles.XVal, 'String', num2str(handles.spPos(handles.CellNoToFind,1)));
set(handles.YVal, 'String', num2str(handles.spPos(handles.CellNoToFind,2)));
set(handles.ZVal, 'String', num2str(handles.spPos(handles.CellNoToFind,3)));
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function CellNoToFind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CellNoToFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Find.
function Find_Callback(hObject, eventdata, handles)
% hObject    handle to Find (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function FTS1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FTS1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate FTS1


% --- Executes on button press in HideCorrPlot.
function HideCorrPlot_Callback(hObject, eventdata, handles)
% hObject    handle to HideCorrPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  



function UpperLimitDisThreInput_Callback(hObject, eventdata, handles)
% hObject    handle to UpperLimitDisThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UpperLimitDisThreInput as text
%        str2double(get(hObject,'String')) returns contents of UpperLimitDisThreInput as a double

handles.DisThreUpperLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function UpperLimitDisThreInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UpperLimitDisThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.DisThreUpperLimit=999;
guidata(hObject,handles);


function LowerLimitDisThreInput_Callback(hObject, eventdata, handles)
% hObject    handle to LowerLimitDisThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowerLimitDisThreInput as text
%        str2double(get(hObject,'String')) returns contents of LowerLimitDisThreInput as a double
handles.DisThreLowerLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function LowerLimitDisThreInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowerLimitDisThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.DisThreLowerLimit=0;
guidata(hObject,handles);


function CorrThreInput_Callback(hObject, eventdata, handles)
% hObject    handle to CorrThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CorrThreInput as text
%        str2double(get(hObject,'String')) returns contents of CorrThreInput as a double
handles.CorrThre=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function CorrThreInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.CorrThre=0.995;
guidata(hObject,handles);

% --- Executes on button press in RandomColorForLinesCheck.
function RandomColorForLinesCheck_Callback(hObject, eventdata, handles)
% hObject    handle to RandomColorForLinesCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RandomColorForLinesCheck
handles.RandomLineColorCheck=get(hObject,'Value');
guidata(hObject,handles);


function UpperLineWidthInput_Callback(hObject, eventdata, handles)
% hObject    handle to UpperLineWidthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UpperLineWidthInput as text
%        str2double(get(hObject,'String')) returns contents of UpperLineWidthInput as a double

handles.LineWidthUpperLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function UpperLineWidthInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UpperLineWidthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.LineWidthUpperLimit=8;
guidata(hObject,handles);


function LowerLineWidthInput_Callback(hObject, eventdata, handles)
% hObject    handle to LowerLineWidthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowerLineWidthInput as text
%        str2double(get(hObject,'String')) returns contents of LowerLineWidthInput as a double
handles.LineWidthLowerLimit=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function LowerLineWidthInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowerLineWidthInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.LineWidthLowerLimit=2;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function RandomColorForLinesCheck_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RandomColorForLinesCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.RandomLineColorCheck=0;
guidata(hObject,handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[roifile,roipath] = uigetfile('.mat','Please select roi .mat file');
cd(roipath);
load([roipath,roifile]);
handles.CurrentRoi = roi;
handles.FluoOfRoi=handles.fts(handles.CurrentRoi,:);
handles.PosOfRoi=handles.spPos(handles.CurrentRoi,:);

guidata(hObject,handles);



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

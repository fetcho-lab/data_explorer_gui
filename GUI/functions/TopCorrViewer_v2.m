function varargout = TopCorrViewer_v2(varargin)
% TOPCORRVIEWER_V2 MATLAB code for TopCorrViewer_v2.fig
%      TOPCORRVIEWER_V2, by itself, creates a new TOPCORRVIEWER_V2 or raises the existing
%      singleton*.
%
%      H = TOPCORRVIEWER_V2 returns the handle to a new TOPCORRVIEWER_V2 or the handle to
%      the existing singleton*.
%
%      TOPCORRVIEWER_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOPCORRVIEWER_V2.M with the given input arguments.
%
%      TOPCORRVIEWER_V2('Property','Value',...) creates a new TOPCORRVIEWER_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TopCorrViewer_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TopCorrViewer_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TopCorrViewer_v2

% Last Modified by GUIDE v2.5 26-Jun-2018 17:31:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TopCorrViewer_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @TopCorrViewer_v2_OutputFcn, ...
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


% --- Executes just before TopCorrViewer_v2 is made visible.
function TopCorrViewer_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TopCorrViewer_v2 (see VARARGIN)

% Choose default command line output for TopCorrViewer_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TopCorrViewer_v2 wait for user response (see UIRESUME)
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
assignin('base','AllPos',handles.spPos);
global AllPosition %Defined a global here!
AllPosition=spPos;
handles.Sc = Sc;
handles.SizeOfDots=20;%Set default dot size here
handles.CurrentCorrType=1;
%---------Load dFF if exist--------
if exist('dFF','var')
    handles.dFF = dFF;
else
    handles.dFF = NaN;
end
%---------------------------------
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
CB1.Label.String='Max Fluorescence Intensity of Cells (Dots)';
grid on

% function handles=plot_pos_maps(handles)
% axes(handles.PosMap), cla
%  
% posAllCell=handles.spPos;
% handles.ColorOfDots=max(handles.fts');
% PosMap=scatter3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),handles.SizeOfDots,handles.ColorOfDots,'.');%,'PickableParts','none'
% hold on
% xlabel(handles.PosMap,'X');
% ylabel(handles.PosMap,'Y');
% zlabel(handles.PosMap,'Z');
% set(gca,'Xcolor','w');
% set(gca,'Ycolor','w');
% set(gca,'Zcolor','w');
% set(gca,'Gridcolor','k');
% % PosMap.XColor='w';
% % set(gca,'buttondownfcn',@clicky);
% % assignin('base','x',posAllCell(:,1));
% % plot3(posAllCell(:,1),posAllCell(:,2),posAllCell(:,3),'.','color',[0.2 0.2 0.2],'hittest','off');
% set(handles.PosMap,'Visible','on');
% colormap jet;
% % caxis(lim1);
% CB1=colorbar;
% CB1.Color='w';
% CB1.Label.String='Fluorescence Intensity of Cells (Dots)';
% 
% grid on


% --- Outputs from this function are returned to the command line.
function varargout = TopCorrViewer_v2_OutputFcn(hObject, eventdata, handles) 
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
[handles.currentfilename,path] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(path);
handles.currentfile = [path,handles.currentfilename];
handles = loadDataset(handles.currentfile,handles,eventdata);
handles.CurrentCorrCalculated=0;

set(handles.FileNameDisplay,'String',handles.currentfile);
guidata(hObject,handles);


% --- Executes on button press in PlotTopCorr.
function PlotTopCorr_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTopCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.PosMap), cla
plot_pos_maps(handles.spPos,handles.fts,handles.SizeOfDots);
hold on

% if isfield(handles,'CurrentCorrCalculated')==0||handles.CurrentCorrCalculated==0
    if handles.CurrentCorrCalculated==0
        [handles.AllMaxcor,handles.AllMaxcorno,handles.Allcor]=TopCorrCellInGroup_corrcoef(handles.fts);
        handles.CurrentCorrCalculated=1;
end
% assignin('base','AllMaxcor',handles.AllMaxcor);
% assignin('base','AllMaxcorno',handles.AllMaxcorno);
assignin('base','Allcor',handles.Allcor);

axes(handles.PosMap);
handles=PlotTopCorr(handles.Allcor,handles.spPos,handles.fts,handles);
handles=PlotDistanceProfile(handles);
guidata(hObject,handles);


function handles=PlotDistanceProfile(handles)
%------------Produce Distance Profile Heatmap----------
figure;
% DistCat=[0:20:700];%Adjust Distance Catagory Here!!!
MaxDist=max(max(handles.CorrDist))
DistCat=[0:5:MaxDist];
CorrNo=zeros(size(DistCat));

for i=1:length(DistCat)-1
    CorrNo(i)=length(find(handles.CorrDist>DistCat(i)&handles.CorrDist<DistCat(i+1)))/2;
%     assignin('base','corrno',CorrNo)
end
bar(DistCat,CorrNo);
title(strcat('Length Distribution of Correlations with coef>' ,num2str(handles.LowerCorrThre),'and <',num2str(handles.HigherCorrThre)));
xlabel('Distance Catagories');
ylabel('Correlation Pair No');
% guidata(hObject,handles);

function handles=PlotTopCorr(Allcor, AllPos, AllFluo,handles)
handles.Dis_thre=[handles.DisThreLowerLimit,handles.DisThreUpperLimit];
handles.Cor_thre=[handles.LowerCorrThre,handles.HigherCorrThre];
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
[DispList_i, DispList_j]=find(Allcor>=handles.LowerCorrThre&Allcor<=handles.HigherCorrThre);
% DispList2=AllMaxcorno(DispList1i);
PlottedCorrelationNo=0;
handles.CorrDist=zeros(size(handles.Allcor));
for i=1:size(DispList_i,1)
    tempx=[AllPos(DispList_i(i),1),AllPos(DispList_j(i),1)];
    tempy=[AllPos(DispList_i(i),2),AllPos(DispList_j(i),2)];
    tempz=[AllPos(DispList_i(i),3),AllPos(DispList_j(i),3)];
    handles.CorrDist(DispList_i(i),DispList_j(i))=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
    if handles.CorrDist(DispList_i(i),DispList_j(i))>handles.Dis_thre(1)
        if handles.CorrDist(DispList_i(i),DispList_j(i))<handles.Dis_thre(2)
            linewidth=(((Allcor(DispList_i(i),DispList_j(i))-handles.Cor_thre)/(1-handles.Cor_thre))*(handles.WidthRangeForLines(2)-handles.WidthRangeForLines(1)))+handles.WidthRangeForLines(1);

            if handles.RandomLineColor==1
                Corrplot=plot3(tempx, tempy, tempz, '-','LineWidth',linewidth);
                PlottedCorrelationNo=PlottedCorrelationNo+1;
            else
                colorR=max(AllFluo(DispList_i(i),:))/handles.FluoLineColorMaxVal;
                    if colorR>1
                        colorR=1;
                    end
                colorG=max(AllFluo(DispList_j(i),:))/handles.FluoLineColorMaxVal;
                    if colorG>1
                        colorG=1;
                    end
                colorB=0;
                Corrplot=plot3(tempx, tempy, tempz, '-','Color',[colorR, colorG, colorB],'LineWidth',linewidth);
                PlottedCorrelationNo=PlottedCorrelationNo+1;
            end
        end
    end
end
    assignin('base','corrdist2',handles.CorrDist);
PlottedCorrelationNo=PlottedCorrelationNo/2
%------------------------------------
% Plotted_Correlation_No=0;
%         for i3=1:size(AllFluo,1)
%             tempx=[AllPos(i3,1),AllPos(AllMaxcorno(i3),1)];
%             tempy=[AllPos(i3,2),AllPos(AllMaxcorno(i3),2)];
%             tempz=[AllPos(i3,3),AllPos(AllMaxcorno(i3),3)];
%             dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
%             if dist>handles.Dis_thre(1)
%                 if dist<handles.Dis_thre(2)
%                     if AllMaxcor(i3)>handles.Cor_thre
%                         if AllMaxcor(i3)<=1
%                         linewidth=(((AllMaxcor(i3)-handles.Cor_thre)/(1-handles.Cor_thre))*(handles.WidthRangeForLines(2)-handles.WidthRangeForLines(1)))+handles.WidthRangeForLines(1);
%                             if handles.RandomLineColor==1
%                                 Corrplot=plot3(tempx, tempy, tempz, '-','LineWidth',linewidth);
%                             else
%                                 colorR=max(AllFluo(i3,:))/handles.FluoLineColorMaxVal;
%                                 if colorR>1
%                                     colorR=1;
%                                 end
%                                 colorG=max(AllFluo(AllMaxcorno(i3),:))/handles.FluoLineColorMaxVal;
%                                 if colorG>1
%                                     colorG=1;
%                                 end
%                                 colorB=0;
%                                 Corrplot=plot3(tempx, tempy, tempz, '-','Color',[colorR, colorG, colorB],'LineWidth',linewidth);
%                                 hold on
%                                 Plotted_Correlation_No=Plotted_Correlation_No+1
%                             end
%                         end
%                     end
%                 end
%             end
%         end 
% guidata(hObject,handles);

function SizeOfDots_Callback(hObject, eventdata, handles)
% hObject    handle to SizeOfDots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SizeOfDots as text
%        str2double(get(hObject,'String')) returns contents of SizeOfDots as a double
handles.SizeOfDots=str2double(get(hObject,'String'));
% handles=plot_pos_maps(handles);

axes(handles.PosMap), cla
plot_pos_maps(handles.spPos,handles.fts,handles.SizeOfDots);
hold on
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

% handles=plot_pos_maps(handles);

axes(handles.PosMap), cla
plot_pos_maps(handles.spPos,handles.fts,handles.SizeOfDots);
hold on

handles=guidata(hObject);

% --- Executes on button press in PlotTopCorrInROIs.
function PlotTopCorrInROIs_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTopCorrInROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.PosMap);
handles=PlotTopCorr(handles.Allcor,handles.PosOfRoi,handles.FluoOfRoi,handles);

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
% uiresume(handles.figure1);
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
handles=PlotFTS1(handles);
guidata(hObject, handles);

function handles=PlotFTS1(handles)
axes(handles.FTS1);
if get(handles.dFFCheck1,'Value')==1
    PlotFTS(handles.dFF((handles.CellNo1),:));
else
    PlotFTS(handles.fts((handles.CellNo1),:));
end
set(handles.FTS1NoInput,'String',num2str(handles.CellNo1));

[CorrValList1 handles.CorrCellNoList1]=sort(handles.Allcor(:,handles.CellNo1),'descend');
DistList1=DistOfVectorsToVectors(handles.spPos(handles.CellNo1,:), handles.spPos(handles.CorrCellNoList1,:));
%     assignin('base','CorrCellNo',CorrCellNo);
    StringToDisplay=strcat(num2str(handles.CorrCellNoList1),'  -',num2str(CorrValList1),' -',num2str(DistList1));
    set(handles.CorrList1,'String',StringToDisplay);


function [Dist]=DistOfVectorsToVectors(A,B)
% A_{N*d}, B_{M*d},Dist_{N*M}
[N d]=size(A);
[M d]=size(B);
Dist = sqrt(sum(A.*A, 2)*ones(1, M) + ones(N, 1) * sum(B.*B, 2)' - 2*A*B');
Dist=Dist';


function FTS2NoInput_Callback(hObject, eventdata, handles)
% hObject    handle to FTS2NoInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FTS2NoInput as text
%        str2double(get(hObject,'String')) returns contents of FTS2NoInput as a double
handles.CellNo2=str2double(get(hObject,'String'));
handles=PlotFTS2(handles);
guidata(hObject, handles);

function handles=PlotFTS2(handles)
axes(handles.FTS2);
if get(handles.dFFCheck2,'Value')==1
    PlotFTS(handles.dFF((handles.CellNo2),:));
else
    PlotFTS(handles.fts((handles.CellNo2),:));
end
set(handles.FTS2NoInput,'String',num2str(handles.CellNo2));

[CorrValList2 handles.CorrCellNoList2]=sort(handles.Allcor(:,handles.CellNo2),'descend');
DistList2=DistOfVectorsToVectors(handles.spPos(handles.CellNo2,:), handles.spPos(handles.CorrCellNoList2,:));
%     assignin('base','CorrCellNo',CorrCellNo);
    StringToDisplay=strcat(num2str(handles.CorrCellNoList2),'  -',num2str(CorrValList2),' -',num2str(DistList2));
    set(handles.CorrList2,'String',StringToDisplay);
    
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


function CellNoToFind_Callback(hObject, eventdata, handles)
% hObject    handle to CellNoToFind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CellNoToFind as text
%        str2double(get(hObject,'String')) returns contents of CellNoToFind as a double
handles.CellNoToFind=str2double(get(hObject,'String'));

handles=DrawCirclePoint(handles);

set(handles.XVal, 'String', num2str(handles.spPos(handles.CellNoToFind,1)));
set(handles.YVal, 'String', num2str(handles.spPos(handles.CellNoToFind,2)));
set(handles.ZVal, 'String', num2str(handles.spPos(handles.CellNoToFind,3)));
guidata(hObject,handles);



function handles=DrawCirclePoint(handles)

axes(handles.PosMap);
hold on
if isfield(handles,'CirclePoint');
    delete(handles.CirclePoint);
end
% assignin('base','spPos',handles.spPos);
% assignin('base','CellNoToFind',handles.CellNoToFind);
handles.CirclePoint=plot3(handles.spPos(handles.CellNoToFind,1), handles.spPos(handles.CellNoToFind,2),handles.spPos(handles.CellNoToFind,3),'ro','linewidth',3,'Markersize',10);

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


function LowerCorrThreInput_Callback(hObject, eventdata, handles)
% hObject    handle to LowerCorrThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowerCorrThreInput as text
%        str2double(get(hObject,'String')) returns contents of LowerCorrThreInput as a double
handles.LowerCorrThre=str2double(get(hObject,'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function LowerCorrThreInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowerCorrThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.LowerCorrThre=0.995;
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
handles.LineWidthUpperLimit=2;
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
handles.LineWidthLowerLimit=1;
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
if get(hObject,'Value') ==2
    handles.CurrentCorrType=2;
    print1='Corr type set to Spearman'
elseif get(hObject,'Value')==1
    handles.CurrentCorrType=1;
    print1='Corr type set to Pearson'
end
guidata(hObject,handles);

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


% --- Executes on button press in LoadRoIFile.
function LoadRoIFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadRoIFile (see GCBO)
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


% --- Executes on selection change in CorrList2.
function CorrList2_Callback(hObject, eventdata, handles)
% hObject    handle to CorrList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrList2
handles.CellNo1=handles.CorrCellNoList2(get(hObject,'Value'));
handles.CellNoToFind=handles.CorrCellNoList2(get(hObject,'Value'));

handles=DrawCirclePoint(handles);
handles=PlotFTS1(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function CorrList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CorrList1.
function CorrList1_Callback(hObject, eventdata, handles)
% hObject    handle to CorrList1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrList1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrList1
% assignin('base','fffffff',get(hObject,'Value'));
handles.CellNo2=handles.CorrCellNoList1(get(hObject,'Value'));
handles.CellNoToFind=handles.CorrCellNoList1(get(hObject,'Value'));

handles=DrawCirclePoint(handles);
handles=PlotFTS2(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function CorrList1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrList1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowCell1Only.
function ShowCell1Only_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCell1Only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CellNoToPlot=find(handles.Allcor(:,handles.CellNo1)>=handles.LowerCorrThre&handles.Allcor(:,handles.CellNo2)<=handles.HigherCorrThre);
CellNoToPlot=[CellNoToPlot;handles.CellNo1];
axes(handles.PosMap), cla
plot_pos_maps(handles.spPos(CellNoToPlot,:),handles.fts(CellNoToPlot,:),handles.SizeOfDots);
hold on

axes(handles.PosMap);
handles=PlotTopCorr(handles.Allcor(CellNoToPlot,CellNoToPlot),handles.spPos(CellNoToPlot,:),handles.fts(CellNoToPlot,:),handles);
hold on  

handles.CellNoToFind=handles.CellNo1;
handles=DrawCirclePoint(handles);

guidata(hObject,handles);


% --- Executes on button press in ShowCell2Only.
function ShowCell2Only_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCell2Only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CellNoToPlot=find(handles.Allcor(:,handles.CellNo2)>=handles.LowerCorrThre&&handles.Allcor(:,handles.CellNo2)<=handles.HigherCorrThre);
CellNoToPlot=[CellNoToPlot;handles.CellNo2];
axes(handles.PosMap), cla
plot_pos_maps(handles.spPos(CellNoToPlot,:),handles.fts(CellNoToPlot,:),handles.SizeOfDots);
hold on

axes(handles.PosMap);
handles=PlotTopCorr(handles.Allcor(CellNoToPlot,CellNoToPlot),handles.spPos(CellNoToPlot,:),handles.fts(CellNoToPlot,:),handles);
hold on  

handles.CellNoToFind=handles.CellNo2;
handles=DrawCirclePoint(handles);

guidata(hObject,handles);


% --- Executes on button press in CalculateCorr_Global.
function CalculateCorr_Global_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateCorr_Global (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if isfield(handles,'CurrentCorrCalculated')==0||handles.CurrentCorrCalculated==0
%----------Decide will we overwrite old corr result---------------
if handles.CurrentCorrCalculated==1
    ans0=inputdlg('Correlation Already Calculated, do you want to overwrite the old corr result? 1 for yes, 0 for no.')
    if str2num(ans0{1,1})==1
        CalculateOrNot=1;
    else
        CalculateOrNot=0;
    end
else
    CalculateOrNot=1;
end
%----------Calculate---------
if CalculateOrNot==1
    answer1 = inputdlg('What data do you want to calculate correlation with? Type dff for dF/F, fluo for fluorescence.','dFF or Fluo');
    if strcmp(answer1{1,1},'dff')||strcmp(answer1{1,1},'dFF')
        if isnan(handles.dFF)
        warndlg('dFF of this dataset is not calculated yet, please calculate dFF first! Or use Fluo for correlation instead.','No dFF calculated')
        else
            if isfield(handles,'CurrentCorrType')==0||handles.CurrentCorrType==1
                print1='Calculating Pearson Correlation on dFF...'
                [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.dFF);
            elseif handles.CurrentCorrType==2
                print1='Calculating Spearman Correlation on dFF...'
                [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Spearman(handles.dFF);
            end
            handles.CurrentCorrCalculated=1;
            handles.dFFUsed=1;
            msgbox('Calculation Done!');
        end
    elseif strcmp(answer1{1,1},'fluo')||strcmp(answer1{1,1},'Fluo')
assignin('base','input', handles.CurrentCorrType);
        if isfield(handles,'CurrentCorrType')==0||handles.CurrentCorrType==1
            print1='Calculating Pearson Correlation on Fluo...'
            [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.fts);
        elseif handles.CurrentCorrType==2
            print1='Calculating Spearman Correlation on Fluo...'
            [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Spearman(handles.fts);
        end
        handles.CurrentCorrCalculated=1;
        handles.dFFUsed=0;
        msgbox('Calculation Done!');
    end
end
% assignin('base','AllMaxcor',handles.AllMaxcor);
% assignin('base','AllMaxcorno',handles.AllMaxcorno);
assignin('base','Allcor',handles.Allcor);

handles=PlotCorrelationHeatMap(handles);
handles=PlotCorrProfile(handles);
guidata(hObject,handles);

function handles=PlotCorrelationHeatMap(handles)
figure;
CorrelationHeatMap=imagesc(handles.Allcor);
colorbar
if isfield(handles,'CurrentCorrType')==0||handles.CurrentCorrType==1
    if handles.dFFUsed==1
        title(strcat('Correlation Heat Map of',{' '},handles.currentfilename,', ','Pearson corr with dFF data'));
    else
        title(strcat('Correlation Heat Map of',{' '},handles.currentfilename,', ','Pearson corr with Fluo data'));
    end
elseif handles.CurrentCorrType==2
    if handles.dFFUsed==1
        title(strcat('Correlation Heat Map of',{' '},handles.currentfilename,', ','Spearman corr with dFF data'));
    else
        title(strcat('Correlation Heat Map of',{' '},handles.currentfilename,', ','Spearman corr with Fluo data'));
    end
end

% colormap jet
% ylim=[0 size(Allfluo,1)];
% guidata(hObject,handles);

function handles=PlotCorrProfile(handles)
%------------Produce Distance Profile Heatmap----------
figure;
% DistCat=[0:20:700];%Adjust Distance Catagory Here!!!

CorCat=[0:0.05:1];
CorrNo=zeros(size(CorCat));

for i=1:length(CorCat)-1
    CorrNo(i)=length(find(handles.Allcor>CorCat(i)&handles.Allcor<CorCat(i+1)))/2;
%     assignin('base','corrno',CorrNo)
end
bar(CorCat,CorrNo);
if isfield(handles,'CurrentCorrType')==0||handles.CurrentCorrType==1
    if handles.dFFUsed==1
        title(strcat('Correlation Profile of',{' '},handles.currentfilename,', ','Pearson corr with dFF data'));
    else
        title(strcat('Correlation Profile of',{' '},handles.currentfilename,', ','Pearson corr with Fluo data'));
    end
elseif handles.CurrentCorrType==2
    if handles.dFFUsed==1
        title(strcat('Correlation Profile of',{' '},handles.currentfilename,', ','Spearman corr with dFF data'));
    else
        title(strcat('Correlation Profile of',{' '},handles.currentfilename,', ','Spearman corr with Fluo data'));
    end
end
    
xlabel('Corr Coef Catagories');
ylabel('Correlation Pair No');
% guidata(hObject,handles);



function HigherCorrThreInput_Callback(hObject, eventdata, handles)
% hObject    handle to HigherCorrThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HigherCorrThreInput as text
%        str2double(get(hObject,'String')) returns contents of HigherCorrThreInput as a double
handles.HigherCorrThre=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function HigherCorrThreInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HigherCorrThreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.HigherCorrThre=1;
guidata(hObject,handles);


% --- Executes on button press in dFFCheck1.
function dFFCheck1_Callback(hObject, eventdata, handles)
% hObject    handle to dFFCheck1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dFFCheck1
handles=PlotFTS1(handles);


% --- Executes on button press in dFFCheck2.
function dFFCheck2_Callback(hObject, eventdata, handles)
% hObject    handle to dFFCheck2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dFFCheck2
handles=PlotFTS2(handles);



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CalculateCorr_Regional.
function CalculateCorr_Regional_Callback(hObject, eventdata, handles)
% hObject    handle to CalculateCorr_Regional (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[roifile,roipath] = uigetfile('.mat','Please select roi .mat file');
cd(roipath);
load([roipath,roifile]);
handles.CurrentRoi = find(roi.members==1);
handles.FluoOfRoi=handles.fts(handles.CurrentRoi,:);
handles.dFFOfRoi=handles.dFF(handles.CurrentRoi,:);
handles.PosOfRoi=handles.spPos(handles.CurrentRoi,:);

%----------Decide will we overwrite old corr result---------------
if handles.CurrentCorrCalculated==1
    ans0=inputdlg('Correlation Already Calculated, do you want to overwrite the old corr result? 1 for yes, 0 for no.')
    if str2num(ans0{1,1})==1
        CalculateOrNot=1;
    else
        CalculateOrNot=0;
    end
else
    CalculateOrNot=1;
end
%----------Calculate---------
if CalculateOrNot==1
    answer1 = inputdlg('What data do you want to calculate correlation with? Type dff for dF/F, fluo for fluorescence.','dFF or Fluo');
    if strcmp(answer1{1,1},'dff')||strcmp(answer1{1,1},'dFF')
        if isnan(handles.dFF)
        warndlg('dFF of this dataset is not calculated yet, please calculate dFF first! Or use Fluo for correlation instead.','No dFF calculated')
        else
            if isfield(handles,'CurrentCorrType')==0||handles.CurrentCorrType==1
                print1='Calculating Pearson Correlation on dFF...'
                [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.dFFOfRoi);
            elseif handles.CurrentCorrType==2
                print1='Calculating Spearman Correlation on dFF...'
                [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Spearman(handles.dFFOfRoi);
            end
            handles.CurrentCorrCalculated=1;
            handles.dFFUsed=1;
            msgbox('Regional Calculation Done!');
        end
    elseif strcmp(answer1{1,1},'fluo')||strcmp(answer1{1,1},'Fluo')
        if isfield(handles,'CurrentCorrType')==0||handles.CurrentCorrType==1
            print1='Calculating Pearson Correlation on Fluo...'
            [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Pearson(handles.FluoOfRoi);
        elseif handles.CurrentCorrType==2
            print1='Calculating Spearman Correlation on Fluo...'
            [handles.Allcor, handles.AllPVal]=TopCorrCellInGroup_Spearman(handles.FluoOfRoi);
        end
        handles.CurrentCorrCalculated=1;
        handles.dFFUsed=0;
        msgbox('Regional Calculation Done!');
    end
end
assignin('base','Allcor_Regional',handles.Allcor);

handles=PlotCorrelationHeatMap(handles);
handles=PlotCorrProfile(handles);

guidata(hObject,handles);

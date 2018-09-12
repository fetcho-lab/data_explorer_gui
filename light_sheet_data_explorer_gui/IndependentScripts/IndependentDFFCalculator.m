clear all

answer=inputdlg({'File Directory','File Name','Basal F Ratio'},'Input Data'); 

% if Allfluo is a nxm matrix, n should be the number of cells, m should be
% number of timepoints.
% BasalFVolume should be the number of data points (among all timepoints'
% data) you want to be taken as basal F.

FileToRead=strcat(answer{1,1},'\',answer{2,1});
BasalFRatio=str2num(answer{3,1});
load(FileToRead);

TimepointNo=size(fluorescence_time_series,2);
BasalFVolume=round(BasalFRatio*TimepointNo);
sortedfluo=sort(fluorescence_time_series,2,'ascend');
basalF= mean(sortedfluo(:,1:BasalFVolume),2);
ZeroBasalFCells_discarded=find(basalF==0);
if size(ZeroBasalFCells_discarded)>0
    warndlg(strcat(num2str(size(ZeroBasalFCells_discarded)),' basal F values contains 0, try to make Basal Volume higher? Cells whose basal F is 0 will be deleted.'));
end
fluorescence_time_series(ZeroBasalFCells_discarded,:)=[];
basalF(ZeroBasalFCells_discarded,:)=[];
spPos(ZeroBasalFCells_discarded,:)=[];
spRadiiXYZ(ZeroBasalFCells_discarded,:)=[];
cellSegmentation(:,ZeroBasalFCells_discarded)=[];

TimepointNo2=size(fluorescence_time_series,2);
basalFexpand=basalF*ones(1,TimepointNo2);
AlldFF=(fluorescence_time_series./basalFexpand)-1;


%------delete Cells with too large dFF here----------

dFFThre=10;
[DisRow DisCol]=find(max(AlldFF,[],2)>dFFThre);

AlldFF(DisRow,:)=[];
fluorescence_time_series(DisRow,:)=[];
basalF(DisRow,:)=[];
spPos(DisRow,:)=[];
spRadiiXYZ(DisRow,:)=[];
cellSegmentation(:,DisRow)=[];
%-----------------------------------------

fluorescence_time_series=AlldFF;

[f1,path1] = uiputfile(strcat(answer{3,1},'TrimeddFFOf_',answer{2,1}));
save([path1,f1],'fluorescence_time_series','basalF','BasalFVolume','dFFThre','spPos','spRadiiXYZ','cellSegmentation','Sc');

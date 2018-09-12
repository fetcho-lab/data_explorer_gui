%Cluster based on distance

GroupNo=50;%Change if need!!

[Filename,path] = uigetfile('.mat','Please select fluorescence time series .mat file');
cd(path);
SeriesFileDirectory= [path,Filename];
load(SeriesFileDirectory);

[OriginalCellNo, TimepointNo]=size(dFF);

% [Mark, C]= kmeans(spPos,GroupNo);
% roi2=[];
% for i=1:GroupNo
%     roi2(i).name=strcat('Clus',num2str(i),'-',num2str(length(find(Mark==i))));
%     roi2(i).members=zeros(size(Mark,1),1);
%     roi2(i).members(find(Mark==i))=1;
%     roi2(i).members=logical(roi2(i).members);
% end
% roi=roi2;
% save(strcat(path,'\DistClustered-',Filename,'.mat'),'roi');
% 
% [Mark, C]= kmeans(dFF,GroupNo);
% roi2=[];
% for i=1:GroupNo
%     roi2(i).name=strcat('Clus',num2str(i),'-',num2str(length(find(Mark==i))));
%     roi2(i).members=zeros(size(Mark,1),1);
%     roi2(i).members(find(Mark==i))=1;
%     roi2(i).members=logical(roi2(i).members);
% end
% roi=roi2;
% save(strcat(path,'\dFFClustered-',Filename,'.mat'),'roi'); 
%% -----------Or, not using kmeans way--------------
%-------Calculate Dist matrix-----------
Mark=zeros(size(dFF,1),1);

Dist=squareform(pdist(spPos));
DistThreshold=15;%Changable!!

CloseCheck=Dist;
CloseCheck(CloseCheck<=DistThreshold)=1;
CloseCheck(CloseCheck>DistThreshold)=0;

j=1;
for i=1:OriginalCellNo
    List=find(CloseCheck(i,:)==1);
    if length(List)>1
        if max(Mark(List))>0
            Mark(List)=max(Mark(List));
        else
            Mark(List)=j;
            j=j+1;
        end
    else
        Mark(i)=j;
        j=j+1;
    end
end

%----------Get Group Dimensions---------------
GroupNo=max(Mark)
GroupSize=zeros(size(GroupNo));
for i2=1:GroupNo
    GroupSize(i2)=length(find(Mark==i2));
end

MaxGroupSize=max(GroupSize)

% %--------------Get Mean dFF of each Group, for later plotting-------
% MeandFFOfCurrentGroup=zeros(GroupNo,TimepointNo);
% for i3=1:GroupNo
%     CellInGroup=find(Mark==i3);
%     MeandFFOfCurrentGroup(i3,:)=mean(dFF(CellInGroup,:));
% end
%-------------Save Cluster ROI files-------------
clear roi
for i6=1:GroupNo
    roi(i6).name=strcat('Clus',num2str(i6),'-',num2str(length(find(Mark==i6))));
    roi(i6).members=zeros(size(Mark,1),1);
    roi(i6).members(find(Mark==i6))=1;
    roi(i6).members=logical(roi(i6).members);
end
save(strcat(path,'\ManualDistClust-',Filename,'.mat'),'roi'); 


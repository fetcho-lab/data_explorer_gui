%-------trying corr (clustering) map in total population------------
clear
clc

CorrThreList=[0.97]%%%Set Corr Thre For Clustering here!!
CorrType='Pearson';

list_of_directories = {...
                'J:\Cui-FilesOnImarisComputer\ProcessedData 062118\ROI\F1T2'...
%                 'J:\Cui\ProcessedData 062118\HindBrainCropData'...
%                 'J:\Cui\ProcessedData 061418\HindBrainCropData'...
%                 'J:\Cui\ProcessedData 060518\HindBrainCropData'...
%                 'J:\Cui\ProcessedData 051718(side)-052118\HindBrainCropData'...
    };

for directory_idx  = 1:numel(list_of_directories)
    cd(list_of_directories{directory_idx});
    
    disp(sprintf('Converting %s',list_of_directories{directory_idx}));
    
    DataList=dir('*mat');

    for fileseq=1:size(DataList,1)
%--------------Load data-------------------
        path=list_of_directories{directory_idx};
        name=DataList(fileseq).name;
        filepath=strcat(path,'\',name);
        Data=load(filepath);

%         assignin('base','filepath',filepath);
        fts = Data.fluorescence_time_series;
        dFF = Data.dFF;
        spPos = Data.spPos;
        Sc = Data.Sc;
        if isfield(Data,'spRadiiXYZ')
            spRadiiXYZ=Data.spRadiiXYZ;
        end

        [OriginalCellNo, TimepointNo]=size(dFF);
        

% %-----------Calculate Correlation-----------
        Mark=zeros(size(dFF,1),1);
        i=1;
        j=1;
        dFF2=dFF;
        HighCorrCheck=zeros(size(dFF,1),size(dFF,1));
        if strcmp(CorrType,'Pearson')
            print1='Start Calculating Correlation Matrix with Pearson...'
            [Allcor, AllPVal]=corr(dFF2','Type','Pearson');
        elseif strcmp(CorrType,'Spearman')
            print1='Start Calculating Correlation Matrix with Spearman...'
            [Allcor, AllPVal]=corr(dFF2','Type','Spearman');
        else
            warndlg('Wrong Corr Type Input! Please enter Pearson or Spearman');
        end
        print1='Calculation Done!'
        
        for corrlisti=1:length(CorrThreList)
            CorrThre=CorrThreList(corrlisti);
%-----------Clustering-----------
            HighCorrCheck=Allcor;
            HighCorrCheck(HighCorrCheck>=CorrThre)=1;
            HighCorrCheck(HighCorrCheck<CorrThre)=0;

            for i=1:OriginalCellNo
                List=find(HighCorrCheck(i,:)==1);
                if length(List)>1
                    if max(Mark(List))>0
                        Mark(List)=max(Mark(List));
                    else
                        Mark(List)=j;
                        j=j+1;
                    end
            %         HighCorrCheck(:,i)=0;
                end
            end
    %----------Get Group Dimensions---------------
            GroupNo=max(Mark)
            GroupSize=zeros(size(GroupNo));
            for i2=1:GroupNo
                GroupSize(i2)=length(find(Mark==i2));
            end

            MaxGroupSize=max(GroupSize)
            %--------------Get Mean dFF of each Group, for later plotting-------
            MeandFFOfCurrentGroup=zeros(GroupNo,TimepointNo);
            for i3=1:GroupNo
                CellInGroup=find(Mark==i3);
                MeandFFOfCurrentGroup(i3,:)=mean(dFF(CellInGroup,:));
            end

    %---------Plot and Save Overall Cluster Figures-----------
%             NewFolderName=strcat(num2str(CorrThre),'ThrClusterPlotOf',name(1:end-4));
            NewFolderName=strcat(num2str(CorrThre),'Thr',CorrType,'ClustOf',name(1:end-4));
            mkdir(NewFolderName);
% 
%             Fig1=figure;
%             scatter3(spPos(:,1),spPos(:,2),spPos(:,3),20,Mark,'.');
%             view(0,90);
%             colormap jet;
%             lim = caxis;
%             title(strcat('All Cluster Plot of',name,'. CorrThre:',num2str(CorrThre),'. CorrType:',CorrType))
%             savefig(Fig1,strcat(path,'\',NewFolderName,'\AllClusterPlot-',CorrType,'.fig'));
% 
%             for i4=1:GroupNo
%                 CellInGroup=find(Mark==i4);
%                 %-------------Save cluster roi----------
%                 roi.members=zeros(OriginalCellNo,1);
%                 roi.members(CellInGroup)=1;
%                 roi.members=logical(roi.members);
%                 roi.name=strcat('Cluster',num2str(i4));
%                 save(strcat(path,'\',NewFolderName,'\Cluster',num2str(i4),'-Size',num2str(length(CellInGroup)),'-roi.mat'),'roi');
%                 %----------Plot cluster position and mean dFF--------
%                 if length(CellInGroup)>5
%                     Fig2=figure;
%                     scatter3(spPos(:,1),spPos(:,2),spPos(:,3),20,'.','MarkerEdgeColor', [0.5 0.5 0.5]);
%                     hold on
%                     scatter3(spPos(CellInGroup,1),spPos(CellInGroup,2),spPos(CellInGroup,3),500,'r.');
%                     view(0,90);
%                     title(strcat('Plot of Cluster ',num2str(i4),' in ', name,'-Size',num2str(length(CellInGroup)),'. CorrThre:',num2str(CorrThre),'. CorrType:',CorrType));
%                     savefig(Fig2,strcat(path,'\',NewFolderName,'\PlotOfCluster',num2str(i4),'-Size',num2str(length(CellInGroup)),'-',CorrType,'.fig'));
% 
%                     Fig3=figure;
%                     plot((1:TimepointNo),MeandFFOfCurrentGroup(i4,:));
%                     title(strcat('Mean dFF of Cluster ',num2str(i4),' in ', name,'-Size',num2str(length(CellInGroup)),'. CorrThre:',num2str(CorrThre),'. CorrType:',CorrType));
%                     savefig(Fig3,strcat(path,'\',NewFolderName,'\MeadDFFOfCluster',num2str(i4),'-Size',num2str(length(CellInGroup)),'-',CorrType,'.fig'));
%                 end
%             end
%             close all
    %------------Save Clustring Result-------------------
%             save(strcat(path,'\',NewFolderName,'\ClusteringResult.mat'),'Mark','MeandFFOfCurrentGroup','HighCorrCheck','GroupNo','MaxGroupSize','CorrThre','CorrType');
            save(strcat(path,'\',NewFolderName,'\ClusteringResult.mat'),'Mark','HighCorrCheck','GroupNo','MaxGroupSize','CorrThre','CorrType');
            %-------------Save Cluster ROI files-------------
            clear roi
            for i6=1:GroupNo
            roi(i6).name=strcat('Clus',num2str(i6),'-',num2str(length(find(Mark==i6))));
            roi(i6).members=zeros(size(Mark,1),1);
            roi(i6).members(find(Mark==i6))=1;
            roi(i6).members=logical(roi(i6).members);
            end
            save(strcat(path,'\',NewFolderName,'\ROI-',NewFolderName,'.mat'),'roi');
            %------------Save Cluster ROI, but only when size>=5 or 1<dFF<5------
            clear roi
            j2=1
            for i7=1:GroupNo
                if length(find(Mark==i7))>0
                MaxdFFofGroup=max(max(dFF(find(Mark==i7),:)));
                    if length(find(Mark==i7))>=5
                    roi(j2).name=strcat('Clus',num2str(i7),'-',num2str(length(find(Mark==i7))));
                    roi(j2).members=zeros(size(Mark,1),1);
                    roi(j2).members(find(Mark==i7))=1;
                    roi(j2).members=logical(roi(j2).members);
                    j2=j2+1;
                    elseif (MaxdFFofGroup<5&&MaxdFFofGroup>1)==1
                    roi(j2).name=strcat('Clus',num2str(i7),'-',num2str(length(find(Mark==i7))));
                    roi(j2).members=zeros(OriginalCellNo,1);
                    roi(j2).members(find(Mark==i7))=1;
                    roi(j2).members=logical(roi(j2).members);
                    j2=j2+1;
                    end
                end
            end
            save(strcat(path,'\',NewFolderName,'\ScreenedROI-',NewFolderName,'.mat'),'roi');
            
        end
        print=strcat('end CorrThreList',num2str(corrlisti))
    end
    
    print=strcat('Finished', name,'!')
end

print='All Done!'

    % for i=1:15
    %     if Mark(i)==0
    %         TargetSeriesNo=i;%max(find(Mark>0))+1;
    %         TargetSeries=dFF2(TargetSeriesNo,:);
    %         [SortedCorrVal, SortedCorrCellNo, AllCorrToSeries, AllPValtoSeries,LowdFFList]=TopkCorr_Pearson(dFF2(:,LowerTimeRange:UpperTimeRange), TargetSeries(:,LowerTimeRange:UpperTimeRange), TargetSeriesNo, k,threshold);
    %         TotalCorrMatrix=[TotalCorrMatrix;AllCorrToSeries'];
    %         HighCorrValCell=find(AllCorrToSeries>CorrThre);
    %         if length(HighCorrValCell)>1
    %             MaxMarkVal=max(Mark(HighCorrValCell));
    %             if MaxMarkVal>0
    %                 Mark(HighCorrValCell)=MaxMarkVal;
    %             else
    %                 Mark(HighCorrValCell)=j;
    %                 j=j+1;
    %             end
    %         dFF2(HighCorrValCell,:)=zeros(length(HighCorrValCell),size(dFF2,2));
    %         end
    %     end
    %     
    %         if mod(i,100)==0
    %             i
    %         end
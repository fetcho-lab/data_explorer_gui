%Formal Kmeans Clustering
clear
clc

CorrThre=0.90;%%%Set Corr Thre For Clustering here!!
CorrType='Pearson';

if strcmp(CorrType,'Pearson')
    print1='Start Calculating Correlation with Pearson...'
elseif strcmp(CorrType,'Spearman')
    print1='Start Calculating Correlation with Spearman...'
else
    warndlg('Wrong Corr Type Input! Please enter Pearson or Spearman');
end
    
list_of_directories = {...
                'J:\Cui-FilesOnImarisComputer\ProcessedData 062118\ROI\F1T2'...
%                 'J:\Cui\ProcessedData 060518\dFFComputed\HindBrainCrop'...
%                 'J:\Cui\ProcessedData 061418-v2\ProcessedData061418 - dFFComputed\HindBrainCrop'...
%                 'J:\Cui\ProcessedData 062118\Data-dFFComputed\HindBrainCrop'...
%                 'E:\Dawnis\IndependenceDay\L07\SptSTructural\_20170704_131820'...
    };


for directory_idx  = 1:numel(list_of_directories)
    cd(list_of_directories{directory_idx});
    
    disp(sprintf('Converting %s',list_of_directories{directory_idx}));
    
    DataList=dir('*mat');

    for fileseq=1:2%size(DataList,1)
%--------------Load data-------------------
        path=list_of_directories{directory_idx};
        name=DataList(fileseq).name;
        filepath=strcat(path,'\',name);
        Data=load(filepath);

%         assignin('base','filepath',filepath);
        fts = Data.fluorescence_time_series;
        dFF = Data.dFF;%(1:1000,:); %Chenge this back!
        spPos = Data.spPos;
        Sc = Data.Sc;
        if isfield(Data,'spRadiiXYZ')
            spRadiiXYZ=Data.spRadiiXYZ;
        end

        [OriginalCellNo, TimepointNo]=size(dFF);
        ClusterNo=round(OriginalCellNo/10);%Set Cluster no here!!

%-----------Calculate Correlation-----------

        Mark=kmeans(dFF, ClusterNo);

%         Centroid=randi([1 OriginalCellNo],1,ClusterNo);%Set initial cluster centroids randomly
%         CentroiddFF=dFF(Centroid,:);
%         
%         for i=1:OriginalCellNo
%             CorToCentroidList=zeros(1,ClusterNo);
%             PValToCentroidList=zeros(1,ClusterNo);
%             for j=1:ClusterNo
%                 if strcmp(CorrType,'Pearson')
%                     [CorToCentroidList(j), PValToCentroidList(j)]=corr(dFF(i,:),CentroiddFF(j,:),'Type','Pearson');
%                 elseif strcmp(CorrType,'Spearman')
%                     [CorToCentroidList(j), PValToCentroidList(j)]=corr(dFF(i,:),CentroiddFF(j,:),'Type','Spearman');
%                 else
%                     warndlg('Wrong Corr Type Input! Please enter Pearson or Spearman');
%                 end
%             end
%             %NOT FINISHED!
            
            
            
            
%% ----------------Old Informal Version---------------
%         Mark=zeros(size(dFF,1),1);
%         i=1;
%         j=1;
%         dFF2=dFF;
%         HighCorrCheck=zeros(size(dFF,1),size(dFF,1));
%         if strcmp(CorrType,'Pearson')
%             print1='Start Calculating Correlation Matrix with Pearson...'
%             [Allcor, AllPVal]=corr(dFF2','Type','Pearson');
%         elseif strcmp(CorrType,'Spearman')
%             print1='Start Calculating Correlation Matrix with Spearman...'
%             [Allcor, AllPVal]=corr(dFF2','Type','Spearman');
%         else
%             warndlg('Wrong Corr Type Input! Please enter Pearson or Spearman');
%         end
%         print1='Calculation Done!'
% %-----------Clustering-----------
%         HighCorrCheck=Allcor;
%         HighCorrCheck(HighCorrCheck>=CorrThre)=1;
%         HighCorrCheck(HighCorrCheck<CorrThre)=0;
% 
%         for i=1:OriginalCellNo
%             List=find(HighCorrCheck(i,:)==1);
%             if length(List)>1
%                 if max(Mark(List))>0
%                     Mark(List)=max(Mark(List));
%                 else
%                     Mark(List)=j;
%                     j=j+1;
%                 end
%         %         HighCorrCheck(:,i)=0;
%             end
%         end
% %----------Get Group Dimensions---------------
        GroupNo=max(Mark)
        GroupSize=zeros(size(GroupNo));
        for i2=1:GroupNo
            GroupSize(i2)=length(find(Mark==i2));
        end

        MaxGroupSize=max(GroupSize)
        MeandFFOfCurrentGroup=zeros(GroupNo,TimepointNo);
        for i3=1:GroupNo
            CellInGroup=find(Mark==i3);
            MeandFFOfCurrentGroup(i3,:)=mean(dFF(CellInGroup,:));
        end
        
% %---------Plot and Save Cluster Figures-----------
        NewFolderName=strcat(num2str(CorrThre),'ThrClusterPlotOf',name(1:end-4),'v2');
        mkdir(NewFolderName);

        Fig1=figure;
        scatter3(spPos(:,1),spPos(:,2),spPos(:,3),20,Mark,'.');
        view(0,90);
        colormap jet;
        lim = caxis;
        title(strcat('All Cluster Plot of',name,'. CorrThre:',num2str(CorrThre),'. CorrType:',CorrType,'v2'))
        savefig(Fig1,strcat(path,'\',NewFolderName,'\AllClusterPlot-',CorrType,'v2','.fig'));

        for i4=1:GroupNo
            CellInGroup=find(Mark==i4);
            if length(CellInGroup)>10
                Fig2=figure;
                scatter3(spPos(:,1),spPos(:,2),spPos(:,3),20,'.','MarkerEdgeColor', [0.5 0.5 0.5]);
                hold on
                scatter3(spPos(CellInGroup,1),spPos(CellInGroup,2),spPos(CellInGroup,3),500,'r.');
                view(0,90);
                title(strcat('Plot of Cluster ',num2str(i4),' in ', name,'-Size',num2str(length(CellInGroup)),'. CorrThre:',num2str(CorrThre),'. CorrType:',CorrType));
                savefig(Fig2,strcat(path,'\',NewFolderName,'\PlotOfCluster',num2str(i4),'-Size',num2str(length(CellInGroup)),'-',CorrType,'v2','.fig'));
                
                Fig3=figure;
                plot((1:TimepointNo),MeandFFOfCurrentGroup(i4,:));
                title(strcat('Mean dFF of Cluster ',num2str(i4),' in ', name,'-Size',num2str(length(CellInGroup)),'. CorrThre:',num2str(CorrThre),'. CorrType:',CorrType));
                savefig(Fig3,strcat(path,'\',NewFolderName,'\MeadDFFOfCluster',num2str(i4),'-Size',num2str(length(CellInGroup)),'-',CorrType,'v2','.fig'));
            end
        end
        close all
%------------Save Clustring Result-------------------
%         save(strcat(path,'\',NewFolderName,'\ClusteringResult.mat'),'Mark','MeandFFOfCurrentGroup','HighCorrCheck','GroupNo','MaxGroupSize','CorrThre','CorrType');
        save(strcat(path,'\',NewFolderName,'\ClusteringResult.mat'),'Mark','MeandFFOfCurrentGroup','GroupNo','MaxGroupSize');
    print=strcat('Finished', name,'!')
    end
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
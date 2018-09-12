clear all
% -----To be adjusted-----
Data=load ('K:\Cui\20180322_180833-GCaMP6s\GCaMP6 fish-no ethnol.mat');%the directory of the file to be processed
Cor_thre=0.995 %thereshold of correlation
disthre= [0,999]; %threshold of distance between correlated cells, a two elements array, first one is lower threshold second one is upper threshold.

use_dFF=1;

Discarded_Fluo_Int_Percent=0;
Discarded_dFF_Percent=0;
weird_plane=0;

SizeRangeForDots=[20,100];%preferably bigger size value for dots than circles, circles are default in scatter3!
WidthRangeForLines=[5,30];
% ----------------------------------

Allpos=Data.spPos;
Allfluo=Data.fluorescence_time_series;
Allscale=Data.Sc;
Cellno=size(Allfluo,1);
TimepointNo=size(Allfluo,2);


%% ---------Calculate dF/F----------------
if use_dFF==1
    %---------To be adjusted-------
    BasalFVolume=20;%This equals how many lowest fluo values (for each cell) will be counted as basal value
    %--------------------------------
    sortedfluo=sort(Allfluo,2,'ascend');
    basalF= mean(sortedfluo(:,1:BasalFVolume),2);
    basalFexpand=basalF*ones(1,TimepointNo);
    deltaF=Allfluo-(basalFexpand);
    AlldFF=deltaF./basalFexpand;
    save ('dFF_of_all_cells','AlldFF','basalF','BasalFVolume');
end


%% ---------Get rid of very dim cells first-----------
fluomin=min(Allfluo,[],2);
fluomax=max(Allfluo,[],2);
MinFluoOfAllCells=min(fluomin);
MaxFluoOfAllCells=max(fluomax);
Fluochangethre=MinFluoOfAllCells+Discarded_Fluo_Int_Percent*(MaxFluoOfAllCells-MinFluoOfAllCells);

fluochange=abs((fluomax-fluomin));
disc_list=find(fluochange<Fluochangethre);
Allpos(disc_list,:)=[];
Allfluo(disc_list,:)=[];
if use_dFF==1
    AlldFF(disc_list,:)=[];
end
% Above: Get rid of cells with too low fluo change (absolute value)

%% ------------Get rid of low dFF cells----------------
if use_dFF==1
    dFFmin=min(AlldFF,[],2);
    dFFmax=max(AlldFF,[],2);
    MaxdFFOfAllCells=max(max(AlldFF));
    MindFFOfAllCells=min(min(AlldFF));
    dFFchangethre=MindFFOfAllCells+Discarded_dFF_Percent*(MaxdFFOfAllCells-MindFFOfAllCells);

    dFFchange=abs((dFFmax-dFFmin));
    disc_list=find(dFFchange<dFFchangethre);
    Allpos(disc_list,:)=[];
    Allfluo(disc_list,:)=[];
    AlldFF(disc_list,:)=[];
end

%% --------Get rid of weird cells in weird plane if needed----------------------

if weird_plane==1
 % --------To be adjusted----------
suspectz=[102.5,92.5,82.5,107.5,127.5];
z_stepsize=5;
Distancethre=20;%if the minimal distance between a cell in suspected plane and any other cell in plane above or below is further than this value, that cell in suspected plane will be discarded.
%----------------------------------
for i6=1:length(suspectz)
    ind_sus=find (Allpos(:,3)==suspectz(i6));
    ind_abo_sus=find (Allpos(:,3)==suspectz(i6)-z_stepsize);
    ind_belo_sus=find (Allpos(:,3)==suspectz(i6)+z_stepsize);
    sus_xyz=Allpos(ind_sus,:);
    abo_sus_xyz=Allpos(ind_abo_sus,:);
    belo_sus_xyz=Allpos(ind_belo_sus,:);
    %calculate distance
    min_dis_to_abo=min(pdist2(abo_sus_xyz,sus_xyz));
    min_dis_to_belo=min(pdist2(belo_sus_xyz,sus_xyz));
    min_of_min=min([min_dis_to_abo; min_dis_to_belo]);
    templist=find(min_of_min>Distancethre);
    disc_list2=ind_sus(templist);
    Allpos(disc_list2,:)=[];
    Allfluo(disc_list2,:)=[];
    if use_dFF==1
        AlldFF(disc_list,:)=[];
    end
end
save ('SuspectedPlane&DistanceThre', 'suspectz','z_stepsize','Distancethre');
end
%% --------Plot Filtered Cell Positions-------
Maxfluoofcells1=max(Allfluo');

totaldifofmaxfluo=max(Maxfluoofcells1)-min(Maxfluoofcells1);
% Maxfluoofcells2=max(Data.fluorescence_time_series');
% Minfluoofcells2=min(Data.fluorescence_time_series');
sizeofpoints=((Maxfluoofcells1-min(Maxfluoofcells1))./totaldifofmaxfluo)*(SizeRangeForDots(2)-SizeRangeForDots(1))+SizeRangeForDots(1);

Allposx=Allpos(:,1);
Allposy=Allpos(:,2);
Allposz=Allpos(:,3);
%-----------To Be Adjusted--------------
scatterplot1=scatter3(Allposx,Allposy,Allposz,sizeofpoints,...
    Maxfluoofcells1,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)%Adjust this line for different plotting settings!%
%--------------------------------------
xlabel('X');
ylabel('Y');
zlabel('Z');
colormap jet;
CB1=colorbar;
CB1.Label.String='Fluorescence intensity';
FluoColorbarLim1=caxis;
hold on
if use_dFF==1
    save ('AllPos&FluoAfterFilter.mat', 'Allfluo', 'Allpos','AlldFF','Discarded_Fluo_Int_Percent','Fluochangethre','Discarded_dFF_Percent','dFFchangethre');
else
    save ('AllPos&FluoAfterFilter.mat', 'Allfluo', 'Allpos','Discarded_Fluo_Int_Percent','Fluochangethre');
end
save('ParametersOfFilteredCellMap','FluoColorbarLim1','SizeRangeForDots');
savefig('FilteredCellMap.fig');


%% -------------Calculate correlation-------------------
if use_dFF==1
    %------------dFF version--------------------
    %----------To be adjusted--------------
    clipno=size(AlldFF,1); %how many nomber of cells you want to process and compare (may be limited by computer RAM)
    %--------------------------------------
    clippeddFF=AlldFF(1:clipno, :);
    %for i=1:Cellno
       %Cor=xcorr(Allfluo(:, 1:3))
       Allcor2=corrcoef(clippeddFF');
       for i2=1:clipno
           Allcor2(i2,i2)=0;
       end
           [Maxcor2, Maxcorno2]=max(Allcor2);
    %------------dFF version end-------------
else
    %-------------AllFluo version----------------
    %----------To be adjusted--------------
    clipno=size(Allfluo,1); %how many nomber of cells you want to process and compare (may be limited by computer RAM)
    %--------------------------------------
    clippedfluo=Allfluo(1:clipno, :);
    %for i=1:Cellno
       %Cor=xcorr(Allfluo(:, 1:3))
       Allcor1=corrcoef(clippedfluo');
       for i2=1:clipno
           Allcor1(i2,i2)=0;
       end
           [Maxcor1, Maxcorno1]=max(Allcor1);
    %------------AllFluo version end------------
end

%% ------------Plot correlation lines------------

    allx=Allpos(:,1);%*Allscale(1,1);
    ally=Allpos(:,2);%*Allscale(2,2);
    allz=Allpos(:,3);%*Allscale(3,3);
    
for i5=1:length(disthre)-1
% %----------plot all cell's position if you like-------------------
% if i5==1
%     Maxfluoofcells2=max(Data.fluorescence_time_series');
%     maxmindif2=max(Maxfluoofcells2)-min(Maxfluoofcells2);
% %-------------To Be Adjusted--------------
%     SizeRangeForDots2=[20,500];%Adjustable!!!
% %-----------------------------------------
%     sizeofdots2=((Maxfluoofcells2-min(Maxfluoofcells2))./maxmindif2)*(SizeRangeForDots2(2)-SizeRangeForDots2(1))+SizeRangeForDots2(1);
% 
%     scatterplot1=scatter3(Data.spPos(:,1),Data.spPos(:,2),Data.spPos(:,3),sizeofdots2,...
%         Maxfluoofcells2,'.');%Adjust this line for different plotting settings!%,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2
%     %     'MarkerFaceColor','b','MarkerEdgeColor','b')<-these
%     %     are old settings.
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     colormap jet;
%     CB2=colorbar;
%     CB2.Label.String='Fluorescence intensity';
%     FluoColorbarLim2=caxis;
% 
%     save('ParametersOfAllCellMap','FluoColorbarLim2','SizeRangeForDots2');
%     savefig('AllCellMap.fig');
%     hold on
% end
% 
% %---------------------------------------------------------------------------
    Plotted_Correlation_No=0;
    if use_dFF==1
        for i3=1:clipno
            tempx=[allx(i3,1),allx(Maxcorno2(i3),1)];
            tempy=[ally(i3,1),ally(Maxcorno2(i3),1)];
            tempz=[allz(i3,1),allz(Maxcorno2(i3),1)];
            dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
            if dist>disthre(i5)
                if dist<disthre(i5+1)
                    if Maxcor2(i3)>Cor_thre
                        if Maxcor2(i3)<1
                        linewidth=(((Maxcor2(i3)-Cor_thre)/(1-Cor_thre))*(WidthRangeForLines(2)-WidthRangeForLines(1)))+WidthRangeForLines(1);
                        colorR=max(AlldFF(i3,:))/MaxdFFOfAllCells;
                        colorG=max(AlldFF(Maxcorno2(i3),:))/MaxdFFOfAllCells;
                        colorB=0;
                        Corrplot=plot3(tempx, tempy, tempz,[colorR colorG colorB],'-','LineWidth',linewidth);
                        hold on
                        Plotted_Correlation_No=Plotted_Correlation_No+1;
                        end
                    end
                end
            end
        end
        hold off
        Plotted_Correlation_No

        save (strcat('Cor_Map_Data_',num2str(disthre(i5)),'-',num2str(disthre(i5+1)),'_',num2str(Cor_thre),'_dFF.mat'), 'Maxcor2', 'Maxcorno2', 'Cor_thre', 'Plotted_Correlation_No','disthre');
        savefig(strcat('Cor_Map_',num2str(disthre(i5)),'-',num2str(disthre(i5+1)),'_',num2str(Cor_thre),'_dFF.fig'));
    else
        for i3=1:clipno
            tempx=[allx(i3,1),allx(Maxcorno1(i3),1)];
            tempy=[ally(i3,1),ally(Maxcorno1(i3),1)];
            tempz=[allz(i3,1),allz(Maxcorno1(i3),1)];
            dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
            if dist>disthre(i5)
                if dist<disthre(i5+1)
                    if Maxcor1(i3)>Cor_thre
                        if Maxcor1(i3)<1
                        linewidth=(((Maxcor1(i3)-Cor_thre)/(1-Cor_thre))*(WidthRangeForLines(2)-WidthRangeForLines(1)))+WidthRangeForLines(1);
                        colorR=max(Allfluo(i3,:))/MaxFluoOfAllCells;
                        colorG=max(Allfluo(Maxcorno1(i3),:))/MaxFluoOfAllCells;
                        colorB=0;
                        Corrplot=plot3(tempx, tempy, tempz, '-','LineWidth',linewidth);
                        hold on
                        Plotted_Correlation_No=Plotted_Correlation_No+1;
                        end
                    end
                end
            end
        end
        hold off
        Plotted_Correlation_No

        save (strcat('Cor_Map_Data_',num2str(disthre(i5)),'-',num2str(disthre(i5+1)),'_',num2str(Cor_thre),'_Fluo.mat'), 'Maxcor1', 'Maxcorno1', 'Cor_thre', 'Plotted_Correlation_No','disthre');
        savefig(strcat('Cor_Map_',num2str(disthre(i5)),'-',num2str(disthre(i5+1)),'_',num2str(Cor_thre),'_Fluo.fig'));
    end
end

% xlswrite('CorrelationMap.xlsx', Maxcor',1,'A1');
% xlswrite('CorrelationMap.xlsx', Maxcorno',2,'A1');
% xlswrite('CorrelationMap.xlsx', corthre,3,'A2');
clear all
% -----To be adjusted-----
Data=load ('M:\Joes_Data_Explorer\ls_fluorescence_time_series_Elb_012717Fish01Trial02_90Hz.mat');%the directory of the file to be processed
Cor_thre=0.995 %thereshold of correlation
Dis_thre= [0,999]; %threshold of distance between correlated cells, a two elements array, first one is lower threshold second one is upper threshold.

use_dFF_For_Cor_Map=0;%set this to 1 means the program will do calculate correlations based on the dFF of each cell, not fluorescence.
use_dFF_For_Cell_Map=0;%set this to 1 means the parameters of cell map (dots size&color) will be calculated based on the dFF of each cell, not fluorescence. 0 is recommend due to the high variability of dFF.
BasalFVolume=20;%This equals how many lowest fluo values (for each cell) will be counted as basal value

Discarded_Fluo_Int_Percent=0;
Discarded_dFF_Percent=0;
weird_plane=0;

Plot_All_Cells_Position_Before_Plot_Cor=0;
SizeRangeForDots=[5,100];%preferably bigger size value for dots than circles, circles are default in scatter3!
TransparencyForDots=0.2;
WidthRangeForLines=[2,10];
dFFLineColorMaxVal=5;%the R value of connecting line's color will be one cell's dFF/this value, and G value will be the other's dFF/this value. if dFF/this value>1 then it will be set to 1. So this is also the upper limit of dFF that can be showned by line color.
% FluoLineColorMaxVal can be set at around line 241.
% ----------------------------------

Allpos=Data.spPos;
Allfluo=Data.fluorescence_time_series;
Allscale=Data.Sc;
Cellno=size(Allfluo,1);
TimepointNo=size(Allfluo,2);


%% ---------Calculate dF/F----------------
if use_dFF_For_Cor_Map==1
    %---------To be adjusted-------
    %--------------------------------
    sortedfluo=sort(Allfluo,2,'ascend');
    basalF= mean(sortedfluo(:,1:BasalFVolume),2);
    ZeroBasalFCells_discarded=find(basalF==0);
    if size(ZeroBasalFCells_discarded)>0
        warndlg(strcat(num2str(size(ZeroBasalFCells_discarded)),' basal F values contains 0, try to make Basal Volume higher? Cells contains 0 basal F value will be deleted.'));
    end
    basalFexpand=basalF*ones(1,TimepointNo);
    deltaF=Allfluo-(basalFexpand);
    AlldFF=deltaF./basalFexpand;
    
    AlldFF(ZeroBasalFCells_discarded,:)=[];
    Allpos(ZeroBasalFCells_discarded,:)=[];
    Allfluo(ZeroBasalFCells_discarded,:)=[];
    save ('dFF_of_all_cells','AlldFF','basalF','BasalFVolume','ZeroBasalFCells_discarded');
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
if use_dFF_For_Cor_Map==1
    AlldFF(disc_list,:)=[];
end
% Above: Get rid of cells with too low fluo change (absolute value)

%% ------------Get rid of low dFF cells----------------
if use_dFF_For_Cor_Map==1
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
        if use_dFF_For_Cor_Map==1
            AlldFF(disc_list,:)=[];
        end
    end
    save ('SuspectedPlane&DistanceThre', 'suspectz','z_stepsize','Distancethre');
end
%% --------Plot Filtered Cell Positions-------
if use_dFF_For_Cell_Map==1
    PlotCellMap(Allpos, use_dFF_For_Cell_Map, AlldFF,SizeRangeForDots, TransparencyForDots);%Sequence of Input Arguments: (Allpos, use_dFF, AllfluoOrAlldFF,SizeRangeForDots, TransparencyForDots). You can set use_dFF to 0 if you want to use fluo for cell map.
else
    PlotCellMap(Allpos, use_dFF_For_Cell_Map, Allfluo,SizeRangeForDots, TransparencyForDots);%Sequence of Input Arguments: (Allpos, use_dFF, AllfluoOrAlldFF,SizeRangeForDots, TransparencyForDots). You can set use_dFF to 0 if you want to use fluo for cell map.
end
hold on
% if use_dFF==1
%     MaxdFFOfCells=max(AlldFF');
%     ColorOfDots=MaxdFFOfCells;
%     TotalDifOfMaxdFF=max(MaxdFFOfCells)-min(MaxdFFOfCells);
%     SizeOfDots=((MaxdFFOfCells-min(MaxdFFOfCells))./TotalDifOfMaxdFF)*(SizeRangeForDots(2)-SizeRangeForDots(1))+SizeRangeForDots(1);
% else
%     MaxFluoOfCells=max(Allfluo');
%     ColorOfDots=MaxFluoOfCells;
%     TotalDifOfMaxFluo=max(MaxFluoOfCells)-min(MaxFluoOfCells);
%     SizeOfDots=((MaxFluoOfCells-min(MaxFluoOfCells))./TotalDifOfMaxFluo)*(SizeRangeForDots(2)-SizeRangeForDots(1))+SizeRangeForDots(1);
% end
% 
% Allposx=Allpos(:,1);
% Allposy=Allpos(:,2);
% Allposz=Allpos(:,3);
% %-----------To Be Adjusted--------------
% scatterplot1=scatter3(Allposx,Allposy,Allposz,SizeOfDots,...
%     ColorOfDots,'MarkerFaceAlpha',TransparencyForDots,'MarkerEdgeAlpha',TransparencyForDots)%Adjust this line for different plotting settings!%
% %--------------------------------------
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% colormap jet;
% CB1=colorbar;
% if use_dFF==1
%     CB1.Label.String='dF/F of cells (dots)';
% else
%     CB1.Label.String='Fluorescence intensity of cells (dots)';
% end
% FluoColorbarLim1=caxis;
% hold on
% if use_dFF==1
%     save('FilteredCellMapPlotParameters_dFF.mat','use_dFF','FluoColorbarLim1','SizeRangeForDots','TransparencyForDots');
%     savefig('FilteredCellMap_dFF.fig');
% else
%     save('FilteredCellMapPlotParameters_Fluo.mat','use_dFF','FluoColorbarLim1','SizeRangeForDots','TransparencyForDots');
%     savefig('FilteredCellMap_Fluo.fig');
% end
%-----------------Save Pos&Fluo data after filter----------------------
if use_dFF_For_Cor_Map==1
    save ('AllPos&FluoAfterFilter_dFF.mat','use_dFF_For_Cell_Map', 'use_dFF_For_Cor_Map','Allfluo', 'Allpos','AlldFF','Discarded_Fluo_Int_Percent','Fluochangethre','Discarded_dFF_Percent','dFFchangethre');
else
    save ('AllPos&FluoAfterFilter_Fluo.mat','use_dFF_For_Cell_Map','use_dFF_For_Cor_Map','Allfluo', 'Allpos','Discarded_Fluo_Int_Percent','Fluochangethre');
end
%% -------------Calculate correlation-------------------
if use_dFF_For_Cor_Map==1
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
    
for i5=1:length(Dis_thre)-1
%----------plot all cell's position(un filtered) if you like-------------------
PlotCellMap (Data.spPos, 0, Data.fluorescence_time_series, SizeRangeForDots, TransparencyForDots);
hold on
% if Plot_All_Cells_Position_Before_Plot_Cor==1
%     if i5==1
%         Maxfluoofcells2=max(Data.fluorescence_time_series');
%         maxmindif2=max(Maxfluoofcells2)-min(Maxfluoofcells2);
%         sizeofdots2=((Maxfluoofcells2-min(Maxfluoofcells2))./maxmindif2)*(SizeRangeForDots2(2)-SizeRangeForDots2(1))+SizeRangeForDots2(1);
% 
%         scatterplot1=scatter3(Data.spPos(:,1),Data.spPos(:,2),Data.spPos(:,3),sizeofdots2,...
%             Maxfluoofcells2,'.');%Adjust this line for different plotting settings!%,'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2
%         %     'MarkerFaceColor','b','MarkerEdgeColor','b')<-these
%         %     are old settings.
%         xlabel('X');
%         ylabel('Y');
%         zlabel('Z');
%         colormap jet;
%         CB2=colorbar;
%         CB2.Label.String='Fluorescence intensity';
%         FluoColorbarLim2=caxis;
% 
%         save('ParametersOfAllCellMap','FluoColorbarLim2','SizeRangeForDots2');
%         savefig('AllCellMap.fig');
%         hold on
%     end
% end
    %---------------------------------------------------------------------------

    Plotted_Correlation_No=0;
    if use_dFF_For_Cor_Map==1
%         MaxdFFOfAllCells2=max(max(AlldFF));
        for i3=1:clipno
            tempx=[allx(i3,1),allx(Maxcorno2(i3),1)];
            tempy=[ally(i3,1),ally(Maxcorno2(i3),1)];
            tempz=[allz(i3,1),allz(Maxcorno2(i3),1)];
            dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
            if dist>Dis_thre(i5)
                if dist<Dis_thre(i5+1)
                    if Maxcor2(i3)>Cor_thre
                        if Maxcor2(i3)<1
                        linewidth=(((Maxcor2(i3)-Cor_thre)/(1-Cor_thre))*(WidthRangeForLines(2)-WidthRangeForLines(1)))+WidthRangeForLines(1);
                        colorR=max(AlldFF(i3,:))/dFFLineColorMaxVal;
                        if colorR>1
                            colorR=1
                        end
                        colorG=max(AlldFF(Maxcorno2(i3),:))/dFFLineColorMaxVal;
                        if colorG>1
                            colorG=1
                        end
                        colorB=0;
                        Corrplot=plot3(tempx, tempy, tempz,'-','Color',[colorR, colorG, colorB],'LineWidth',linewidth);
                        hold on
                        Plotted_Correlation_No=Plotted_Correlation_No+1;
                        end
                    end
                end
            end
        end
        hold off
        Plotted_Correlation_No

        save (strcat('CorMapData',num2str(Dis_thre(i5)),'-',num2str(Dis_thre(i5+1)),'_',num2str(Cor_thre),'_dFF.mat'),'use_dFF_For_Cell_Map','use_dFF_For_Cor_Map', 'Maxcor2', 'Maxcorno2', 'Cor_thre','Dis_thre','Plotted_Correlation_No');
        save (strcat('CorMapPlotParameters',num2str(Dis_thre(i5)),'-',num2str(Dis_thre(i5+1)),'_',num2str(Cor_thre),'_dFF.mat'),'use_dFF_For_Cell_Map','use_dFF_For_Cor_Map','WidthRangeForLines','dFFLineColorMaxVal', 'Cor_thre','Dis_thre','Plotted_Correlation_No');
        savefig(strcat('Cor_Map_',num2str(Dis_thre(i5)),'-',num2str(Dis_thre(i5+1)),'_',num2str(Cor_thre),'_dFF.fig'));
    else
        FluoLineColorMaxVal=max(max(Allfluo));%set FluoLineColorMaxVal here!!!!!!
        for i3=1:clipno
            tempx=[allx(i3,1),allx(Maxcorno1(i3),1)];
            tempy=[ally(i3,1),ally(Maxcorno1(i3),1)];
            tempz=[allz(i3,1),allz(Maxcorno1(i3),1)];
            dist=sqrt((tempx(2)-tempx(1))^2+(tempy(2)-tempy(1))^2+(tempz(2)-tempz(1))^2);
            if dist>Dis_thre(i5)
                if dist<Dis_thre(i5+1)
                    if Maxcor1(i3)>Cor_thre
                        if Maxcor1(i3)<1
                        linewidth=(((Maxcor1(i3)-Cor_thre)/(1-Cor_thre))*(WidthRangeForLines(2)-WidthRangeForLines(1)))+WidthRangeForLines(1);
                        colorR=max(Allfluo(i3,:))/FluoLineColorMaxVal;
                        if colorR>1
                            colorR=1;
                        end
                        colorG=max(Allfluo(Maxcorno1(i3),:))/FluoLineColorMaxVal;
                        if colorG>1
                            colorG=1;
                        end
                        colorB=0;
                        Corrplot=plot3(tempx, tempy, tempz, '-','Color',[colorR, colorG, colorB],'LineWidth',linewidth);
                        hold on
                        Plotted_Correlation_No=Plotted_Correlation_No+1;
                        end
                    end
                end
            end
        end
        hold off
        Plotted_Correlation_No

        save (strcat('CorMapData',num2str(Dis_thre(i5)),'-',num2str(Dis_thre(i5+1)),'_',num2str(Cor_thre),'_Fluo.mat'),'use_dFF_For_Cell_Map','use_dFF_For_Cor_Map','Maxcor1', 'Maxcorno1', 'Cor_thre', 'Dis_thre','Plotted_Correlation_No');
        save (strcat('CorMapPlotParameters',num2str(Dis_thre(i5)),'-',num2str(Dis_thre(i5+1)),'_',num2str(Cor_thre),'_Fluo.mat'),'use_dFF_For_Cell_Map','use_dFF_For_Cor_Map','WidthRangeForLines','FluoLineColorMaxVal', 'Cor_thre', 'Dis_thre','Plotted_Correlation_No');
        savefig(strcat('Cor_Map_',num2str(Dis_thre(i5)),'-',num2str(Dis_thre(i5+1)),'_',num2str(Cor_thre),'_Fluo.fig'));
    end
end

% xlswrite('CorrelationMap.xlsx', Maxcor',1,'A1');
% xlswrite('CorrelationMap.xlsx', Maxcorno',2,'A1');
% xlswrite('CorrelationMap.xlsx', corthre,3,'A2');
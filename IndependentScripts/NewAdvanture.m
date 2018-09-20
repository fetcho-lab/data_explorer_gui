clear all


answer=inputdlg({'File Directory','File Name'},'Input Data'); 
FileToRead=strcat(answer{1,1},'\',answer{2,1});
Data=load(FileToRead);

Discarded_dFF_Percent=0.2;
Discarded_Fluo_Percent=0;

Allpos=Data.spPos;
Allfluo=Data.fluorescence_time_series;
Allscale=Data.Sc;
Cellno=size(Allfluo,1);
TimepointNo=size(Allfluo,2);

BasalFVolume=round(Cellno*0.1);

for i=2:size(Allfluo,2)
AllFluoDif(:,i-1)=Allfluo(:,i)-Allfluo(:,i-1);
end

clims=[min(min(AllFluoDif)),max(max(AllFluoDif))];
imagesc(AllFluoDif,clims);
colorbar
colormap jet
ylim=[0 size(Allfluo,1)];

List=zeros(size(AllFluoDif));
for i2=1:size(AllFluoDif,1)
    for j2=1:size(AllFluoDif,2)
        if AllFluoDif(i2,j2)<5
            if Allfluo(i2,j2)>min(Allfluo(i2,:))+0.5*(max(Allfluo(i2,:))-min(Allfluo(i2,:)));
                List(i2,j2)=1;
            end
        end
    end
end
            

%% ----------------Merge Cells---------------------
% [AllMerFluo, AllMerPos]=CellMerger(Allfluo, Allpos,10000,5);
% 
% spPos=AllMerPos;
% fluorescence_time_series=AllMerFluo;
% Sc=Data.Sc;
% save(strcat(answer{1,1},'\merged',answer{2,1}) ,'spPos','fluorescence_time_series','Sc');

%% ------------HeatMapping Brain----------------
    Fluomin=min(Allfluo,[],2);
    Fluomax=max(Allfluo,[],2);
    MaxFluoOfAllCells=max(max(Allfluo));
    MinFluoOfAllCells=min(min(Allfluo));
      
    Fluochange=abs(Fluomax-Fluomin);
    Fluochangethre=min(Fluochange)+Discarded_Fluo_Percent*(max(Fluochange)-min(Fluochange));

    disc_list_Fluo=find(Fluochange<Fluochangethre);
    
    Allpos2=Allpos;
    Allfluo2=Allfluo;

    Allpos2(disc_list_Fluo,:)=[];
    Allfluo2(disc_list_Fluo,:)=[];

    clims=[MinFluoOfAllCells,MaxFluoOfAllCells];
    imagesc(Allfluo2(1:end,:),clims);
    colorbar
    colormap jet
    ylim=[0 size(Allfluo,1)];

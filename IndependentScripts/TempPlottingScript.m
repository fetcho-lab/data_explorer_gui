% Please make sure dFF/Fluo data has already loadded into base workspace!
% [DataFileName,DataFilePath] = uigetfile('.mat','Please select dffData .mat file');
% wholeDataf=strcat(DataFilePath,DataFileName);
% load(wholeDataf);

[ROIFileName,ROIFilePath] = uigetfile('.mat','Please select roi.mat file');
wholeROIf=strcat(ROIFilePath,ROIFileName)
load(wholeROIf);
ROIdFF=dFF(find(roi.members),:);
ROIFluo=fluorescence_time_series(find(roi.members),:);

[StimFileName,StimFilePath] = uigetfile('.mat','Please select Stim.mat file');
wholeStimf=strcat(StimFilePath,StimFileName);
load(wholeStimf);
Stim=StimTimeSeries;

Stackpersec=1.28;%Change this!!
xlim1=[0, size(dFF,2)];

for i=1:3
    xlength=size(dFF,2);

    pos1=[0.1 0.2 0.8 0.1];
    F1=figure;
    left_color = [0 0 0];
    right_color = [1 0 0];
    set(F1,'defaultAxesColorOrder',[left_color; right_color]);

    x1=(1:xlength);
    ax1=subplot('Position',pos1);
    ax1.XGrid='on';
    ax1.YGrid='on';

    plot(x1,Stim);
    xlim(xlim1);%change x limit here!!
    xt=[0:Stackpersec*100:550];
    xticks(xt);
    xticklabels({'0','100','200','300','400'});%Change this label!!
    xlabel('Time(s)');
    ylabel('Stim Amplitude');
    % ax1=gca;
    % ax1pos=ax1.Position;
    % ax2=axes('Position',ax1pos,...
    %     'XAxisLocation','top',...
    %     'YAxisLocation','right',...
    %     'Color','none');


    pos2=[0.1 0.3 0.8 0.6]
    ax2=subplot('Position',pos2);
    ax2.XGrid='on';
    ax2.YGrid='on';
    if i==1
    %--------Plot mean Fluo---------
        plot(x1,mean(ROIFluo));
        xlim(xlim1);
        xt=[0:Stackpersec*100:550];
        xticks(xt);
        xticklabels({'0','100','200','300','400'});%Change this label!!

        yyaxis right
        plot(x1,Stim,'r');
%         title('Mean Fluo of ROI cells');
        
        savefig(F1,strcat(wholeROIf,'-MeanFluo.fig'));

    elseif i==2
    %--------Plot mean dFF---------
        plot(x1,mean(ROIdFF));
        xlim(xlim1);
        xt=[0:Stackpersec*100:550];
        xticks(xt);
        xticklabels({'0','100','200','300','400'});%Change this label!!

        yyaxis right
        plot(x1,Stim,'r');
%         title('Mean dFF of ROI cells');
        
        savefig(F1,strcat(wholeROIf,'-MeandFF.fig'));
    elseif i==3
        %-------Heat Map Part----------
        imagesc(ROIdFF);
        hold on
        colormap jet
        cl=caxis;
        caxis([-1 4]);
        xlim(xlim1);
        xt=[0:Stackpersec*100:550];
        xt=[0:0];
        xticks(xt);
        xticklabels({});%Change this label!!
        xticklabels({'0','100','200','300','400'});%Change this label!!
%         xlabel('Frame Count');
        ylabel('dF/F');
        ylim([-0.1 0.4]);
        ylabel('Cell #');
        ylim([0 size(ROIFluo,1)+1]);
        
        yyaxis right
        plot(x1,Stim,'r');
%         title('dFF of ROI cells');
        
        savefig(F1,strcat(wholeROIf,'-dFF.fig'));
%     elseif i==4
%         %-------Heat Map Part----------
%         imagesc(ROIFluo);
%         hold on
%         colormap jet
%         cl=caxis;
%         caxis([-1 4]);
%         xlim([0,500]);
%         xt=[0:Stackpersec*100:500];
%         xt=[0:0];
%         xticks(xt);
%         xticklabels({});%Change this label!!
%         xticklabels({'0','100','200','300','400'});%Change this label!!
%         xlabel('Frame Count');
%         ylabel('dF/F');
%         ylim([-0.1 0.4]);
%         ylabel('Cell #');
%         ylim([0 size(ROIFluo,1)]);
%         
%         yyaxis right
%         plot(x1,Stim,'r');
%         title('dFF of ROI cells');
    end
end

print='Done!'
% ---------------Same plot but based on time---------
% pos1=[0.1 0.1 0.8 0.2];
% F1=figure;
% 
% time=round(xlength/1.28);
% x2=(1:time);
% for i=1:time
% %     if Stim(round(time*1.28))>0
%         Stim2(i)=Stim(floor(i*1.28));
%         dFF2(i)=dFF(floor(i*1.28));
% %     end
% end
% 
% subplot('Position',pos1);
% plot(x2,Stim2);
% xlim([0,500]);
% xlabel('Time(s)');
% ylabel('Stim Amplitude');
% 
% pos2=[0.1 0.4 0.8 0.5]
% subplot('Position',pos2);
% plot(x2,dFF2);
% xlim([0,500]);
% % xlabel('Frame Count');
% ylabel('dF/F');
% print='Done!'
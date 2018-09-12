load('J:\Cui\ProcessedData 060518\ForeBrainCrop\20180605_151325-F1T1-6f-TS1.mat-Cropped.mat');
dFF=dFF;
load('J:\Cui\ProcessedData 060518\StimTimeSeries(use this one!) 060518\StimTimeSeriesOf20180605_151325-F1T1-6f-TS1.mat')
Stim=StimTimeSeries;

F1=figure;
% left_color = [0 0 0];
% right_color = [1 0 0];
% set(F1,'defaultAxesColorOrder',[left_color; right_color]);

% imagesc(dFF);
% colormap jet
% cl=caxis;
% caxis([-1 4]);
% 
% CB=colorbar;
% CB.Label.String='dF/F';

xt=[0:Stackpersec*100:650];
xticks(xt);
xticklabels({'0','100','200','300','400','500'});%Change this label!!
xlabel('Time(s)');
ylabel('Cell #');

hold on 

% yyaxis right
x1=(1:size(dFF,2));
plot(x1,Stim);
xlim([0 650]);
% ylim([0,3.5]);
ylabel('Stim Amplitude');

% title('Heatmap of Forebrain Activity During Stimulation');


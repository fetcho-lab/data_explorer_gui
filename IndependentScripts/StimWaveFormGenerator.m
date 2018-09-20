%Waveform generator
clear all
[Curfile,Curpath] = uigetfile('.txt','Please select stimulus file');
cd(Curpath);
AllCurFileName= [Curpath,Curfile]
% [StimPattern,StimIntv]=textread(AllCurFileName,'%f');
StimFile=textread(AllCurFileName,'%f');
StimPattern=StimFile(1:end-1);
StimInterval=StimFile(length(StimFile));

% answer=inputdlg('Total Time Point No.:','Input Data'); 
% TimePointNo=str2num(answer{1,1});

[Datafile,Datapath] = uigetfile('.mat','Please select data file');
cd(Datapath);
AllCurDataFileName= [Datapath,Datafile]
CurrentData=load(AllCurDataFileName);
TimePointNo=size(CurrentData.fluorescence_time_series,2)

StimTimeSeries=zeros(1,TimePointNo);

for i=1:length(StimPattern)
    if StimPattern(i)>0
%         StimTimeSeries(floor(i*StimInterval*0.8))=StimPattern(i);%Use
%         this if the frame per stack was wrongly set to be 36 instead of
%         45.
        StimTimeSeries(i*StimInterval)=StimPattern(i);
        StimTimeSeries(i*StimInterval+1)=StimPattern(i);
    end
end

fig1=plot(1:length(StimTimeSeries),StimTimeSeries);
title=strcat('StimTimeSeriesOf',Datafile);

save(strcat('StimTimeSeriesOf',Datafile),'StimTimeSeries');
savefig(strcat(Curpath,'StimTimeSeriesOf',Datafile,'.fig'));
%Move Stack 1 frame back

clear all
[Curfile,Curpath] = uigetfile('.mat','Please select stimulus file');
cd(Curpath);
AllCurFileName= [Curpath,Curfile]

load(AllCurFileName);

StimTimeSeries2=StimTimeSeries;
StimTimeSeries=zeros(1,length(StimTimeSeries2));

for i=1:length(StimTimeSeries2)
    StimTimeSeries(i+2)=StimTimeSeries2(i);
end


save(strcat('Amended',Curfile(1:end-4)),'StimTimeSeries');


fig1=plot(1:length(StimTimeSeries),StimTimeSeries);
title=strcat('Amended',Curfile(1:end-4));

savefig(strcat(Curpath,'Amended',Curfile(1:end-4),'.fig'));

print1='Done!'
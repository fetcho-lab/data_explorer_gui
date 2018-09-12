%Batch Cell Merger
clear all
TargetCellNo=10000;%change target cell no here!!!

list_of_directories = {...
                'C:\CUI\ProcessedData 062118\Data-dFFComputed-Merged\t',...
    };

for directory_idx  = 1:numel(list_of_directories)
    cd(list_of_directories{directory_idx});
    
    disp(sprintf('Converting %s',list_of_directories{directory_idx}));
    filelist=dir('*.mat');
    for i=1:length(filelist)
        curdata=load(filelist(i).name);
        [fluorescence_time_series,spPos1]=CellMerger(curdata.fluorescence_time_series,curdata.spPos,TargetCellNo,5);
        [dFF,spPos]=CellMerger(curdata.dFF,curdata.spPos,TargetCellNo,5);
%         mkdir MergedData
        savefilename=strcat(num2str(TargetCellNo),'Merged-',filelist(i).name)
        
%         fluorescence_time_series=curdata.fluorescence_time_series;
%         dFF=curdata.dFF;
%         spPos=curdata.spPos;
        Sc=curdata.Sc;
        cellSegmentation=curdata.cellSegmentation;
        extractParams=curdata.extractParams;
        spRadiiXYZ=curdata.spRadiiXYZ;
        save (savefilename,'fluorescence_time_series','dFF','spPos','Sc','cellSegmentation','extractParams','spRadiiXYZ');
        display(strcat(num2str(i),'/',num2str(length(filelist)),'in currentfolder Done!'));
    end
    display(strcat(num2str(directory_idx),'/',num2str(numel(list_of_directories)),'foldr is Done!'));
end
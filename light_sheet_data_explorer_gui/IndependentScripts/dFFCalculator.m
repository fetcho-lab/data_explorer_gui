function [AlldFF,basalF,ZeroBasalFCells_discarded]=dFFCalculator(Allfluo,BasalFVolume)
% if Allfluo is a nxm matrix, n should be the number of cells, m should be
% number of timepoints.
% BasalFVolume should be the number of data points (among all timepoints'
% data) you want to be taken as basal F.
    TimepointNo=size(Allfluo,2);
    sortedfluo=sort(Allfluo,2,'ascend');
    basalF= mean(sortedfluo(:,1:BasalFVolume),2);
    ZeroBasalFCells_discarded=find(basalF==0);
    if size(ZeroBasalFCells_discarded)>0
        warndlg(strcat(num2str(size(ZeroBasalFCells_discarded)),' basal F values contains 0, try to make Basal Volume higher? Cells whose basal F is 0 will be deleted.'));
    end
    basalFexpand=basalF*ones(1,TimepointNo);
    deltaF=Allfluo-(basalFexpand);
    AlldFF=deltaF./basalFexpand;
    
    AlldFF(ZeroBasalFCells_discarded,:)=[];
end


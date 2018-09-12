function [AllMerFluo, AllMerPos]=CellMerger(Allfluo, Allpos, TargetCellNo, ZStepsize)
%This merging function will (likely?) randomly discard some data points, maybe especially cells in the periphral.
    [OriginalCellNo, TimepointNo]=size(Allfluo);
    allz=Allpos(:,3);
    ZLayerNo=((max(allz)-min(allz))/ZStepsize)+1;
    ratio=ceil(OriginalCellNo/TargetCellNo);
    
    MergedCellNo=1
    for i=1:ZLayerNo
        CurrentLayer=(i-1)*ZStepsize+min(allz);
        CurrentCellList=find(allz==CurrentLayer);
        CurrentCellPos=Allpos(CurrentCellList,:);
        for i2=1:floor(length(CurrentCellList)/ratio)
            KernelCellNo=CurrentCellList(1);
            KernelCellPos=CurrentCellPos(1,:);
            DistanceToOthers=pdist2(KernelCellPos,CurrentCellPos);
            [SortedDis,SortedDisCellList]=sort(DistanceToOthers);
            CellNoToMerge=CurrentCellList(SortedDisCellList(1:ratio));
            
            AllMerPos(MergedCellNo,:)=mean(Allpos(CellNoToMerge,:));
            AllMerFluo(MergedCellNo,:)=mean(Allfluo(CellNoToMerge,:));
            MergedCellNo=MergedCellNo+1;
            
            CurrentCellList(SortedDisCellList(1:ratio),:)=[];
            CurrentCellPos(SortedDisCellList(1:ratio),:)=[];
        end
        MergedCellNo
    end
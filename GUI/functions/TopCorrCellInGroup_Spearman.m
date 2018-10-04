function [Allcor,AllPVal]=TopCorrCellInGroup_Spearman(AllFluoOrAlldFF)
    CellNo=size(AllFluoOrAlldFF,1);
    [Allcor, AllPVal]=corr(AllFluoOrAlldFF','Type','Spearman');
    for i=1:CellNo
        Allcor(i,i)=0;
    end
%         [Maxcor, Maxcorno]=max(Allcor);
end
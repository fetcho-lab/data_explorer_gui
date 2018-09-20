function [Allcor,AllPVal]=TopCorrCellInGroup_Pearson(AllFluoOrAlldFF)
    CellNo=size(AllFluoOrAlldFF,1);
    [Allcor, AllPVal]=corr(AllFluoOrAlldFF','Type','Pearson');
    for i=1:CellNo
        Allcor(i,i)=0;
    end
%         [Maxcor, Maxcorno]=max(Allcor);
end
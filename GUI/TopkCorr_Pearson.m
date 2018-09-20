function [TopCorVal, TopCorNo, AllCor, AllPVal,LowdFFList] = TopkCorr_Pearson(AllFluoOrAlldFF, inputSeries, inputCellNo, k,FilterThreshold)
%AllCor(i,j) is the correlation between i th input series and j th cell in
%all fluo series.
    AllCellNo=size(AllFluoOrAlldFF,1);
    TimePointNo=size(AllFluoOrAlldFF,2);
    inputSeriesNo=size(inputSeries,1);
    
%     AllCor=zeros(inputSeriesNo, AllCellNo);
    
       [AllCor, AllPVal]=corr(AllFluoOrAlldFF',inputSeries','Type','Pearson');
       [Nanx,Nany]=find(isnan(AllCor)==1);
       for i2=1:length(Nanx)
            AllCor(Nanx(i2),Nany(i2))=0;
       end
    
    assignin('base','AllCor',AllCor);
    assignin('base','AllPVal',AllPVal);
    [LowdFFList]=LowMaxFluoOrdFFDetector(AllFluoOrAlldFF,FilterThreshold);
    assignin('base','LowdFFList',LowdFFList);
	AllCor(LowdFFList,:)=0;
    AllPVal(LowdFFList,:)=NaN;
    
    assignin('base','AllCor2',AllCor);
    assignin('base','AllPVal2',AllPVal);
%            AllCor(i,inputCellNo(i))=0;%Make each input cell's self correlation to 0 so they can be excluded in max searching result.
    
%---------------Sorting part----------------
    if k==1
        [TopCorVal, TopCorNo]=max(AllCor,[],1);
    else
        [SortedAllCor,SortedAllCorNo]=sort(AllCor,1, 'descend');
        TopCorVal=SortedAllCor(1:k,:);
        TopCorNo=SortedAllCorNo(1:k,:);
        assignin('base','TopCorVal',TopCorVal);
    end
end

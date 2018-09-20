function [ MaxCor, MaxCorNo] = CalculateCorrcoef(AllFluoOrAlldFF, inputSeries, inputCellNo, ResultNo)
    AllCellNo=size(AllFluoOrAlldFF,1);
    TimePointNo=size(AllFluoOrAlldFF,2);
    inputSeriesNo=size(inputSeries,1);
    AllCor=zeros(inputSeriesNo, AllCellNo);
    for i=1:inputSeriesNo
           for j=1:AllCellNo
               TempCor=corrcoef(AllFluoOrAlldFF(j,:),inputSeries(i,:));
               AllCor(i,j)=TempCor(1,2);
           end
           AllCor(i,inputCellNo(i))=0;%Make each input cell's self correlation to 0 so they can be excluded in max searching result.
    end
    if ResultNo==1
        [MaxCor, MaxCorNo]=max(Allcor,[],2);
    else
        [SortedAllCor,SortedAllCorNo]=sort(AllCor, 'descend');
        MaxCor=SortedAllCor(1:ResultNo,:);
        MaxCorNo=SortedAllCorNo(1:ResultNo,:);
    end
    %------------AllFluo version end------------
end

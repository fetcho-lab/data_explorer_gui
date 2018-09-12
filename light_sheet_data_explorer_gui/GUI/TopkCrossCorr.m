function[TopCorVal, TopCorNo, AllCor, AllLag] = TopkCrossCorr(FluoOrdFF, inputSeries, TargetSeriesNo, k, LagRange);

    AllCellNo=size(FluoOrdFF,1);
    TimePointNo=size(FluoOrdFF,2);
    inputSeriesNo=size(inputSeries,1);
    
    AllCor=zeros(inputSeriesNo, AllCellNo);
    AllLag=zeros(inputSeriesNo, AllCellNo);
    
    for i=1:AllCellNo
        parfor j=1:inputSeriesNo
            [AllCor(i,j),AllLag(i,j)]=max(xcorr(FluoOrdFF(i,:),inputSeries(j,:),LagRange));
            AllLag=AllLag-(LagRange+1);
        end
        if mod(i,1000)==0
            i
        end
    end
    
    if k==1
        [TopCorVal, TopCorNo]=max(AllCor,[],1);
    else
        [SortedAllCor,SortedAllCorNo]=sort(AllCor,1, 'descend');
        TopCorVal=SortedAllCor(1:k,:);
        TopCorNo=SortedAllCorNo(1:k,:);
    end
end
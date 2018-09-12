function [ListToDiscard]=LowMaxFluoOrdFFDetector(FluoOrdFF,threshold)
    maxoffluo=max(FluoOrdFF,[],2);
    ListToDiscard=find(maxoffluo<threshold);
end
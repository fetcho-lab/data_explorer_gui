function [ output_args ] = PlotCellMap( Allpos,use_dFF, AllfluoOrAlldFF,SizeRangeForDots,TransparencyForDots)
%Sequence of Input Arguments: (Allpos, use_dFF, AllfluoOrAlldFF,SizeRangeForDots, TransparencyForDots). You can set use_dFF to 0 if you want to use fluo for cell map.
%   Sequence of Input Arguments: (Allpos, use_dFF, AllfluoOrAlldFF,SizeRangeForDots, TransparencyForDots). You can set use_dFF to 0 if you want to use fluo for cell map.
if use_dFF==1
    MaxdFFOfCells=max(AllfluoOrAlldFF');
    ColorOfDots=MaxdFFOfCells;
    TotalDifOfMaxdFF=max(MaxdFFOfCells)-min(MaxdFFOfCells);
    SizeOfDots=((MaxdFFOfCells-min(MaxdFFOfCells))./TotalDifOfMaxdFF)*(SizeRangeForDots(2)-SizeRangeForDots(1))+SizeRangeForDots(1);
else
    MaxFluoOfCells=max(AllfluoOrAlldFF');
    ColorOfDots=MaxFluoOfCells;
    TotalDifOfMaxFluo=max(MaxFluoOfCells)-min(MaxFluoOfCells);
    SizeOfDots=((MaxFluoOfCells-min(MaxFluoOfCells))./TotalDifOfMaxFluo)*(SizeRangeForDots(2)-SizeRangeForDots(1))+SizeRangeForDots(1);
end

Allposx=Allpos(:,1);
Allposy=Allpos(:,2);
Allposz=Allpos(:,3);
%-----------To Be Adjusted--------------
scatterplot1=scatter3(Allposx,Allposy,Allposz,SizeOfDots,...
    ColorOfDots,'MarkerFaceAlpha',TransparencyForDots,'MarkerEdgeAlpha',TransparencyForDots)%Adjust this line for different plotting settings!%
%--------------------------------------
xlabel('X');
ylabel('Y');
zlabel('Z');
colormap jet;
CB1=colorbar;
if use_dFF==1
    CB1.Label.String='dF/F of cells (dots)';
else
    CB1.Label.String='Fluorescence intensity of cells (dots)';
end
FluoColorbarLim1=caxis;

if use_dFF==1
    save('FilteredCellMapPlotParameters_dFF.mat','use_dFF','FluoColorbarLim1','SizeRangeForDots','TransparencyForDots');
    savefig('FilteredCellMap_dFF.fig');
else
    save('FilteredCellMapPlotParameters_Fluo.mat','use_dFF','FluoColorbarLim1','SizeRangeForDots','TransparencyForDots');
    savefig('FilteredCellMap_Fluo.fig');
end
hold off

end


function genFIJI_ROIs_CuiEdited(SpotPositions, SpotIDs, radii)
%genFIJI_ROIs(SpotPositions, SpotIds, radii)
%   Generates a .txt file for each Z and creates the information necessary
%   to map ROIs using the drawROIs.py script. SpotIDs is a column vector
%   and radii is a 1x3 vector. 

pixelConversion = diag([0.41,0.41,5]);
pxPositions = SpotPositions*pixelConversion^-1;
pxRadii = radii*pixelConversion^-1;

zRange = [min(pxPositions(:,3)) max(pxPositions(:,3))];
zEdges = [0:ceil(zRange(2))];
zCenters = zEdges+0.5;
[zHist,zSteps] = histc(pxPositions(:,3),[0:ceil(zRange(2))]);

mkdir('spotROIs');

for Z=1:numel(zHist)-1
    
    zSpotsIdx = abs(pxPositions(:,3) - zCenters(Z)) <= pxRadii(3);
    
    if sum(zSpotsIdx)
        zOffsets = pxPositions(zSpotsIdx,3) - zCenters(Z);
        offSetRadii = pxRadii(1)*sqrt(1 - zOffsets.^2/pxRadii(3)^2);

        writeList = [SpotIDs(zSpotsIdx) pxPositions(zSpotsIdx,:) offSetRadii];

        writeName = sprintf('spotROIs/ROILocations_Plane%02.f',zEdges(Z)+1);
        writeFile = fopen(writeName,'w');
        fprintf(writeFile,'%4.0f\t %f\t %f\t %f\t %f\n',writeList');
        fclose(writeFile);

    end

end
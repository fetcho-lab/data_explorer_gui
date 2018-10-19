function lookup_table_pixels = calculate_exact_cell_radii_by_plane(SpotPositions, radii, Sc)
%lookup_table = calculate_exact_cell_radii_by_plane(SpotPositions, radii)
%   Calulates the radius of each spot in increasing z-planes based on its
%   center and x,y,z radii (x and y are assumed to be equal). Results are
%   returned in pixels.

pixelConversion = Sc;
pxPositions = SpotPositions*pixelConversion^-1;
pxRadii = radii*pixelConversion^-1;

zRange = [min(pxPositions(:,3)) max(pxPositions(:,3))];
zEdges = [0:ceil(zRange(2))];
zCenters = zEdges+0.5;

lookup_table_pixels = zeros(size(SpotPositions,1), length(zCenters));

for Z=1:length(zCenters)
    
    zSpotsIdx = abs(pxPositions(:,3) - zCenters(Z)) <= pxRadii(:,3);
    
    if sum(zSpotsIdx)
        zOffsets = pxPositions(zSpotsIdx,3) - zCenters(Z);
        offSetRadii = pxRadii(zSpotsIdx,1).*sqrt(1 - zOffsets.^2./pxRadii(zSpotsIdx,3).^2);

        lookup_table_pixels(zSpotsIdx,Z) = offSetRadii; %writes the projected radii onto each slice 1:Z
    end

end
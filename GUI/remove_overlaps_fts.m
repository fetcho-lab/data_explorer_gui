function [ toKeep, toRemove ] = remove_overlaps_fts( intensity, spPos_um, micron_threshold, adjSliceOnly )
% [ toRemove ] = remove_overlaps_fts( intensity, spPos_um, micron_threshold, adjSliceOnly )
%   This function removes any neurons that are within the threshold
%   distance of each other by only keeping the brightest of the pair in
%   intensity. It returns a logical of which neurons should be removed. 

toRemove = zeros(size(spPos_um,1), 1, 'logical');
% adjSliceOnly = 1;
disp('Removing overlaps...');
tic
m=0;
while m<length(toRemove) %a while loop is implemented to allow re-check f
    
    m = m+1; 
%     disp(m)
    if toRemove(m)
       continue 
    end
    
    manhattan_distance = sum( abs( bsxfun(@minus, spPos_um, spPos_um(m,:) ) ), 2);
    
    if adjSliceOnly
        inZ = abs( spPos_um(:,3) - spPos_um(m,3) ) < 5;
        within_threshold = [manhattan_distance <  sqrt(2)*micron_threshold] & ~ inZ;
    else
        within_threshold = manhattan_distance < 2*sqrt(2)*micron_threshold;
        within_threshold(m) = 0;
    end
    
    
    
    if sum(within_threshold & ~toRemove) > 0
        xCells = find(within_threshold & ~toRemove);
        exact_distance = pdist2(spPos_um(m,:), spPos_um(xCells,:));
        
        for k=1:length(exact_distance)
            
            if exact_distance(k) < micron_threshold
                cells_comparing = [m xCells(k)];
                [br, selected] = max( [intensity(m), intensity(xCells(k))] );
                
                toRemoveIdx = cells_comparing( mod(selected, 2) + 1);
                toRemove(toRemoveIdx) = 1;
                
                m = m - 1; 
                break;
            end
        end
        
    end
end
toc

toKeep = ~toRemove;
fprintf('Selected %4.0f cells for removal\n', sum(toRemove)) ;




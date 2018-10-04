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
while m<length(toRemove) %a while loop is implemented to allow re-check of the current neuron after each removal
    
    m = m+1; 
%     disp(m)
    if toRemove(m) %if current neuron has been flagged already, skip this one. 
       continue 
    end
    
    manhattan_distance = sum( abs( bsxfun(@minus, spPos_um, spPos_um(m,:) ) ), 2); %use manhattan distance to narrow down a region to check to save computational time. 
    
    if adjSliceOnly %ONLY check slices that aren't the same as the candidate. this deals strictly with z-overlaps.
        inZ = abs( spPos_um(:,3) - spPos_um(m,3) ) < 5;
        within_threshold = [manhattan_distance <  sqrt(2)*micron_threshold] & ~ inZ;
    else %otherwise check everywhere. 
        within_threshold = manhattan_distance < 2*sqrt(2)*micron_threshold;
        within_threshold(m) = 0;
    end
    
    
    
    if sum(within_threshold & ~toRemove) > 0 %if all candidates haven't been removed...
        xCells = find(within_threshold & ~toRemove);
        exact_distance = pdist2(spPos_um(m,:), spPos_um(xCells,:));%now, check the exact distances to find those below threshold
        
        for k=1:length(exact_distance)
            
            if exact_distance(k) < micron_threshold
                cells_comparing = [m xCells(k)]; %sets up a structure to easily access the brightest of the two being compared below threshold
                [br, selected] = max( [intensity(m), intensity(xCells(k))] ); 
                
                toRemoveIdx = cells_comparing( mod(selected, 2) + 1); %selects the dimmer one for removal. maybe min would be more straightforward above. 
                toRemove(toRemoveIdx) = 1;
                
                m = m - 1; %if something is removed, we must recheck the current neuron to see if any other candidates are also below threshold. 
                break;
            end
        end
        
    end
end
toc

toKeep = ~toRemove;
fprintf('Selected %4.0f cells for removal\n', sum(toRemove)) ;




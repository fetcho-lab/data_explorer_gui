function [ toRemove, toKeep ] = remove_overlaps_fts( intensity, spPos_um, micron_threshold )
% [ toRemove ] = remove_overlaps_fts( intensity, spPos_um, micron_threshold )
%   This function removes any neurons that are within the threshold
%   distance of each other by only keeping the brightest of the pair in
%   intensity. It returns a logical of which neurons should be removed. 

toRemove = zeros(size(spPos_um,1), 1, 'logical');

disp('Removing overlaps...');
tic
m=0;
while m<length(toRemove) %a while loop is implemented to allow re-check of candidate neurons in case more than 1 is within the threshold distance
    
    m = m+1; 
    
    if toRemove(m) %if the current candidate has already been flagged, continue to the next
       continue 
    end
    
    %manhattan distance is used to speed up computation time. Euclidean is used for those that fall within a rough distance of the candidate
    manhattan_distance = sum( abs( bsxfun(@minus, spPos_um, spPos_um(m,:) ) ), 2); 
    within_threshold = manhattan_distance < 2*sqrt(2)*micron_threshold;
    within_threshold(m) = 0; %the candidate itself will of course be within threshold and should be excluded
    
    if sum(within_threshold) > 0
        xCells = find(within_threshold & ~toRemove); %the exact distance is calculated for a smaller number of cells. xCells is the row number in spPos_um. those that already have been removed are excluded
        exact_distance = pdist2(spPos_um(m,:), spPos_um(xCells,:));
        
        for k=1:length(exact_distance)
            
            if exact_distance(k) < micron_threshold
                cells_comparing = [m xCells(k)]; %this bit of code selects the brightest pair for the first pair that is under the distance threshold
                [br, selected] = max( [intensity(m), intensity(xCells(k))] );
                
                toRemoveIdx = cells_comparing( mod(selected, 2) + 1);
                toRemove(toRemoveIdx) = 1; %this flags that particular neuron for removal. 
                
                m = m - 1; %re-check the candidate to eat up all pairs that are too close. 
                break;
            end
        end
        
    end
end
toc

toKeep = ~toRemove;
fprintf('Selected %4.0f cells for removal\n', sum(toRemove)) ;




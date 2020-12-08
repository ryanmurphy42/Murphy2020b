function [possible, res_points] = function_DISCRETE_Voronoi_onerespointpossible(x_current)

%% Determine whether it is possible to define an Interval-Voronoi partition using only a residence point in each cell

% Inputs:
% x_current = cell boundaries
% N_cells = number of cells


% Calculate the cell lengths
cell_lengths = diff(x_current);
N_cells = length(cell_lengths);

% create vector of possible regions in each cell
vor_pos_region_min = zeros(N_cells, 1);
vor_pos_region_max = zeros(N_cells, 1);
vor_pos_region_size = zeros(N_cells, 1);
res_points = zeros(N_cells, 1); % store the resident points of the cells

% First cell
vor_pos_region_min(1) = x_current(1);
vor_pos_region_max(1) = x_current(2);
vor_pos_region_size(1) = vor_pos_region_max(1)-vor_pos_region_min(1);


% set a variable for possible - 0 - not possible, 1 - is possible
possible = 1;

% For second cell to final cell
for jj= 2:N_cells
    
    if possible == 1
        if vor_pos_region_size(jj-1) < 1e-10
            possible =0;
        end
        
        previous_vor_pos_region_max_reflect = (x_current(jj) - vor_pos_region_max(jj-1)) + x_current(jj); % this is the left boundary in the current cell - due to reflection
        previous_vor_pos_region_min_reflect = (x_current(jj) - vor_pos_region_min(jj-1)) + x_current(jj);  % this is the right boundary in the current cell - due to reflection
        
        %Intersect this with the cell
        
        
        if (previous_vor_pos_region_max_reflect <= x_current(jj+1)) && (previous_vor_pos_region_min_reflect <= x_current(jj+1))
            % Possibility 1 - the reflected region is entirely inside the next cell
            vor_pos_region_min(jj) = previous_vor_pos_region_max_reflect;
            vor_pos_region_max(jj) = previous_vor_pos_region_min_reflect;
            vor_pos_region_size(jj) = vor_pos_region_max(jj)-vor_pos_region_min(jj);
            
        elseif (previous_vor_pos_region_max_reflect <= x_current(jj+1)) && (previous_vor_pos_region_min_reflect > x_current(jj+1))
            % Possibility 2 - the reflected region one edge is inside the cell and one edge outside
            vor_pos_region_min(jj) = previous_vor_pos_region_max_reflect;
            vor_pos_region_max(jj) = min(x_current(jj+1),previous_vor_pos_region_min_reflect);
            vor_pos_region_size(jj) = vor_pos_region_max(jj)-vor_pos_region_min(jj);
            
        elseif (previous_vor_pos_region_max_reflect > x_current(jj+1)) %&& (previous_vor_pos_region_min_reflect > x_current(jj+1))
            % Possibility 3 - the reflected region both edges are outside the next cell
            vor_pos_region_min(jj) =  min(x_current(jj+1),previous_vor_pos_region_max_reflect);
            vor_pos_region_max(jj) =  min(x_current(jj+1),previous_vor_pos_region_min_reflect);
            vor_pos_region_size(jj) = vor_pos_region_max(jj)-vor_pos_region_min(jj);
            
        end
        
        
    end
    
    
    
end

% also check possible based on final cell
if vor_pos_region_size(jj) < 1e-10
    possible =0;
end

if possible == 1
    %when it is possible to define the voronoi partition for one compartment per cell
    %work backwards to calculate the Voronoi partition
    
    % choose the centre of the possible region for the final cell
    res_points(end) = 0.5*(vor_pos_region_min(end) + vor_pos_region_max(end));
    
    % work backwards through the cells determining the resident points
    for jj=1:(N_cells-1)
        res_points(N_cells-jj) = x_current(N_cells+1-jj) - (res_points(N_cells +1 - jj)- x_current(N_cells+1-jj));
    end
    
else
    % it is not possible - will have to subdivide the cells  - in different script
    res_points = 0.5*(x_current(1:end-1) + x_current(2:end));
%     fprintf('\n Error: Voronoi \n');
end


end
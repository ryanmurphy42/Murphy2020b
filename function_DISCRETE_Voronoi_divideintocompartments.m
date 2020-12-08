function  [res_points,interval_lengths_compartments,compartments_per_cell,possible] = function_DISCRETE_Voronoi_divideintocompartments(x_current)

%% Determine whether it is possible to define an Interval-Voronoi partition using only a residence point in each cell

% Inputs:
% x_current = cell boundaries
% N_cells = number of cells


x_compartment = x_current;

% Calculate the cell lengths
cell_lengths = diff(x_current);
N_cells = length(cell_lengths);
N_compartments = length(cell_lengths); % initially set this equal to the cell boundaries - updated below


compartments_per_cell = ones(length(cell_lengths),1); % intially set the number of compartments per cell equal to one

% create vector of possible regions in each cell
vor_pos_region_min = zeros(N_compartments, 1);
vor_pos_region_max = zeros(N_compartments, 1);
vor_pos_region_size = zeros(N_compartments, 1);
res_points = zeros(N_compartments, 1); % store the resident points of the cells


% First cell
vor_pos_region_min(1) = x_compartment(1);
vor_pos_region_max(1) = x_compartment(2);
vor_pos_region_size(1) = vor_pos_region_max(1)-vor_pos_region_min(1);


% set a variable for possible - 0 - not possible, 1 - is possible
possible = 1;

cell_counter =1;
interval_counter=1;

% For second cell to final cell
while cell_counter < N_cells
    cell_counter = cell_counter + 1;
    interval_counter=interval_counter+1;
    if possible == 1
        if vor_pos_region_size(interval_counter-1) < 1e-10
            possible =0;
        end
        
        previous_vor_pos_region_max_reflect = (x_compartment(interval_counter) - vor_pos_region_max(interval_counter-1)) + x_compartment(interval_counter); % this is the left boundary in the current cell - due to reflection
        previous_vor_pos_region_min_reflect = (x_compartment(interval_counter) - vor_pos_region_min(interval_counter-1)) + x_compartment(interval_counter);  % this is the right boundary in the current cell - due to reflection
        
        %Intersect this with the cell
        
        
        if (previous_vor_pos_region_max_reflect <= x_compartment(interval_counter+1)) && (previous_vor_pos_region_min_reflect <= x_compartment(interval_counter+1))
            % Possibility 1 - the reflected region is entirely inside the next cell
            vor_pos_region_min(interval_counter) = previous_vor_pos_region_max_reflect;
            vor_pos_region_max(interval_counter) = previous_vor_pos_region_min_reflect;
            vor_pos_region_size(interval_counter) = vor_pos_region_max(interval_counter)-vor_pos_region_min(interval_counter);
            
            
        elseif (previous_vor_pos_region_max_reflect <= x_compartment(interval_counter+1)) && (previous_vor_pos_region_min_reflect > x_compartment(interval_counter+1))
            % Possibility 2 - the reflected region one edge is inside the cell and one edge outside
            vor_pos_region_min(interval_counter) = previous_vor_pos_region_max_reflect;
            vor_pos_region_max(interval_counter) = min(x_compartment(interval_counter+1),previous_vor_pos_region_min_reflect);
            vor_pos_region_size(interval_counter) = vor_pos_region_max(interval_counter)-vor_pos_region_min(interval_counter);
            
        elseif (previous_vor_pos_region_max_reflect > x_compartment(interval_counter+1)) %&& (previous_vor_pos_region_min_reflect > x_compartment(interval_counter+1))
            % Possibility 3 - the reflected region both edges are outside the current cell
            
            %%      Divide the previous cell into two compartments
            
            % Obtain the previous cells possible region
            prev_cell_pos_region_min =  vor_pos_region_min(interval_counter-1);
            prev_cell_pos_region_max =  vor_pos_region_max(interval_counter-1);
            
            % Choose the location of the compartment boundary in the cell, x_compartment_new_temp
            
            % First choose x_compartment_new_temp so that reflecting the left boundary of the existing possible region for the previous cell equals the right boundary of the cell
            
            x_compartment_new_temp = 0.5*(prev_cell_pos_region_min + x_compartment(interval_counter));
            
            % Second update x_compartment_new_temp so that it does not overlap with the previous possible region of the cell
            if x_compartment_new_temp < prev_cell_pos_region_max
                x_compartment_new_temp = prev_cell_pos_region_max;
            end
            
            
            
            
            
            %%  update the possible region of second compartment of the previous cell
            
            vor_pos_region_min(interval_counter) =  x_compartment_new_temp + (x_compartment_new_temp -prev_cell_pos_region_max);
            vor_pos_region_max(interval_counter) =  x_compartment(interval_counter);
            vor_pos_region_size(interval_counter) = vor_pos_region_max(interval_counter)-vor_pos_region_min(interval_counter);
            
            % update x_compartment by introducing a compartment boundary at the midpoint of the previous cell
            x_compartment =     [x_compartment(1:(interval_counter-1));...
                x_compartment_new_temp;...
                x_compartment((interval_counter):end)];
            
            % update compartments per cell vector
            compartments_per_cell(cell_counter-1) = 2;
            
            
            %% Update the possible region of the current cell
            
            interval_counter = interval_counter + 1;
            N_compartments = N_compartments + 1;
            % Now only two possibilities
            
            previous_vor_pos_region_max_reflect = (x_compartment(interval_counter) - vor_pos_region_max(interval_counter-1)) + x_compartment(interval_counter); % this is the left boundary in the current cell - due to reflection
            previous_vor_pos_region_min_reflect = (x_compartment(interval_counter) - vor_pos_region_min(interval_counter-1)) + x_compartment(interval_counter);  % this is the right boundary in the current cell - due to reflection
            
            %Intersect this with the cell
            
            
            if (previous_vor_pos_region_max_reflect <= x_compartment(interval_counter+1)) && (previous_vor_pos_region_min_reflect <= x_compartment(interval_counter+1))
                % Possibility 1 - the reflected region is entirely inside the next cell
                vor_pos_region_min(interval_counter) = previous_vor_pos_region_max_reflect;
                vor_pos_region_max(interval_counter) = previous_vor_pos_region_min_reflect;
                vor_pos_region_size(interval_counter) = vor_pos_region_max(interval_counter)-vor_pos_region_min(interval_counter);
                
                
            elseif (previous_vor_pos_region_max_reflect <= x_compartment(interval_counter+1)) && (previous_vor_pos_region_min_reflect > x_compartment(interval_counter+1))
                % Possibility 2 - the reflected region one edge is inside the cell and one edge outside
                vor_pos_region_min(interval_counter) = previous_vor_pos_region_max_reflect;
                vor_pos_region_max(interval_counter) = min(x_compartment(interval_counter+1),previous_vor_pos_region_min_reflect);
                vor_pos_region_size(interval_counter) = vor_pos_region_max(interval_counter)-vor_pos_region_min(interval_counter);
                
                
            end
            
            
        end
        
        
    end
    
end

% also check possible based on final cell
if vor_pos_region_size(interval_counter) < 1e-10
    possible =0;
end


if possible == 1
    %when it is possible to define the voronoi partition for one compartment per cell
    %work backwards to calculate the Voronoi partition
    
    % choose the centre of the possible region for the final cell
    res_points(N_compartments) = 0.5*(vor_pos_region_min(end) + vor_pos_region_max(end));
    
    % work backwards through the cells determining the resident points
    for jj=1:(N_compartments-1)
        res_points(N_compartments-jj) = x_compartment(N_compartments+1-jj) - (res_points(N_compartments +1 - jj)- x_compartment(N_compartments+1-jj));
    end
    
else
    % it is not possible - will have to subdivide the cells  - in different script
    res_points = 0.5*(x_compartment(1:end-1) + x_compartment(2:end));
    fprintf('\n Error: Voronoi \n');
end


interval_lengths_compartments = diff(x_compartment); % initially set the compartment lengths equal to the cell lengths


end


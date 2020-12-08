function [res_points,interval_lengths_compartments,compartments_per_cell,possible] =  function_DISCRETE_Voronoiparition_of_interval(interval_lengths)

%% See if Voronoi partition possible current postions of cell boundaries.

x_current = [0;cumsum(interval_lengths)];

interval_lengths_compartments = interval_lengths;

compartments_per_cell = ones(length(interval_lengths),1);

[possible, res_points] = function_DISCRETE_Voronoi_onerespointpossible(x_current);

%% If a Voronoi partition is not possible with the current interval lengths then introduce compartments per cell.

if possible == 0
    
    [res_points,interval_lengths_compartments,compartments_per_cell,possible] = function_DISCRETE_Voronoi_divideintocompartments(x_current);
    
end

end

   
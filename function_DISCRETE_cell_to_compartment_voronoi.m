function [res_points,interval_lengths_compartments, compartments_per_cell_overall, interval_source_rates_compartments,onerespossible] = ...
    function_DISCRETE_cell_to_compartment_voronoi(interval_lengths, interval_source_rates)

%% Determine resident points using the Interval-Voronoi method.

compartments_per_cell_1 = ones(length(interval_lengths),1); % set compartments per cell equal to one to begin. This is then updated if necesssary.
[res_points,interval_lengths_compartments,compartments_per_cell_2,onerespossible] =  function_DISCRETE_Voronoiparition_of_interval(interval_lengths);

% if a Voronoi partition can not be defined, or if a voronoi partition can be defined but the distances between resident points is too far then update the resident points with compartments per cell

gamma=0.8;

max_loops =5;
interval_lengths_min = interval_lengths;
voronoi_loop_counter=1;
while (abs(log10(max(diff(res_points))) - log10(min(diff(res_points)))) > gamma) && (voronoi_loop_counter <= max_loops)
    
    
    max_length =0.5*min(interval_lengths_min); % set maximum length equal to half of minimum of current compartment lengths
    
	% divide cells into compartments so below threshold
    [interval_lengths_below_size_threshold, compartments_per_cell_1] =  function_DISCRETE_compartments_below_size_threshold(interval_lengths, max_length);
    
	% determine resident points`
    [res_points,interval_lengths_compartments,compartments_per_cell_2,onerespossible] =  function_DISCRETE_Voronoiparition_of_interval(interval_lengths_below_size_threshold);
    
	% update minimium length
    interval_lengths_min = interval_lengths_compartments;
    
	% increase loop counter so dont exceed maximum number of loops
    voronoi_loop_counter= voronoi_loop_counter + 1;
end

if voronoi_loop_counter > max_loops
   fprintf('Voronoi loop max reached without satisfying condition') 
end

%% Number of compartments per cell

% calculate the number of compartments per cell

compartments_per_cell_overall = zeros(length(compartments_per_cell_1),1);
compartments_per_cell_2_counter = 1;

for mm = 1:length(compartments_per_cell_1)
    running_compartment_total=0;
    for nn=1:compartments_per_cell_1(mm)
        running_compartment_total = running_compartment_total + compartments_per_cell_2(compartments_per_cell_2_counter);
        compartments_per_cell_2_counter = compartments_per_cell_2_counter + 1;
    end
    compartments_per_cell_overall(mm) = running_compartment_total;
end


%% Source rates

% update the source rate for the compartments

interval_source_rates_compartments = zeros(sum(compartments_per_cell_overall),1);

% only the final compartment in the final real cell is supplied with chemical
interval_source_rates_compartments(sum(compartments_per_cell_overall(1:end-1))) = interval_source_rates(end-1);

end
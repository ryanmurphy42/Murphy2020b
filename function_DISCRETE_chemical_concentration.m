function      [ chem_1_m] = function_DISCRETE_chemical_concentration(D_chem_1,chem_1_m,x_current_m,x_prev_m,cell_source_rates_chem_1,dt)

%% Update the chemical concentration for the cells

interval_lengths = diff(x_current_m); % calculate the cell lengths

%% Calculate the residence points and source/sink rates for chem1

[res_points,interval_lengths_compartments, compartments_per_cell_overall, compartment_source_rates_chem_1,onerespossible] = ...
    function_DISCRETE_cell_to_compartment_voronoi(interval_lengths, cell_source_rates_chem_1);



%% Check to see if the number of compartments changed

compartments_per_cell_overall_prev = compartments_per_cell_overall;

if max(compartments_per_cell_overall) > 1
    chem_1_m_compartments = zeros(sum(compartments_per_cell_overall_prev),1);
    chem_1_m_compartments_counter=1;
    for ii=1:length(compartments_per_cell_overall_prev)
        for jj = 1:compartments_per_cell_overall_prev(ii)
            chem_1_m_compartments(chem_1_m_compartments_counter) = chem_1_m(ii);
            chem_1_m_compartments_counter = chem_1_m_compartments_counter + 1;
        end
    end
else
    chem_1_m_compartments = chem_1_m;
end


%% Calculate the propensities - chemical 1

x_end_final_cell =  x_current_m(end-1);

[prop_for_nonum_chem_1, prop_bac_nonum_chem_1] = function_DISCRETE_calculate_propnonum(D_chem_1, res_points,x_end_final_cell);

%% Calculate concentration in each compartment chemical 1

x_current_m_cells = x_current_m; % to use to calculate x_prev_m.

%use the same number of compartments to calculate x_prev_m residence points to be used to calculate dilution term in chemical concentration.
if max(compartments_per_cell_overall) > 1
    
    x_current_m =  [0;0.5*(res_points(1:end-1) + res_points(1:end-1));res_points(end) + 0.5*(res_points(end)-res_points(end-1))] ; % calculate x_current_m based on the residence points
    
    %use the same number of compartments to calculate x_prev_m residence points and then update x_prev_m
    [x_prev_m_compartments] = function_DISCRETE_usenumcompartments_find_x_prev_m(x_current_m, x_prev_m, x_current_m_cells,compartments_per_cell_overall);
    
    x_prev_m =  x_prev_m_compartments;
end

[chem_1_m_compartments] = function_DISCRETE_lawmassactioneqns_chem(chem_1_m_compartments,prop_for_nonum_chem_1, prop_bac_nonum_chem_1, ...
    interval_lengths_compartments, compartment_source_rates_chem_1, dt,x_current_m,x_prev_m);

%% Return the concentration in each cell chemical 1

if max(compartments_per_cell_overall) > 1
    % Calculate the concentration of the cell given the concentration of compartments per cell
    [chem_1_m] = function_DISCRETE_concentration_compartment_to_cell(compartments_per_cell_overall, chem_1_m_compartments,interval_lengths_compartments);
else
    chem_1_m = chem_1_m_compartments;
end

end
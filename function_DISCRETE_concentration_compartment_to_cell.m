function [chem_1_m] = function_DISCRETE_concentration_compartment_to_cell(compartments_per_cell_overall, chem_1_m_compartments,interval_lengths_compartments)

%% Convert the concentration of a compartment to the concentration of a cell

len_compartments_per_cell_overall = length(compartments_per_cell_overall);

cell_lengths = zeros(len_compartments_per_cell_overall,1);
mass_per_cell = zeros(len_compartments_per_cell_overall,1);


comparment_counter=1;
for ii = 1:len_compartments_per_cell_overall

    for jj=1:compartments_per_cell_overall(ii)
        
        cell_lengths(ii) = cell_lengths(ii) + interval_lengths_compartments(comparment_counter);
        mass_per_cell(ii)  =       mass_per_cell(ii) + interval_lengths_compartments(comparment_counter)*chem_1_m_compartments(comparment_counter);
        comparment_counter = comparment_counter + 1;
    end

end

chem_1_m = mass_per_cell./cell_lengths;


end
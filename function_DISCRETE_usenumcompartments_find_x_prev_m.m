function [x_prev_m_compartments] = function_DISCRETE_usenumcompartments_find_x_prev_m(x_current_m, x_prev_m, x_current_m_cells,compartments_per_cell_overall)


% determine x_prev_m_compartments given the current positions of compartment boundaries and resident points

compartment_counter = 1;

diff_x_current_m = diff(x_current_m);
diff_x_current_m_cells = diff(x_current_m_cells);
diff_x_prev_m= diff(x_prev_m);

x_prev_m_compartments = zeros(length(x_current_m),1);

for ii = 1:length(compartments_per_cell_overall)-1
    ratio_vec = zeros(compartments_per_cell_overall(ii),1);
    for jj = 1:compartments_per_cell_overall(ii)
        ratio_vec(jj) = diff_x_current_m(compartment_counter)/diff_x_current_m_cells(ii);
        
        x_prev_m_compartments(compartment_counter+1) = x_prev_m_compartments(compartment_counter) +  ratio_vec(jj)*diff_x_prev_m(ii);

        compartment_counter = compartment_counter + 1;
    end
end

x_prev_m_compartments(end) = x_prev_m(end-1);

end


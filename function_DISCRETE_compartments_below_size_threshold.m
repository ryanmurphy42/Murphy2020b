function [interval_lengths_below_size_threshold, compartments_per_cell] = function_DISCRETE_compartments_below_size_threshold(interval_lengths, max_length)


%% Function to calculate how many compartments to divide the intervals into to be below the max length threshold

compartments_per_cell = zeros(length(interval_lengths),1);

for ii=1:(length(interval_lengths)-1)
    
    int_length = interval_lengths(ii);
    int_length_div_max = int_length/max_length;
    ceil_int_length_div_max = ceil(int_length_div_max);
    compartments_per_cell(ii) =  ceil_int_length_div_max;
    
    if ii == 1
        interval_lengths_below_size_threshold = (int_length/ceil_int_length_div_max)*ones(ceil_int_length_div_max,1);
    else
        interval_lengths_below_size_threshold = [interval_lengths_below_size_threshold;(int_length/ceil_int_length_div_max)*ones(ceil(int_length_div_max),1)];
    end
    
end

% final interval is not included as it is a fake cell
 ii = (length(interval_lengths));
    compartments_per_cell(ii) = 1;
    interval_lengths_below_size_threshold = [interval_lengths_below_size_threshold;interval_lengths(ii)];

end
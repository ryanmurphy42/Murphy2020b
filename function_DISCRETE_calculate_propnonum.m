function [prop_for_nonum, prop_bac_nonum] = function_DISCRETE_calculate_propnonum(D, res_points,x_end_final_cell)

%% Calculate forward (left) and backward (right) transition rates

num_cells = length(res_points);

prop_for_nonum = zeros(num_cells,1);
prop_bac_nonum = zeros(num_cells,1);

if num_cells > 2
    %% forward reactions
    
    prop_for_nonum(1) = (2*D)/( (res_points(2)-res_points(1))*( res_points(2) + res_points(1) )  );
    for ii=2:num_cells-2
        prop_for_nonum(ii) =   (2*D)/( (res_points(ii+1)-res_points(ii))*( res_points(ii+1) - res_points(ii-1) )  );
    end
    %final cell does not have a forward reaction
    
    
    %% backward reactions
    %first cell does not have a backward reaction
    for ii=2:num_cells-2
        prop_bac_nonum(ii) = (2*D)/( (res_points(ii)-res_points(ii-1))*( res_points(ii+1) - res_points(ii-1) )  );
    end
    %final cell
    prop_bac_nonum(num_cells-1) = (2*D)/( (res_points(num_cells-1) - res_points(num_cells-2))*(  2*(x_end_final_cell - res_points(num_cells-2)) - res_points(num_cells-1)  + res_points(num_cells-2)   )  );
    
    
end

prop_for_nonum(prop_for_nonum<0) =0;
prop_bac_nonum(prop_bac_nonum<0) =0;


end
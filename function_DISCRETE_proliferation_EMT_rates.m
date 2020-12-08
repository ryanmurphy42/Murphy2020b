function [G_m,G_prolif_m,G_EMT_m] = ....
    function_DISCRETE_proliferation_EMT_rates(N_cells,m_springs_per_cell, prolif_law,beta_prolif_m,x_current_m, EMT_law, EMT_param1, EMT_param2, chem_conc_finalcell,a_m)

%Calculate proliferation and death propensities 
% Prolif mechanism: 1 - cell-length-independent, 2 - linearly cell-length-dependent
% EMT mechanism: 1 - chemically-independent, 2 - chemically dependent

%% Proliferation law
G_prolif_m=zeros(N_cells*m_springs_per_cell,1);

if prolif_law == 1 %constant
    G_prolif_m =  beta_prolif_m;
elseif prolif_law == 2 %linear
    spring_sizes = diff(x_current_m);
    N_springs = length(spring_sizes);
    for ii=1:(N_springs-m_springs_per_cell) %for all cells until the final cell
        G_prolif_m(ii) = beta_prolif_m(ii)*m_springs_per_cell*spring_sizes(ii)/a_m(ii);
    end
end

%Final cell is not a real cell and cannot undergo proliferation
G_prolif_m(end)=0;

%% EMT law

G_EMT_m =zeros(N_cells*m_springs_per_cell,1);

if EMT_law == 1 %EMT at constant rate
    %EMT_param1 - constant rate
    if N_cells > 1
    G_EMT_m = EMT_param1;
    end
elseif EMT_law == 2 %EMT at chemical threshold
    %EMT param 2 - concentration threshold
    if N_cells > 1
       if  chem_conc_finalcell < EMT_param2 %if the concentration of the final cell is below the threshold then no chance of EMT
           G_EMT_m=0;
       else %if above the threshold then there is a rate of EMT
           G_EMT_m = EMT_param1;
       end
    end
    
end


%% Combination
G_m = [G_prolif_m; G_EMT_m];


end
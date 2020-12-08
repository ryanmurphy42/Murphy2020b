function [chem_1_m_tdt] = function_DISCRETE_lawmassactioneqns_chem(chem_1_m,prop_for_nonum, prop_bac_nonum,...
    interval_lengths, interval_source_rates,dt,x_current_m,x_prev_m)

%% Update the chemical concentration for each compartment.

% ode for the concentration of morphogen in each compartment, the uptake in each compartment i.e. the sink of the morphogen.
len_u=length(chem_1_m);

n=len_u;
chem_1_m_tdt = zeros(1,len_u);

%% Morphogen concentration

chem_1_m_tdt(1) = chem_1_m(1) + dt*( (1/interval_lengths(1))*(  prop_bac_nonum(2)*chem_1_m(2)*interval_lengths(2) - prop_for_nonum(1)*chem_1_m(1)*interval_lengths(1) )... %diffusion
    +interval_source_rates(1))... %source
    + ( -(chem_1_m(1)/(x_current_m(2)-x_current_m(1)))*( (x_current_m(2)-x_current_m(1)) - (x_prev_m(2)-x_prev_m(1))  ))... %dilution/concentrating
    ;

for i=2:n-2
    chem_1_m_tdt(i) = chem_1_m(i) + dt*((1/interval_lengths(i))*( prop_for_nonum(i-1)*interval_lengths(i-1)*chem_1_m(i-1)...
        -( prop_for_nonum(i) +  prop_bac_nonum(i) )*interval_lengths(i)*chem_1_m(i)...
        +prop_bac_nonum(i+1)*interval_lengths(i+1)*chem_1_m(i+1)  )... %diffusion
        +interval_source_rates(i))... %source)
        + ( -(chem_1_m(i)/(x_current_m(i+1)-x_current_m(i)))*( (x_current_m(i+1)-x_current_m(i)) - (x_prev_m(i+1)-x_prev_m(i))  ))... %dilution/concentrating
        ;
end

if n > 2
    chem_1_m_tdt(n-1) = chem_1_m(n-1) + dt*((1/interval_lengths(n-1))*( prop_for_nonum(n-2)*chem_1_m(n-2)*interval_lengths(n-2) - prop_bac_nonum(n-1)*chem_1_m(n-1)*interval_lengths(n-1))... %diffusion
        +interval_source_rates(n-1))... %source
        + ( -(chem_1_m(n-1)/(x_current_m(n)-x_current_m(n-1)))*( (x_current_m(n)-x_current_m(n-1)) - (x_prev_m(n)-x_prev_m(n-1))  ))... %dilution/concentrating
        ;
end

chem_1_m_tdt = chem_1_m_tdt';

end

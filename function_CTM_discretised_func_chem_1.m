function [g] = function_CTM_discretised_func_chem_1(q,k,a,nodesz,dt,dz,eta,L_old,L,prolif_law,beta_prolif,chem_1, chem_1_old, D_chem_1,node_velocity,omega,phi)

%% Chemical concentration, c

z=0:dz:1;

r = dt/(eta*dz^2*L^2);
dL_dt = (L-L_old)/dt;
g = zeros(nodesz,1);

%% First node
g(1) = chem_1(2) - chem_1(1); % x0 BC

%% Interior nodes
for i = 2:nodesz-1
    
    if node_velocity(i) >= 0 %upwinding
        
        g(i) = - chem_1(i) + chem_1_old(i) ...
            -(dt/dz)*node_velocity(i)*(chem_1(i)-chem_1(i-1))...
            -dt*chem_1(i)*(1/(2*eta*(L^2)*dz^2))*(  (1/q(i+1) + 1/q(i)   )*( k(i+1)*(1/q(i+1) - a(i+1)) - k(i)*(1/q(i) - a(i)) )...
            -(1/q(i) + 1/q(i-1)   )*( k(i)*(1/q(i) - a(i)) - k(i-1)*(1/q(i-1) - a(i-1)) )   )...
            +(dt/dz^2)*(D_chem_1/L^2)*( chem_1(i+1) - 2*chem_1(i) + chem_1(i-1))...
            ;
        
        
    elseif node_velocity(i) < 0 %downwinding
        
        g(i) = - chem_1(i) + chem_1_old(i) ...
            -(dt/dz)*node_velocity(i)*(chem_1(i+1)-chem_1(i))...
            -dt*chem_1(i)*(1/(2*eta*(L^2)*dz^2))*(  (1/q(i+1) + 1/q(i)   )*( k(i+1)*(1/q(i+1) - a(i+1)) - k(i)*(1/q(i) - a(i)) )...
            -(1/q(i) + 1/q(i-1)   )*( k(i)*(1/q(i) - a(i)) - k(i-1)*(1/q(i-1) - a(i-1)) )   )...
            +(dt/dz^2)*(D_chem_1/L^2)*( chem_1(i+1) - 2*chem_1(i) + chem_1(i-1))...
            ;
        
    end
end


% Final node

% source 
i = nodesz;
g(i) = (- chem_1(i) + chem_1(i-1)) + (dz*L)*(50*omega/phi)/D_chem_1; % source at the boundary corresponds to S.



end

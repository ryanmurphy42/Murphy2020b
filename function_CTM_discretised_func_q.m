function [g] = function_CTM_discretised_func_q(q,q_old,k,a,nodesz,dt,dz,eta,L_old,L,prolif_law,beta_prolif)

% Density, q

z=0:dz:1;

r = dt/(eta*dz^2*L^2);
dL_dt = (L-L_old)/dt;
g = zeros(nodesz,1);

%% First node

g(1) = k(2)*(1/q(2)-a(2)) - k(1)*(1/q(1)-a(1)); % x0 BC

%% Interior nodes

if prolif_law == 1
    for i = 2:nodesz-1
        C = z(i)*dt/(2*L*dz);
        
        g(i) = - q(i) + q_old(i) ...
            -r*(k(i-1)*(1/q(i-1)-a(i-1)) - 2*k(i)*(1/q(i)-a(i)) + k(i+1)*(1/q(i+1)-a(i+1)))...
            + C*dL_dt*(q(i+1)-q(i-1))...
            + dt*beta_prolif(i)*q_old(i);
    end
    
elseif prolif_law == 2
    
    for i = 2:nodesz-1
        
        C = z(i)*dt/(2*L*dz);
        
        g(i) = - q(i) + q_old(i) ...
            -r*(k(i-1)*(1/q(i-1)-a(i-1)) - 2*k(i)*(1/q(i)-a(i)) + k(i+1)*(1/q(i+1)-a(i+1)))...
            + C*dL_dt*(q(i+1)-q(i-1))...
            + dt*beta_prolif(i)/a(i);
    end
end

%% Final node
i = nodesz;
g(i) = 1/(2*q(i)*L*dz)*(k(i)*(1/q(i)-a(i)) - k(i-1)*(1/q(i-1)-a(i-1))) ...
    + k(i)*(1/q(i)-a(i));


end
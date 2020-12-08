function [L_diag, D_diag, U_diag] = function_CTM_jacobian_chem_1(q,k,nodesz,dt,dz,eta,a,L_old,L,chem_1, chem_1_old,D_chem_1, node_velocity)

% Jacobian for chemical concentration, c

L_diag = zeros(nodesz,1);
D_diag = zeros(nodesz,1);
U_diag = zeros(nodesz,1);

%% First node
D_diag(1) = -1;
U_diag(1) = 1;

%% Internal nodes


for i = 2:nodesz-1
    
    if node_velocity(i) >= 0 %upwinding
        
        L_diag(i) =   -(dt/dz)*node_velocity(i)*(-1)...
            +(dt/dz^2)*(D_chem_1/L^2)*( 1);
        
        D_diag(i) = -1  -(dt/dz)*node_velocity(i) ...
             -dt*(1/(2*eta*(L^2)*dz^2))*(  (1/q(i+1) + 1/q(i)   )*( k(i+1)*(1/q(i+1) - a(i+1)) - k(i)*(1/q(i) - a(i)) )...
            -(1/q(i) + 1/q(i-1)   )*( k(i)*(1/q(i) - a(i)) - k(i-1)*(1/q(i-1) - a(i-1)) )   )...
             +(dt/dz^2)*(D_chem_1/L^2)*( - 2 );        
        
        U_diag(i) =  (dt/dz^2)*(D_chem_1/L^2);
        
        
    elseif node_velocity(i) < 0 %downwinding
        
        L_diag(i) =  (dt/dz^2)*(D_chem_1/L^2)*( 1);
        
        D_diag(i) = -1  -(dt/dz)*node_velocity(i)*(-1) ...
             -dt*(1/(2*eta*(L^2)*dz^2))*(  (1/q(i+1) + 1/q(i)   )*( k(i+1)*(1/q(i+1) - a(i+1)) - k(i)*(1/q(i) - a(i)) )...
            -(1/q(i) + 1/q(i-1)   )*( k(i)*(1/q(i) - a(i)) - k(i-1)*(1/q(i-1) - a(i-1)) )   )...
             +(dt/dz^2)*(D_chem_1/L^2)*( - 2 );        
        
        U_diag(i) =  (dt/dz^2)*(D_chem_1/L^2)  -(dt/dz)*node_velocity(i);
        
    end
    
    
    
end

%% Final node

% constant source
i = nodesz;
L_diag(i) = 1 ;
D_diag(i) =-1 ;

end
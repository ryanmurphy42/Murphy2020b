function [L] = function_CTM_boundary_position(L_old,q,k,a,dt,eta,nodesz,EMT1)

  %% Update the boundary position, L 
 
    
    L = L_old + dt*(1/eta)*(  -2*k(nodesz)*(1/q(nodesz)-a(nodesz)) ) - dt*EMT1/q(nodesz) ;

end
function [L_diag, D_diag, U_diag] = function_CTM_jacobian_q(q,k,nodesz,dt,dz,eta,a,L_old,L)

% Jacobian for density, q

z=0:dz:1;
L_diag = zeros(nodesz,1);
D_diag = zeros(nodesz,1);
U_diag = zeros(nodesz,1);

r = dt/(eta*dz^2*L^2);

dL_dt = (L-L_old)/dt;


%% First node
D_diag(1) = k(1)/q(1)^2;
U_diag(1) = -k(2)/q(2)^2;

%% Internal nodes
for i = 2:nodesz-1
    C = z(i)*dt/(2*L*dz);
    D_diag(i) = -1 - 2*r*k(i)*1/q(i)^2;
    L_diag(i) = r*k(i-1)*1/q(i-1)^2 - C*dL_dt;
    U_diag(i) = r*k(i+1)*1/q(i+1)^2 + C*dL_dt;
end

%% Final node

i = nodesz;
D_diag(i) = 1/(2*L*dz)*(-1/q(i)^2*(k(i)*(1/q(i)-a(i)) - k(i-1)*(1/q(i-1)-a(i-1))) ...
    - 1/q(i)*k(i)/q(i)^2) - k(i)/q(i)^2;
L_diag(i) = 1/(2*q(i)*L*dz)*k(i-1)/q(i-1)^2;



end
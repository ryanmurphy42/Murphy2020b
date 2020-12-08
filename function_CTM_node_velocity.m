function node_velocity = function_CTM_node_velocity(q,k,a,eta,nodesz,dz,L_old,L,dt)

% Update node velocity

node_velocity = zeros(nodesz,1);
z=(0:dz:1);

i = 1;
node_velocity(i) = 1/eta*1/q(i)*1/L^2*(k(i+1)*(1/q(i+1)-a(i+1)) - k(i)*(1/q(i)-a(i)))/(dz)...
        - z(i)/L*(L-L_old)/dt;

for i=2:nodesz-1
    node_velocity(i) = 1/eta*1/q(i)*1/L^2*(k(i+1)*(1/q(i+1)-a(i+1)) - k(i-1)*(1/q(i-1)-a(i-1)))/(2*dz)...
        - z(i)/L*(L-L_old)/dt;
end

i = nodesz;
node_velocity(i) = 1/eta*1/q(i)*1/L^2*(k(i)*(1/q(i)-a(i)) - k(i-1)*(1/q(i-1)-a(i-1)))/(dz)...
        - z(i)/L*(L-L_old)/dt;

end


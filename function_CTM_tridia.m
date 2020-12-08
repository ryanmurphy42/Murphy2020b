function  [ conc ] = function_CTM_tridia( nodesz,a,b,c,d )

% FUNCTION TRIDI
%
% solves the tridiagonal matrix system using the Thomas algorithm. Returns
% the solution to the tridiagonal system.
%
% nodes == number of node points
% a     == upper diagonal
% b     == diagonal
% c     == lower diagonal
% d     == rhs

nodesz = size(b,1);

% initialise solution
conc = zeros(nodesz, 1);

% pass diagonal and rhs
bb = b;
dd = d;

for i = 2:nodesz
    ff = a(i)/bb(i-1);
    bb(i) = bb(i) - c(i-1)*ff;
    dd(i) = dd(i) - dd(i-1)*ff;
end

% perform back substitution
conc(nodesz) = dd(nodesz)/bb(nodesz);
for i = 1:nodesz-1
    j = nodesz-i;
    conc(j) = (dd(j)-c(j)*conc(j+1))/bb(j);
end





end

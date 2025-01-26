function [quad_nodes, quad_weights ] = ...
    Gauss_Jacobi_Radial_Weights_and_Nodes(N_order,p_weight,Radius)
%
% This function outputs the Gauss-Jacobi weights and nodes for a weight of 
% r^p, over the interval [0,R]. N_order is the order of the quadrature
% and is equal to the number of quadrature nodes (and weights).
% To integrate f(r) r^p dr over the interval [0,R], evaluate f(r) at the
% points contained in the vector quad_nodes and then compute the dot
% product of this vector with the quad_weights vector.
% 
% Based on the function spherequad.m available on Mathworks
% written by by Greg von Winckel gregvw(at)math(dot)unm(dot)edu


N = N_order;
k = p_weight;

k1 = k + 1; k2 = k + 2; n = 1:N;  nnk = 2*n+k;
A = [k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
n = 2:N; nnk = nnk(n);
B1 = 4*k1/(k2*k2*(k+3)); nk = n+k; nnk2 = nnk.*nnk;
B = 4*(n.*nk).^2./(nnk2.*nnk2-nnk2);
ab = [A' [(2^k1)/k1; B1; B']]; s = sqrt(ab(2:N,2));
[V,X] = eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I] = sort(diag(X));    
x = (X+1)/2; w = (1/2)^(k1)*ab(1,2)*V(1,I)'.^2;

quad_weights = w * (Radius^(p_weight + 1));
quad_nodes = x * Radius;
end


function [Radial_Basis_Function_Table, r_nodes, r_weights] = ...
    Create_Radial_Basis_Function(S)

Radial_Basis_Function_Table = zeros(S.N_r, (2*S.N_max+1)*S.K_max);
[r_nodes, r_weights] = Gauss_Jacobi_Radial_Weights_and_Nodes(S.N_r,1,S.R);

const = (sqrt(S.N/(pi*S.Tau)))*(1/S.R);

for n = -S.N_max:1:S.N_max
    for k = 1:1:S.K_max
        col_index = (n+S.N_max)*S.K_max + k;
        a = const/Get_Bessel_Evals(n,k,S);
        Radial_Basis_Function_Table(:,col_index) = ...
            a * besselj((n*S.N),(1/S.R)*Get_Bessel_Roots(n,k,S).*r_nodes');
    end
end

% Note that as formulated above, the radial basis functions (with weight r)
% have norm^2 = N/(2*pi*Tau)

end
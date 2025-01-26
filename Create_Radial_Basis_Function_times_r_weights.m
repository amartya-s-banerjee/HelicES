function [Radial_Basis_Function_Table_times_weights] = ...
    Create_Radial_Basis_Function_times_r_weights(S)

Radial_Basis_Function_Table_times_weights = zeros(S.N_r, (2*S.N_max+1)*S.K_max);

const = (sqrt(S.N/(pi*S.Tau)))*(1/S.R);

for n = -S.N_max:1:S.N_max
    for k = 1:1:S.K_max
        col_index = (n+S.N_max)*S.K_max + k;
        a = const/Get_Bessel_Evals(n,k,S);
        Radial_Basis_Function_Table_times_weights(:,col_index) = ...
            a*besselj(n*S.N,(1/S.R)*Get_Bessel_Roots(n,k,S).*S.r_nodes').*S.r_weights';
    end
end

end
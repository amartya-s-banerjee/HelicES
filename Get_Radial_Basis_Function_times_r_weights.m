function [Radial_Basis_Function_times_r_weights] = Get_Radial_Basis_Function_times_r_weights(n,k,S)

col_index = (n+S.N_max)*S.K_max + k;
Radial_Basis_Function_times_r_weights = S.Radial_Basis_Function_Table_times_weights(:,col_index);

end
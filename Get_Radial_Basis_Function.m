function [Radial_Basis_Function] = Get_Radial_Basis_Function(n,k,S)

col_index = (n+S.N_max)*S.K_max + k;
Radial_Basis_Function = S.Radial_Basis_Function_Table(:,col_index);

end
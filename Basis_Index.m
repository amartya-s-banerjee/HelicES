function [i] = Basis_Index(m,n,k,S)

i = (m + S.M_max)*S.K_max*(2*S.N_max + 1) + (n + S.N_max)*S.K_max + k;

end

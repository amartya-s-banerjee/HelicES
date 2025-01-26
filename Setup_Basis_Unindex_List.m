function [Basis_Unindex_List] = Setup_Basis_Unindex_List(S)

Basis_Unindex_List = zeros(S.Num_Basis, 3);

for m = -S.M_max:1:S.M_max
    for n = -S.N_max:1:S.N_max
        for k = 1:1:S.K_max
            i = Basis_Index(m,n,k,S);
            Basis_Unindex_List(i,1) = m;
            Basis_Unindex_List(i,2) = n;
            Basis_Unindex_List(i,3) = k;
        end
    end
end

end
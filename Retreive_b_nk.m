function b_nk = Retreive_b_nk(n,k,matrix)

columns = size(matrix, 2);
N_max_table = (columns - 1)/2;

b_nk = matrix(k, n + N_max_table + 1);

end
function I_Er_grad_nkkprime = Retreive_E_r_grad_innerprod_nkkprime(nprime,S)

I_Er_grad_nkkprime = S.E_r_grad_innerprod_table.E_r_grad_innerprod_table(...
    1:S.K_max, nprime+S.E_r_grad_innerprod_table.N_Er_grad_table_m+1:nprime+...
    S.E_r_grad_innerprod_table.N_Er_grad_table_m+S.K_max);

end
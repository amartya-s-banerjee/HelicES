function I_nkkprime = Retreive_I_nkkprime(nprime,S)

I_nkkprime = S.I_nkkprime_table.I_nkkprime_table(1:S.K_max, ...
    (nprime * S.N + S.I_nkkprime_table.N_rsquared_table_m)*...
    S.I_nkkprime_table.K_rsquared_table_m + 1:(nprime * S.N + ...
    S.I_nkkprime_table.N_rsquared_table_m) * S.I_nkkprime_table.K_rsquared_table_m...
    + S.K_max);

end
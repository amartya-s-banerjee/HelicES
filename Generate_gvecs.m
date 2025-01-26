function gvecs = Generate_gvecs(S)

gvecs = zeros(S.Num_Basis, S.num_k_pts);

for ii = 1:size(S.kptlist,1)
    
    eta = S.kptlist(ii,1);
    nu = S.kptlist(ii,2);
    
    
    const_a = 2*pi*pi*(nu*S.alpha*(2*eta-nu*S.alpha)-eta*eta)/(S.Tau^2);
    const_b = 2*1i*pi*(nu*S.alpha-eta)/(S.Tau^2);
    const_c = 2*1i*pi*S.alpha*(eta-nu*S.alpha)/(S.Tau^2);
    
    ind = 1;
    
    
    if(nu == 0)       
        for m_prime = -S.M_max:S.M_max
            for n_prime = -S.N_max:S.N_max
                for k_prime = 1:S.K_max
                    gvecs(ind, ii) = 0.5*(S.lambda0_mnk(ind)) ...
                        - const_a - const_b*(1i*2*pi*m_prime) ...
                        - const_c * (1i*2*pi*n_prime*S.N);                        
                    ind = ind + 1;
                end
            end
        end
        
    else
        
        % Non-zero values of nu
        for m_prime = -S.M_max:S.M_max
            for n_prime = -S.N_max:S.N_max
                table_portion = Retreive_I_nkkprime(n_prime,S);
                for k_prime = 1:S.K_max
                    gvecs(ind, ii) = 0.5*(S.lambda0_mnk(ind)) ...
                        - const_a - const_b*(1i*2*pi*m_prime) ...
                        - const_c * (1i*2*pi*n_prime*S.N) ...
                        + (nu*nu/(S.R*S.R))*table_portion(k_prime, k_prime) ...
                        + (1i*2*pi*n_prime*S.N)*(1i*nu/(S.R*S.R*pi))* ...
                        table_portion(k_prime, k_prime);
                    ind = ind + 1;
                end
            end
        end
    end
    
end

end
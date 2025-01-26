function Y_hat = Apply_KE_Operator(S,kpt_num,phi_hat)

eta = S.kptlist(kpt_num,1);
nu = S.kptlist(kpt_num,2);

const_a = 2*pi*pi*(nu*S.alpha*(2*eta-nu*S.alpha)-eta*eta)/(S.Tau^2);
const_b = 2*1i*pi*(nu*S.alpha-eta)/(S.Tau^2);
const_c = 2*1i*pi*S.alpha*(eta-nu*S.alpha)/(S.Tau^2);

% Number of bands in the y-direction is the width of phi_hat
num_bands = size(phi_hat, 2);

Y_hat = zeros(S.Num_Basis,num_bands);

for band_iter = 1:num_bands
    
    phi_hat_jth = phi_hat(:,band_iter);
    
    
    phi_hat_jth_times_nprime_vect = S.nprime_vect.*phi_hat_jth;
    
    Y_hat(:,band_iter) = 0.5*(S.lambda0_mnk.*phi_hat_jth) ...
        -const_a*phi_hat_jth - const_b*S.mprime_vect.*phi_hat_jth...
        -const_c*phi_hat_jth_times_nprime_vect;
    
    ind = 1;
    
    % Handle the nu dependent 1/r^2 terms
    if(nu ~= 0)
        
        for m_prime = -S.M_max:S.M_max
            for n_prime = -S.N_max:S.N_max
                table_portion = Retreive_I_nkkprime(n_prime,S);
                
                Y_hat(ind:ind+S.K_max-1,band_iter) = ...
                    Y_hat(ind:ind+S.K_max-1,band_iter) + ...
                    table_portion*((nu*nu*(1/(S.R*S.R))*...
                    phi_hat_jth(ind:ind+S.K_max-1)) + ...
                    ((1i*nu/(S.R*S.R*pi))*...
                    phi_hat_jth_times_nprime_vect(ind:ind+S.K_max-1)));
                
                ind = ind + S.K_max;
            end
        end
    end
    

end


end
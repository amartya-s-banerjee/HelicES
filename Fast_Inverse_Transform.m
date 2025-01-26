function f = Fast_Inverse_Transform(f_hat,S)
% f_hat is a (2M_max+1)*(2N_max+1)*K_max vector whereas 
% f would come out to be a N_r*N_theta_1*N_theta_2 size vector

f = zeros(S.N_theta_1*S.N_theta_2*S.N_r,1);

H_mn_start_index = 1;
f_hat_mn_start_index = 1;

for m = -S.M_max:S.M_max
    for n = -S.N_max:S.N_max
        
        R_kr = S.Radial_Basis_Function_Table(:,...
            (n+S.N_max)*S.K_max + 1 : (n+S.N_max)*S.K_max + S.K_max);
        
        H_mn(H_mn_start_index:H_mn_start_index+S.N_r-1) = ...
            R_kr*f_hat(f_hat_mn_start_index:f_hat_mn_start_index+S.K_max-1);
        
        f_hat_mn_start_index = f_hat_mn_start_index + S.K_max;
        H_mn_start_index = H_mn_start_index + S.N_r;
        
    end
end

for r_loop = 1:S.N_r
    mn_vec = H_mn(r_loop:S.N_r:end);
    f(r_loop:S.N_r:end) = FFT_Theta_1_Theta_2_Inverse_Transform(mn_vec,S);
end

end
function f_hat = Fast_Forward_Transform(f,S)
% f is a (N_theta_1)*(N_theta_2)*N_r sized vector whereas f_hat
% should come out to be (2M_max+1)*(2N_max+1)*K_max sized vector

f_hat = zeros((2*S.M_max+1)*(2*S.N_max+1)*S.K_max, 1);

G_hat_mnr = zeros(S.N_r, (2*S.M_max+1)*(2*S.N_max+1));

for r_loop = 1:S.N_r
    theta_1_theta_2_vec = f(r_loop:S.N_r:end);
    G_hat_mnr(r_loop,:) = ...
        FFT_Theta_1_Theta_2_Forward_Transform(theta_1_theta_2_vec,S);
end

G_mn_index = 1;
f_hat_start_index = 1;

% inc_f = S.K_max*(2*S.N_max+1);
% inc_G = 2*S.N_max + 1;
% 
% for m = -S.M_max:S.M_max
%     
%     temp = S.Radial_Basis_Function_Table_times_weights_transpose*...
%         G_hat_mnr(:, G_mn_start_index : G_mn_start_index + inc_G - 1);
%     temp2 = G_hat_mnr(:, G_mn_start_index : G_mn_start_index + inc_G - 1);
%     whos
%     f_hat(f_hat_start_index : f_hat_start_index + inc_f - 1) = temp(:);
%     
%     G_mn_start_index = G_mn_start_index + inc_G;
%     f_hat_start_index = f_hat_start_index + inc_f;
%     
% end

for m = -S.M_max:S.M_max
    for n = -S.N_max:S.N_max
        
        R_w_kr = S.Radial_Basis_Function_Table_times_weights_transpose(...
            (n+S.N_max)*S.K_max + 1 : (n+S.N_max)*S.K_max + S.K_max, :);
        
        f_hat(f_hat_start_index : f_hat_start_index + S.K_max-1) = R_w_kr*...
            G_hat_mnr(:, G_mn_index);
        
        f_hat_start_index = f_hat_start_index + S.K_max;
        G_mn_index = G_mn_index + 1;
        
    end
end

% Fix the overall constant 
f_hat = f_hat * (2*pi*S.Tau/S.N);

end
function Y_mat = TPA_preconditioner_generalized(S, kpt_num, X_mat)


% Number of bands in the y-direction is the width of phi_hat
num_bands = size(X_mat, 2);
Y_mat = zeros(S.Num_Basis,num_bands);

coeffs = [27.0 18.0 12.0 8.0]; % Original TPA preconditioner
% coeffs = [243.0 162.0 108.0 72.0 48.0 32.0];
% coeffs = [2187.0 1458.0 972.0 648.0 432.0 288.0 192.0 128.0];

order = length(coeffs);
last_coeff = 2.0^order;

for band_iter = 1:num_bands
    
    X_mat_jth = X_mat(:, band_iter);
    KE_op_on_Xmat_jth = Apply_KE_Operator(S, kpt_num, X_mat_jth);
    
    KE_jth = sum(KE_op_on_Xmat_jth.* conj(X_mat_jth));
    g_vec =  S.gvecs(:,kpt_num) / (KE_jth);
    g_vec_pow = g_vec;
    
    g_vec_fac = coeffs(1) * ones(S.Num_Basis,1);
    for ii=2:order
        g_vec_fac = g_vec_fac + coeffs(ii) * g_vec_pow;
        g_vec_pow = g_vec_pow .* g_vec;
    end
    
    
    Y_mat(:, band_iter) = g_vec_fac ./ (g_vec_fac + ...
        last_coeff * g_vec_pow);
    
    Y_mat(:, band_iter) = Y_mat(:, band_iter) .* X_mat_jth;
end

% % Y_mat = X_mat; % This disables the preconditiponer

end
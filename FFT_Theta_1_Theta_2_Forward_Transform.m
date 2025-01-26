function reciprocal_space_fun_1D_fft = FFT_Theta_1_Theta_2_Forward_Transform ...
    (real_space_fun_1D, S)

% Input data is a 1D array -- this is the form we expect to get our data in
% after the radial part is taken care of. In this array each element is
% associated with a theta_1 and a theta_2 value with theta_2 varying faster
% and theta_1 varying slower.


real_space_fun_2D = reshape(real_space_fun_1D, S.N_theta_2, S.N_theta_1);

reciprocal_space_fun_2D = fft2(real_space_fun_2D) / (S.N_theta_2 * S.N_theta_1);

shifted_reciprocal_space_fun = fftshift(reciprocal_space_fun_2D);

mid_big_grid_theta_2 = S.N_theta_2_factor * S.N_max + 1;
mid_big_grid_theta_1 = S.N_theta_1_factor * S.M_max + 1;

smaller_reciprocal_space_fun = shifted_reciprocal_space_fun ...
    ((mid_big_grid_theta_2 - S.N_max):(mid_big_grid_theta_2 + S.N_max), ...
    (mid_big_grid_theta_1 - S.M_max):(mid_big_grid_theta_1 + S.M_max));

reciprocal_space_fun_1D_fft = smaller_reciprocal_space_fun(:);

end
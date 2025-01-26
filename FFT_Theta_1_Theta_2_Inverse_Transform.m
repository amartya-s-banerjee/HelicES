function real_space_fun_1D_fft = FFT_Theta_1_Theta_2_Inverse_Transform ...
    (reciprocal_space_fun_1D, S)

% Input data is a 1D array -- this is the form we expect to get our data in
% after the radial part is taken care of. In this array each element is
% associated with two indices (m,n) with n varying faster (analogously, in
% real space theta_2 varies faster than theta_1. So reshaping this to 2D
% should keep this into consideration.

reciprocal_space_fun_2D = reshape(reciprocal_space_fun_1D, ...
    (2 * S.N_max + 1), (2 * S.M_max + 1));

bigger_reciprocal_space_grid = zeros(S.N_theta_2, S.N_theta_1);

mid_big_grid_theta_2 = S.N_theta_2_factor * S.N_max + 1;
mid_big_grid_theta_1 = S.N_theta_1_factor * S.M_max + 1;

bigger_reciprocal_space_grid(...
    (mid_big_grid_theta_2 - S.N_max):(mid_big_grid_theta_2 + S.N_max), ...
    (mid_big_grid_theta_1 - S.M_max):(mid_big_grid_theta_1 + S.M_max) ) ...
= reciprocal_space_fun_2D;

shifted_bigger_reciprocal_space_grid = ...
    ifftshift(bigger_reciprocal_space_grid);

real_space_fun_2D_fft = ifft2(shifted_bigger_reciprocal_space_grid) * ...
(S.N_theta_1*S.N_theta_2);

real_space_fun_1D_fft = real_space_fun_2D_fft(:);

end
function Y_hat = MATVEC_LOCAL(S,kpt_num, V_LOCAL,phi_hat)

% This MATVEC routine handles the nu dependent 1/r^2 terms spectrally
% by calling Apply_KE_Operatorwhich uses the I(n,k,k') table.
% The local potential is handled pseudo-spectrally of course.


% Apply the KE operator
Y_hat = Apply_KE_Operator(S,kpt_num,phi_hat);

% Number of bands in the y-direction is the width of phi_hat
num_bands = size(phi_hat, 2);


for band_iter = 1:num_bands
    
    phi_hat_jth = phi_hat(:,band_iter);
    
% %     if(S.apply_Ecut_mask == 1)
% %       phi_hat_jth(S.gvecs_Ecut_flags(:,kpt_num)) = 0.0;
% %     end
    
    Inv_transform_phi_hat_jth = Fast_Inverse_Transform(phi_hat_jth,S);
    Y_hat(:,band_iter) = Y_hat(:,band_iter)...
        + Fast_Forward_Transform(Inv_transform_phi_hat_jth.*V_LOCAL,S);
    
    if(S.apply_Ecut_mask == 1)
      Y_hat(S.gvecs_Ecut_flags(:,kpt_num), band_iter) = 0.0;
    end
    
    
end





end
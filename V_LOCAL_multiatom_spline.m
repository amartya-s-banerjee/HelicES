function [all_Cart_coords, V_final, Pseudo_splines_and_deriv] = ...
    V_LOCAL_multiatom_spline(S)

% Hydrogen
% V = -1./(sqrt(X_1D.^2 + Y_1D.^2 + (Z_1D - S.Tau * 0.5).^2));

V_final = zeros(S.Num_Real_Space_Points,S.n_typ);
all_Cart_coords = cell(S.n_typ,2);
Pseudo_splines_and_deriv = cell(S.n_typ,3);

for jj = 1:S.n_typ
    
    m_image_max = ceil(S.Atm(jj).rb/S.Tau) + 1;
    % m_image_max = 0;

    if m_image_max ~= 0

        unit_cell = S.Atm(jj).coords;
        helical_coords = zeros(S.Atm(jj).n_atm_typ*S.N*(2*m_image_max+1),3);
        num_atoms_unit_cell = S.Atm(jj).n_atm_typ;
        index = 1;
        for m = -m_image_max:m_image_max
             for n = 0:S.N-1
               shifted_unit_cell = unit_cell + ...
                   repmat([0 m n/S.N], num_atoms_unit_cell, 1);
               helical_coords(index:index+num_atoms_unit_cell-1,:) = ...
                   shifted_unit_cell;

               index = index + num_atoms_unit_cell;
             end
        end
    else

        helical_coords = S.Atm(jj).coords;      

    end
    
    
    Cart_coords = ConvertHelicalToCartersian(helical_coords,S.alpha,S.Tau);
    all_Cart_coords{jj,1} = S.Atm(jj).typ;
    all_Cart_coords{jj,2} = Cart_coords;
    
    V = zeros(S.Num_Real_Space_Points,1);
    
    pseudo_file_name = strcat('./Pseudopotentials/', ...
        S.Atm(jj).typ, '_empirical_pseudopotential.mat');
    pseudo_struct = load(pseudo_file_name);
    pseudo_arr = pseudo_struct.pseudo;
    
    pseudo_spline = spline(pseudo_arr(:,1), pseudo_arr(:,2));
    pseudo_spline_deriv = fnder(pseudo_spline);
    Pseudo_splines_and_deriv{jj,1} = S.Atm(jj).typ;
    Pseudo_splines_and_deriv{jj,2} = pseudo_spline;
    Pseudo_splines_and_deriv{jj,3} = pseudo_spline_deriv;

    % Add a check for the cutoff for consistency of pseudopotential
    % Not really that important, but can prevent users from using 
    % wrong pseudopotentials and/or parameters
    assert((pseudo_arr(end,1) == S.Atm(jj).rb), ...
        'Local pseudopotential cutoff specified by user in input file does not match pseudopotential file data !');
             
    for ii = 1:size(Cart_coords,1)
        r_vecs = sqrt((S.X_1D-Cart_coords(ii,1)).^2 + ...
            (S.Y_1D-Cart_coords(ii,2)).^2 + ...
            (S.Z_1D-Cart_coords(ii,3)).^2);
        cutoff_mask = r_vecs < S.Atm(jj).rb;
        V = V + fnval(pseudo_spline, r_vecs).*cutoff_mask;
    end
    
    V_final(:,jj) = V;
    
end

V_final = sum(V_final,2);

end
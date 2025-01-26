clc; clear all;
format longg;

t_total = tic;

lastn = maxNumCompThreads(1);

num_worker_heuristic = 4;


%fname = './Input_Files/GNR_armchair_passivated';
fname = './Input_Files/Input_C6,6_no_cylic_alpha_0_test';
%fname = './Input_Files/Input_diagnostic';
%fname = './Input_Files/Input_test_forces';
%fname = './Input_Files/Moire_Retry_1';
%fname = './Input_Files/Input_C16,16_no_cylic_alpha_0';

diary_fname = sprintf('%s_diary.txt',fname);
if (exist(diary_fname, 'file'))
    delete(diary_fname);
end
% [~,~] = system(sprintf('rm -f %s_diary.txt',fname));

diary(diary_fname);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%              HelicES                %%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%% (c) Ab Initio Simulations Lab, UCLA %%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%  Supported by DOE Grant #DE-SC0023432 %%%%%%%%%%%%%%')	
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

fprintf('\n Input file is : %s\n', fname);
dt = datestr(now,'mmmm dd, yyyy HH:MM:SS AM');
fprintf(' Running on %s .\n',dt);
v_ = version; c_ = computer;
fprintf(' Environment is MATLAB %s on %s .\n',v_, c_);

if (~ispc)
    [~, hostname] = system('hostname');
else
    hostname = getenv('COMPUTERNAME');
end

if (isdeployed)
    fprintf(' Running in deployed mode (mcc compiled) on %s', hostname);
else
    fprintf(' Running in MATLAB session on %s', hostname);
end

fprintf('\n');

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

fprintf('\n Input file dump: \n');
type(fname); % Dump the entire input file
%fprintf('\n');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% Run the input parser and initialize all inputs
t1 = tic;
fprintf('\n Parsing input file... ');
S = Input_Parser(fname);
fprintf('\n Input File parsing completed (%f s).', toc(t1));


% Get the basis unindexing list and add it to S
t1 = tic;
fprintf('\n Generating basis un-index list ... ');
S.Basis_Unindex_List = Setup_Basis_Unindex_List(S);
fprintf('Done (%f s).', toc(t1));

% Get the bessel roots and values table and add it to S
t1 = tic;
fprintf('\n Generating Bessel Roots and Evals ... ');
[S.Bessel_Roots, S.Bessel_Evals] = Setup_Bessel_Roots_and_Evals(S);
fprintf('Done (%f s).', toc(t1));

% Get the radial basis function table and add it to S,
% along with the nodes and weights
t1 = tic;
fprintf('\n Generating radial basis funcs. and quad. weights/nodes ... ');
[S.Radial_Basis_Function_Table, S.r_nodes, S.r_weights] = ...
    Create_Radial_Basis_Function(S);
fprintf('Done (%f s).', toc(t1));

% Get the radial basis function times radial weights table
t1 = tic;
fprintf('\n Generating radial basis func. table modifications ... ');
S.Radial_Basis_Function_Table_times_weights = ...
    Create_Radial_Basis_Function_times_r_weights(S);

% Get the transpose of the above matrix
S.Radial_Basis_Function_Table_times_weights_transpose = ...
    Create_Radial_Basis_Function_times_r_weights_transpose(S);
fprintf('Done (%f s).', toc(t1));


% Setup real space mesh
t1 = tic;
fprintf('\n Generating real space mesh  and weights ... ');
[S.theta_1_nodes, S.theta_2_nodes, S.mesh_3D_r_1D, ...
    S.mesh_3D_theta_1_1D, S.mesh_3D_theta_2_1D, ...
    S.X_1D, S.Y_1D, S.Z_1D, S.weights_3D] = Setup_Real_Space_Mesh(S);
fprintf('Done (%f s).', toc(t1));

% Set up the k points using a Monkhorst-Pack grid
t1 = tic;
fprintf('\n Generating reciprocal space mesh ... ');
[S.num_k_pts, S.weights_k_pts, S.kptlist] = ...
    Create_List_and_Weights_of_k_points(S);
fprintf('Done (%f s).', toc(t1));

% Adding MATVEC prereqs lambda0_mnk, mprime vect and nprime vect to S
t1 = tic;
fprintf('\n Generating MATVEC prereqs ... ');
[S.lambda0_mnk, S.mprime_vect, S.nprime_vect] = MATVEC_prereq(S);
% Needed in pseudospactral handling of 1/r^2 terms
S.one_by_rsquared = 1./(S.mesh_3D_r_1D.*S.mesh_3D_r_1D);
fprintf('Done (%f s).', toc(t1));



% Compute the local potential
t1 = tic;
fprintf('\n Generating local potential ... ');
[S.all_Cart_cords, S.V, S.Pseudo_splines_and_deriv] = ...
    V_LOCAL_multiatom_spline(S);

fprintf('Done (%f s).', toc(t1));

% Set up the 'g' vectors for TPA preconditioner, energy cutoff, etc.
t1 =tic;
fprintf('\n Computing g vectors and mask ... ');
S.gvecs = Generate_gvecs(S);
% Compute mask
S.gvecs_Ecut_flags = (S.gvecs > S.E_cut);
fprintf('Done (%f s).', toc(t1));

% Only engage the Ecut mask under specific well tested cases
% i.e., EIGS with no cyclic symmetry (or diagnostic mode).
% Note that with eta points some eigenstates near 0 might not come out
% correct, so we need to issue a warning.
S.apply_Ecut_mask = 0;
if((S.flag_E_cut == 1) && (S.EigSolver <= 0 ) && (S.N == 1))
    S.apply_Ecut_mask = 1;
    fprintf('\n Energy cutoff (%f Ha) based mask to be imposed on wavefunctions.', ...
        S.E_cut);
    if(S.num_k_pts > 1)
        fprintf('\n\n Note: Number of eta points > 1 => Ecut might cause issues.');
        fprintf('\n Eigenstates close to zero may not be always calculated correctly.');
    end
    
    % Print a list of effective number of basis functions
    fprintf('\n\n Effective number of basis functions at each k-point:')
    for ii=1:S.num_k_pts
        num_eff_basis_ii = S.Num_Basis - sum(S.gvecs_Ecut_flags(:,ii));
        fprintf('\n %d. (eta = %f, nu = %f) --> %d.',ii, ...
            S.kptlist(ii,1), S.kptlist(ii,2), num_eff_basis_ii);
    end
    fprintf('\n');
    
else
    fprintf('\n Energy cutoff based mask will not be imposed on wavefunctions.\n');
    if(S.flag_E_cut == 1)
        if( S.N > 1 )
            fprintf('\n Note: Mask disabled as problem has cyclic symmetry.');
        end
        if(S.EigSolver ~= 0)
            fprintf('\n Note: Mask disabled as Eigensolver is not EIGS.');
        end
    end
end


% % Modification to default values of number of electrons
% % S.Nelectron = S.n_atm;

fprintf('\n \n');

% % Output some simulation details
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
if(S.EigSolver < 0)
    fprintf('\n Running in diagnostic mode...')
end
fprintf('\n Some simulation details:')
fprintf('\n M_max = %d, N_max = %d, K_max = %d.',S.M_max,S.N_max,S.K_max);
fprintf('\n Total basis set size (without mask) =  %d.',S.Num_Basis);
fprintf('\n N_r = %d, N_theta_1 = %d, N_theta_2 =  %d.',...
    S.N_r, S.N_theta_1, S.N_theta_2);
fprintf('\n Number of real space points =  %d.',S.Num_Real_Space_Points);
fprintf('\n Reciprocal space mesh has %d k-points.',S.num_k_pts);
fprintf('\n Number of atoms in fundamental domain =  %d.',S.n_atm);
fprintf('\n Number of electrons in fundamental domain =  %d.',S.Nelectron);
fprintf('\n Memory required for eigenfunction storage = %f MB.',...
    S.Num_Basis * S.NEV * S.num_k_pts * 16/(1024*1024));

if(S.EigSolver < 0)
    
    % Diagnostic mode
    % Evaluate forward and inverse transform time and errors
    N_run = 10;
    t_inv = 0.0;  t_fwd = 0.0; transform_err = 0.0;
    for run_i = 1:N_run
        f_hat = rand(S.Num_Basis,1) + (1i)*rand(S.Num_Basis,1);
        tic_transform_inv = tic;
        f = Fast_Inverse_Transform(f_hat,S);
        t_inv = t_inv + toc(tic_transform_inv);
        
        tic_transform_fwd = tic;
        f_hat_hat = Fast_Forward_Transform(f,S);
        t_fwd = t_fwd + toc(tic_transform_fwd);
        
        transform_err = transform_err +norm(f_hat_hat - f_hat);
        
    end
    fprintf('\n ------------------------------------- \n');
    fprintf('\n Full basis results:');
    fprintf('\n Average forward transform time = %f s.', t_fwd / N_run );
    fprintf('\n Average inverse transform time = %f s.', t_inv / N_run );
    fprintf('\n Average transform error (L2 norm) = %.4g \n', ...
        transform_err / N_run );
    
    
    if(S.apply_Ecut_mask == 1)
        fprintf('\n\n Masked basis results:');
        for ii =  1:size(S.kptlist,1)
            eta = S.kptlist(ii,1);
            nu = S.kptlist(ii,2);
            fprintf('\n %d. eta = %f, nu = %f:',ii, eta,nu);
            N_run = 10;
            t_inv = 0.0;  t_fwd = 0.0; transform_err = 0.0;
            
            for run_i = 1:N_run
                f_hat = rand(S.Num_Basis,1) + (1i)*rand(S.Num_Basis,1);
                f_hat(S.gvecs_Ecut_flags(:,ii)) = 0.0; % Apply mask
                
                tic_transform_inv = tic;
                f = Fast_Inverse_Transform(f_hat,S);
                t_inv = t_inv + toc(tic_transform_inv);
                
                tic_transform_fwd = tic;
                f_hat_hat = Fast_Forward_Transform(f,S);
                t_fwd = t_fwd + toc(tic_transform_fwd);
                transform_err = transform_err +norm(f_hat_hat - f_hat);
            end
            
            fprintf('\n Average forward transform time = %f s.', t_fwd / N_run );
            fprintf('\n Average inverse transform time = %f s.', t_inv / N_run );
            fprintf('\n Average transform error (L2 norm) = %.4g \n', ...
                transform_err / N_run );
            
        end
        
    end
    
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\n Diagnostic run completed. \n');
    
else
    
    
    if(S.EigSolver == 0)
        eigsolver = 'EIGS';
    elseif(S.EigSolver == 1)
        eigsolver = 'LOBPCG with TPA Preconditioner';
    else
        eigsolver = '#########';
    end
    fprintf('\n Eigensolver is %s.',eigsolver);
    fprintf('\n Eigensolver will compute %d eigenstates at each k-point.',S.NEV);
    fprintf('\n Eigensolver convergence tolerance is %.3g .',S.EigSolverTol);
    
    fprintf('\n \n');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\n \n');
    
    
    % Clean up older pools
    delete(gcp('nocreate'));
    
    % Launch new pool
    fprintf(' ');
    poolobj = parpool(num_worker_heuristic,'IdleTimeout', Inf ) ;
    
    
    % Diagonalize
    fprintf('\n -----------------------------------');
    % fprintf('\n Note: Memory required for eigenfunction storage = %f MB.',...
    %    S.Num_Basis * S.NEV * S.num_k_pts * 16 /(1024*1024));
    fprintf('\n Performing diagonalization at every eta,nu point:\n');
    energy_table = zeros(S.NEV, S.num_k_pts);
    all_eigen_vectors = zeros(S.Num_Basis, S.NEV, S.num_k_pts);
    
    
    t1 = tic;
    parfor (ii = 1:size(S.kptlist,1), num_worker_heuristic)
        % %for ii = 1:size(S.kptlist,1)
        eta = S.kptlist(ii,1);
        nu = S.kptlist(ii,2);
        
        % Definition of matrix-vector product
        Hfun = @(x_vec_mat) MATVEC_LOCAL(S, ii, S.V,x_vec_mat);
        
        if(S.EigSolver == 0)
            % % % EIGS
            v0 = rand(S.Num_Basis,1) + (1i)*rand(S.Num_Basis,1);
            
            % This seems to add to computation time but appears necessary
            % to compute eigenstates near zero correctly
            if(S.apply_Ecut_mask == 1)
                v0(S.gvecs_Ecut_flags(:,ii)) = 0.0;
            end
            
            [blockVectorX, eig_vals, conv_flag] = ...
                eigs(Hfun, S.Num_Basis, S.NEV, 'smallestreal', ...
                'IsFunctionSymmetric', 1, ...
                'Tolerance',S.EigSolverTol, 'MaxIterations', 5000,...
                'StartVector', v0, 'FailureTreatment', 'keep');
            
            % Save eigenvalues and eigenvectors
            energy_table(:,ii) = diag(eig_vals);
            all_eigen_vectors(:,:, ii) = blockVectorX;
            
            
            if(conv_flag == 0)
                fprintf('\n Kpt Iteration %d of %d (eta = %f, nu = %f) --> EIGS converged.\n',...
                    ii, size(S.kptlist,1), eta,nu);
            else
                fprintf('\n Kpt Iteration %d of %d (eta = %f, nu = %f) --> EIGS did not converge.\n',...
                    ii, size(S.kptlist,1), eta,nu);
            end
            
        else
            % % % LOBPCG
            precfun =  @(x_vec_mat) TPA_preconditioner_generalized(S, ii, x_vec_mat);
            blockVectorX = rand(S.Num_Basis,S.NEV) + 1i*rand(S.Num_Basis,S.NEV);
            [blockVectorX, eig_vals, failureFlag, lambdaHistory,...
                residualNormsHistory] =  lobpcg(blockVectorX, Hfun, [], ...
                precfun, S.EigSolverTol, 5000, 0);
            sz_residual = size(residualNormsHistory);
            fprintf('\n Kpt Iteration %d of %d (eta = %f, nu = %f) --> %d LOBPCG iters.\n',...
                ii, size(S.kptlist,1), eta,nu, sz_residual(2));
            
            % Save eigenvalues and eigenvectors
            energy_table(:,ii) = (eig_vals(1:S.NEV));
            all_eigen_vectors(:,:, ii) = blockVectorX;
            
        end
        
    end
    fprintf('\n -----------------------------------');
    fprintf('\n\n Total diagonalization time = %f s. \n',toc(t1));
    
    % energy_table
    mys=sort(reshape(energy_table,S.NEV * size(S.kptlist,1),1));
    
    
    
    % Solve for Fermi energy at T_electron
    t1 = tic;
    fprintf('\n Calculating Fermi level and occupations ... ');
    T_electron_ref = 315.77;
    beta = 27.211384523 / ((8.617343e-5)*T_electron_ref);
    
    FermiEnergyEvaluator = ...
        @(lambda_fg)TotalNumberOfElectrons(lambda_fg, energy_table, ...
        S.Nelectron, beta, S.num_k_pts, S.weights_k_pts);
    
    lambda_f = fzero(FermiEnergyEvaluator,0);
    
    % Calculate occupations at T_electron
    occ = 1./(1+exp(beta*(energy_table-lambda_f)));
    S.occ = occ;
    
    fprintf('Done. (%f s)', toc(t1));
    
    
    % Calculate the electron density
    fprintf('\n Calculating electron density ... ');
    
% %     % Serial routine
% %     t1 = tic;
% %     S.rho = zeros(S.Num_Real_Space_Points, 1);
% %     for kpt=1:size(S.kptlist,1)
% %         for band_iter = 1:S.NEV
% %             abs_Inv_transform_wave_fun = ...
% %                 abs(Fast_Inverse_Transform(all_eigen_vectors(:,band_iter,kpt) ,S));
% %             
% %             S.rho = S.rho + ...
% %                 2.0 * S.weights_k_pts(kpt) * S.occ(band_iter, kpt) ...
% %                 * (abs_Inv_transform_wave_fun .* abs_Inv_transform_wave_fun);
% %         end
% %         
% %     end
% %     
% %     fprintf('Done. (%f s)\n', toc(t1));
    
    % Parallel routine
    t1 = tic;
    transformed_wave_fun_sum = zeros(S.Num_Real_Space_Points, size(S.kptlist,1));
    
    parfor (kpt = 1:size(S.kptlist,1), num_worker_heuristic)
        temp_eigvecs = all_eigen_vectors(:,:,kpt);
        for band_iter = 1:S.NEV
            abs_Inv_transform_wave_fun = ...
                abs(Fast_Inverse_Transform(temp_eigvecs(:,band_iter) ,S));
            
            transformed_wave_fun_sum(:,kpt) = ...
                transformed_wave_fun_sum(:,kpt) + ...
                2.0 * S.weights_k_pts(kpt) * S.occ(band_iter, kpt) ...
                * (abs_Inv_transform_wave_fun .* abs_Inv_transform_wave_fun);
        end        
    end
    
    S.rho = sum(transformed_wave_fun_sum,2);
    transformed_wave_fun_sum = [];
    
    
    fprintf('Done. (%f s)\n', toc(t1));
    num_electron = sum(S.rho .* S.weights_3D);
    fprintf('\n Integral of rho (= no. of electrons in fundamental domain) = %1.14f .\n',...
        num_electron);
    
    
    %plot_vector_real_space(S.rho,0.1,S);
    
    % Calculate the forces (only local component here)
    fprintf('\n Calculating atomic forces  ... ');
    t1 = tic;
    
    S.FD_atoms_Cart_Coords = cell(S.n_typ,2);
    S.FD_atoms_forces = cell(S.n_typ, 2);
    
    for atom_type = 1:S.n_typ
        
        S.FD_atoms_Cart_Coords{atom_type,1} = S.Atm(atom_type).typ;
        S.FD_atoms_forces{atom_type,1} = S.Atm(atom_type).typ;
        S.FD_atoms_Cart_Coords{atom_type,2} = ...
            ConvertHelicalToCartersian(S.Atm(atom_type).coords,S.alpha,S.Tau);
        FD_atoms_type_forces = zeros(S.Atm(atom_type).n_atm_typ,3);
        
        for atom_iter = 1:S.Atm(atom_type).n_atm_typ
            
            FD_atom_type_Cart_Coords = S.FD_atoms_Cart_Coords{atom_type,2};
            diff_x = S.X_1D - FD_atom_type_Cart_Coords(atom_iter,1);
            diff_y = S.Y_1D - FD_atom_type_Cart_Coords(atom_iter,2);
            diff_z = S.Z_1D - FD_atom_type_Cart_Coords(atom_iter,3);
            
            r_vecs = sqrt(diff_x .* diff_x + diff_y .* diff_y + ...
                diff_z .* diff_z);
            vloc_part = fnval(S.Pseudo_splines_and_deriv{atom_type,3}, ...
                r_vecs)./r_vecs;
            rho_vloc = S.rho .* vloc_part;
            
            integrand_force_x = rho_vloc .* diff_x;
            integrand_force_y = rho_vloc .* diff_y;
            integrand_force_z = rho_vloc .* diff_z;
            
            FD_atoms_type_forces(atom_iter, 1) = ...
                sum(integrand_force_x .* S.weights_3D);
            FD_atoms_type_forces(atom_iter, 2) = ...
                sum(integrand_force_y .* S.weights_3D);
            FD_atoms_type_forces(atom_iter, 3) = ...
                sum(integrand_force_z .* S.weights_3D);
            
        end
        
        S.FD_atoms_forces{atom_type,2} = FD_atoms_type_forces;
        
    end
    
    fprintf('Done. (%f s)', toc(t1));
    
    S.FD_atoms_forces{:,2}
    
    
    % Calculate HOMO, LUMO, and BandGap
    % Do a rough calculation of the HOMO, LUMO and bandgap
    fprintf('\n \n');
    kptlist = S.kptlist;
    weights_k_pts = S.weights_k_pts;
    
    sz_eigvals = size(energy_table);
    energy_tablea = reshape(energy_table, sz_eigvals(1) * sz_eigvals(2), 1);
    energy_tablea_with_fermi = [energy_tablea; lambda_f];
    energy_tablea_with_fermi_sorted = sort(energy_tablea_with_fermi);
    ind_fermi = find(energy_tablea_with_fermi_sorted == lambda_f);
    
    eigval_HOMO = energy_tablea_with_fermi_sorted(ind_fermi - 1);
    eigval_LUMO = energy_tablea_with_fermi_sorted(ind_fermi + 1);
    
    ind_HOMO = find(energy_tablea == eigval_HOMO);
    ind_LUMO = find(energy_tablea == eigval_LUMO);
    
    sz_occ = size(occ);
    occ_a = reshape(occ, sz_occ(1) * sz_occ(2), 1);
    
    occ_HOMO = occ_a(ind_HOMO);
    occ_LUMO = occ_a(ind_LUMO);
    
    [row_HOMO, col_HOMO] = find(energy_table == eigval_HOMO);
    [row_LUMO, col_LUMO] = find(energy_table == eigval_LUMO);
    
    % Compute the band energy
    Eband = 0;
    for kpt=1:size(S.kptlist,1)
        Eband = Eband + ...
            weights_k_pts(kpt) * 2.0 * sum(energy_table(:,kpt).*occ(:,kpt));
    end
    
    
    fprintf('\n Number of electrons = %1.3f, Number of states = %1.3f .', ...
        S.Nelectron, S.NEV);
    fprintf('\n At electronic temperature = %1.3f :', T_electron_ref);
    fprintf('\n Electronic band energy = %1.14f Ha.\n',Eband);
    
    fprintf('\n HOMO eigenvalue = %1.9f Ha, occ = %1.9f .', ...
        eigval_HOMO, occ_HOMO);
    fprintf('\n HOMO eigenvalue occurs at eta = %1.9f, nu = %1.9f,',...
        kptlist(col_HOMO,1), kptlist(col_HOMO,2));
    fprintf('\n                           eigenvalue no. = %d .\n', row_HOMO);
    fprintf('\n Fermi level = %1.9f Ha\n', ...
        energy_tablea_with_fermi_sorted(ind_fermi));
    fprintf('\n LUMO eigenvalue = %1.9f Ha, occ = %1.9f .', ...
        eigval_LUMO, occ_LUMO);
    fprintf('\n LUMO eigenvalue occurs at eta = %1.9f, nu = %1.9f,',...
        kptlist(col_LUMO,1), kptlist(col_LUMO,2));
    fprintf('\n                           eigenvalue no. = %d .\n', row_LUMO);
    
    band_gap = eigval_LUMO - eigval_HOMO;
    
    fprintf('\n Band-gap = %1.9f Ha = %1.9f eV.\n\n', band_gap, (band_gap * 27.2114));
    
    fprintf('\n Saving eigenvalues, occupations, etc. ...');
    eigval_fname = strcat(fname,'_alpha_',num2str(S.alpha),...
        '_Eigval_Data.mat');
    
    Num_electron = S.Nelectron;
    
    save(eigval_fname,'energy_table', 'kptlist','weights_k_pts','occ','lambda_f',...
        'eigval_LUMO', 'eigval_HOMO', 'band_gap', 'occ_HOMO', 'occ_LUMO',...
        'Num_electron', 'T_electron_ref');
    fprintf(' Done.\n\n');
    
    fprintf(' ');
    delete(poolobj);
    
end

fprintf('\n \n Total elapsed time = %f s. \n\n',toc(t_total));
diary off;

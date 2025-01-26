function [S] = Input_Parser(filename)

fileID = fopen(filename, 'r');

textscan(fileID,'%s',1,'delimiter','\n');
R_tau = textscan(fileID,'%f %f',1,'delimiter','\n');
alpha_N = textscan(fileID,'%f %f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
table_n_table_k = textscan(fileID,'%f %f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
flag_E_cut = textscan(fileID,'%f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
E_cut = textscan(fileID,'%f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
m_n_k = textscan(fileID,'%f %f %f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
N_r_theta1_theta2_factor = textscan(fileID,'%f %f %f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
EigSolver = textscan(fileID,'%f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
NEV_EigSolverTol = textscan(fileID,'%f %f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
input_n_k_pts = textscan(fileID,'%f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
flag_time_reversal_sym = textscan(fileID,'%f',1,'delimiter','\n');

textscan(fileID,'%s',1,'delimiter','\n');
textscan(fileID,'%s',1,'delimiter','\n');
atom_types = textscan(fileID,'%f',1,'delimiter','\n');

n_typ = atom_types{1,1};
n_atm=0; Nelectron = 0;
%Atoms = [];

Atm = repmat(struct('typ','C','Z','4',...
    'coords',[1 0.5 0.5],'n_atm_typ',1,'rb', 15), n_typ,1) ; % Preallocate
% Atom types info stored in structure Atm with fields typ, Z, coords
for ii = 1:n_typ
    
    textscan(fileID,'%s',1,'delimiter','\n');
    textscan(fileID,'%s',1,'delimiter','\n');
    
    symbol = textscan(fileID,'%s',1,'delimiter','\n');
    chg = textscan(fileID,'%f',1,'delimiter','\n');
    natm = textscan(fileID,'%f',1,'delimiter','\n');
    
    textscan(fileID,'%s',1,'delimiter','\n');
    textscan(fileID,'%s',1,'delimiter','\n');
    
    rb = textscan(fileID,'%f',1,'delimiter','\n');
    rb = rb{1,1};
    
    textscan(fileID,'%s',1,'delimiter','\n');
    textscan(fileID,'%s',1,'delimiter','\n');
    
    C_atom = textscan(fileID,'%f %f %f',natm{1,1},'delimiter','\n');
    
    if (feof(fileID)~=1 && ii==n_typ) || (feof(fileID)==1 && ii~=n_typ)
        error('Error: Inconsistency between no. of types of atoms and the ones actually specified!')
    end
    
    typ = symbol{1,1};
    Z = chg{1,1};
    n_atm_typ = natm{1,1};
    n_atm = n_atm + natm{1,1};
    Nelectron = Nelectron + Z*n_atm_typ;
    coords = [C_atom{:,1},C_atom{:,2},C_atom{:,3}];
    % Atoms = [Atoms;coords]; % full list of atom coords of all types
    Atm(ii)=struct('typ',typ,'Z',Z,'coords',coords,'n_atm_typ',n_atm_typ,'rb',rb);
    
end

fclose(fileID);

R = R_tau{1,1};
Tau = R_tau{1,2};
alpha = alpha_N{1,1};
N = 1/alpha_N{1,2};

if mod(N,1) ~= 0
    N = round(N,0);
end

% Modify the atomic helical coordinates to be within fundamental domain
mod_flag = 0;
for ii=1:n_typ
    coords = Atm(ii).coords;
    
    for jj=1:size(coords,1)
        
        while(coords(jj,2) > 1.0)
            coords(jj,2) =  coords(jj,2) - 1.0;
            mod_flag = 1;
        end
        
        while(coords(jj,2) < 0.0)
            coords(jj,2) =  coords(jj,2) + 1.0;
            mod_flag = 1;
        end
        
        while(coords(jj,3) > 1.0/N)
            coords(jj,3) =  coords(jj,3) - 1.0/N;
            mod_flag = 1;
        end
        
        while(coords(jj,3) < 0.0)
            coords(jj,3) =  coords(jj,3) + 1.0/N;
            mod_flag = 1;
        end
        
    end
    
    Atm(ii).coords = coords;
end

if(mod_flag == 1)
    fprintf('\n \n Some atomic coordinates were mapped back into the fundamental domain.');
    fprintf('\n Modified helical coordinates are:\n');
    
    for ii=1:n_typ
        fprintf('\n Type: %s, Z_v: %f',Atm(ii).typ, Atm(ii).Z);
        coords = Atm(ii).coords;
        
        for jj=1:size(coords,1)
            fprintf('\n %1.14f\t%1.14f\t%1.14f',coords(jj,1), coords(jj,2), ...
                coords(jj,3));
        end
    end
    
    fprintf('\n\n');
end



% Highest values used to generate b_nk_table, I_nkkprime_table, or
% E_r_grad_table
% % save_table_M_max = table_m_table_n_table_k{1,1};
save_table_N_max = table_n_table_k{1,1};
save_table_K_max = table_n_table_k{1,2};



% Load I_nkkprime_table
if(N ~= 1)
    fname_I_nkkprime_table = sprintf('I_nkkprime_%i_%i_table.mat',...
        save_table_N_max, save_table_K_max);
    I_nkkprime_table = load(fname_I_nkkprime_table);
else
    I_nkkprime_table = []; % Not used in this case !
end


% Load E_r_grad_innerprod_table
% fname_E_r_grad_innerprod = ...
%     sprintf('E_r_grad_innerprod_%i_%i_table.mat', save_table_N_max,...
%     save_table_K_max);
% E_r_grad_innerprod_table = load(fname_E_r_grad_innerprod);
E_r_grad_innerprod_table = [];

flag_E_cut = flag_E_cut{1,1};
E_cut = E_cut{1,1};



M_max = m_n_k{1,1};
N_max = m_n_k{1,2};
K_max = m_n_k{1,3};

N_max_mod = N_max * N; % This is done to take care of cyclic symmetries

if(flag_E_cut == 1)
    
    % Load b_nk_table for Ecut based limits
    fname_b_nk_table = sprintf('b_nk_%i_%i_table.mat', save_table_N_max, ...
        save_table_K_max);
    b_nk_table_loaded = load(fname_b_nk_table);
    % If Ecut based, choose/modify M_max, N_max, K_max here.
    save_table_M_max = save_table_N_max;
    for m = -save_table_M_max:save_table_M_max
        for n = -save_table_N_max:save_table_N_max
            for k = 1:save_table_K_max
                lamda_0_mnk = ...
                    abs(Retreive_b_nk(n,k, b_nk_table_loaded.b_nk_table)/R)^2 ...
                    + abs((2*pi/Tau)*(m-alpha*n))^2; % We don't use the factor of N here !
                if lamda_0_mnk < E_cut
                    if (abs(m) > M_max)
                        M_max = abs(m);
                        % % fprintf('\n \n Found m : %d %d %d \n',m, n, k);
                    end
                    
                    if(abs(n) > N_max_mod)
                        N_max_mod = abs(n);
                        % % fprintf('\n \n Found n : %d %d %d \n',m, n, k);
                    end
                    if (abs(k) > K_max)
                        K_max = abs(k);
                        % % fprintf('\n \n Found k : %d %d %d \n',m, n, k);
                    end
                end
            end
        end
    end
    
end

N_max = fix(N_max_mod/N); % Account for cyclic symmetry here

if((M_max ~= m_n_k{1,1}) || (N_max ~= m_n_k{1,2}) || (K_max ~= m_n_k{1,3}))
    fprintf('\n Ecut (= %f Ha) imposition changes basis set maxes.', E_cut)
    fprintf('\n Revised M_max = %d, N_Max = %d, K_max = %d.', ...
        M_max, N_max, K_max);
end


Num_Basis = (2*M_max+1)*(2*N_max+1)*K_max;




N_theta_1_factor = N_r_theta1_theta2_factor{1,1};
N_theta_2_factor = N_r_theta1_theta2_factor{1,2};
N_r_factor = N_r_theta1_theta2_factor{1,3};

N_theta_1 = N_theta_1_factor*2*M_max+1;
N_theta_2 = N_theta_2_factor*2*N_max+1;
N_r = N_r_factor*K_max;
Num_Real_Space_Points = N_r * N_theta_1 * N_theta_2;

NEV = NEV_EigSolverTol{1,1};
EigSolverTol = NEV_EigSolverTol{1,2};

ip_num_k_pts = input_n_k_pts{1,1};
flag_time_reversal_sym = flag_time_reversal_sym{1,1};



% Creating structure of all the variables
S = struct('R',R,'Tau',Tau,'alpha',alpha,'N',N, ...
    'M_max',M_max,'N_max', N_max,'K_max',K_max,'Num_Basis',Num_Basis,...
    'N_r',N_r,'N_theta_1',N_theta_1,...
    'N_theta_2',N_theta_2,'Num_Real_Space_Points', Num_Real_Space_Points, ...
    'N_r_factor',N_r_factor,'N_theta_1_factor',N_theta_1_factor, ...
    'N_theta_2_factor',N_theta_2_factor, ...
    'EigSolver', EigSolver, 'NEV',NEV, 'EigSolverTol', EigSolverTol, ...
    'Input_num_k_pts',ip_num_k_pts,...
    'Flag_Time_Reversal_Sym',flag_time_reversal_sym, ...
    'flag_E_cut',flag_E_cut,...
    'E_cut',E_cut, ...
    'save_table_N_max', save_table_N_max, ...
    'save_table_K_max',save_table_K_max,...
    'I_nkkprime_table',I_nkkprime_table,'E_r_grad_innerprod_table',...
    E_r_grad_innerprod_table,...
    'n_atm',n_atm,...
    'Nelectron',Nelectron,...
    'n_typ', n_typ,'Atm',Atm,'input_fname',filename);

end
function plot_vector_real_space(vec_real_space,spacing_h, S)

plot_vec = vec_real_space;

% Create the X,Y,Z nodes
theta_1_nodes = (0:S.N_theta_1 - 1)' / S.N_theta_1;
% theta_1_weights = ones(S.N_theta_1, 1) / S.N_theta_1;
theta_2_nodes = (1/S.N) * ((0:S.N_theta_2 - 1)' / S.N_theta_2);
% theta_2_weights = (1/S.N) * (ones(S.N_theta_2, 1) / S.N_theta_2);

% Set up the nodes (in 3 dimensions)
mesh_3D_r_1D = repmat(S.r_nodes,S.N_theta_1*S.N_theta_2,1);

temp_mesh_3D_theta_2_1D = repelem(theta_2_nodes,S.N_r);
mesh_3D_theta_2_1D = repmat(temp_mesh_3D_theta_2_1D,S.N_theta_1,1);

mesh_3D_theta_1_1D = repelem(theta_1_nodes,S.N_r*S.N_theta_2);

X_1D = mesh_3D_r_1D.*cos(2*pi*(S.alpha*mesh_3D_theta_1_1D + ...
    mesh_3D_theta_2_1D));
Y_1D = mesh_3D_r_1D.*sin(2*pi*(S.alpha*mesh_3D_theta_1_1D + ...
    mesh_3D_theta_2_1D));
Z_1D = mesh_3D_theta_1_1D * S.Tau;

% Need to interpolate this 3D Data
interp = scatteredInterpolant(X_1D,Y_1D,Z_1D, plot_vec);
[X_3D,Y_3D,Z_3D] = meshgrid(-S.R:spacing_h:S.R,...
    -S.R:spacing_h:S.R,0:spacing_h:S.Tau);
vq = interp(X_3D,Y_3D,Z_3D);

isosurface(X_3D,Y_3D,Z_3D, vq);
daspect([1 1 1]);
end
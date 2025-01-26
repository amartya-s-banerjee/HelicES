function [theta_1_nodes, theta_2_nodes, mesh_3D_r_1D, ...
    mesh_3D_theta_1_1D, mesh_3D_theta_2_1D, X_1D, Y_1D, Z_1D, ...
    weights_3D] = Setup_Real_Space_Mesh(S)

% Setup real space mesh 
 
theta_1_nodes = (0:S.N_theta_1 - 1)' / S.N_theta_1;
% theta_1_weights = ones(S.N_theta_1, 1) / S.N_theta_1;
theta_2_nodes = (1/S.N) * ((0:S.N_theta_2 - 1)' / S.N_theta_2);
% theta_2_weights = (1/S.N) * (ones(S.N_theta_2, 1) / S.N_theta_2);

% Set up the nodes (in 3 dimensions)

% Radial nodes at each point in 3D mesh
mesh_3D_r_1D = repmat(S.r_nodes,S.N_theta_1*S.N_theta_2,1);

% Theta2 nodes at each point in 3D mesh
temp_mesh_3D_theta_2_1D = repelem(theta_2_nodes,S.N_r);
mesh_3D_theta_2_1D = repmat(temp_mesh_3D_theta_2_1D,S.N_theta_1,1);

% Theta1 nodes at each point in 3D mesh
mesh_3D_theta_1_1D = repelem(theta_1_nodes,S.N_r*S.N_theta_2);

% X,Y,Z nodes at each point in 3D mesh
X_1D = mesh_3D_r_1D.*cos(2*pi*(S.alpha*mesh_3D_theta_1_1D + ...
    mesh_3D_theta_2_1D));
Y_1D = mesh_3D_r_1D.*sin(2*pi*(S.alpha*mesh_3D_theta_1_1D + ...
    mesh_3D_theta_2_1D));
Z_1D = mesh_3D_theta_1_1D * S.Tau;


% Set up the quadrature weights in 3D : Only the radial weights vary,
% the theta_1 and theta_2 weights are constant
% weights_3D below integrates to volume of FD = pi*R*R/(N*tau)
r_weights_3D = repmat(S.r_weights,S.N_theta_1*S.N_theta_2,1);
weights_3D = r_weights_3D * (2 * pi * S.Tau) / ...
    (S.N_theta_1 * S.N_theta_2 * S.N);


end

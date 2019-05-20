function [R_dqp_12, R_dqp_3, R_xyz, S_0lab_dir] = plugin_APS_34IDC(theta_bl, chi_bl, phi_bl, delta_bl, gamma_bl, rocking_increment, rocking_angle)
% converts 34-ID-C angles to right-handed convention and creates matrices
fprintf('plugin_APS_34IDC');
% convert reflection SPEC angles to a right-handed coordinate system for new reflection
fprintf('\n...converting angles...');
theta = theta_bl;
chi = 90 - chi_bl;
phi = phi_bl;
delta = delta_bl;
gamma = -gamma_bl;

% rotation matrices
fprintf('\n...creating rotation matrices...');
R_dqp_12 = roty(delta)*rotx(gamma);
R_xyz = roty(theta)*rotz(chi)*rotx(phi); % rotation matrix to rotate a vector in sample coordiantes into lab coordinates
if strcmp(rocking_angle, 'dtheta')
    R_dqp_3 = roty(-rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
elseif strcmp(rocking_angle, 'dphi')
    R_dqp_3 = rotx(-rocking_increment); % it's the negative in rocking increment because we scan from negative to positive
end

% S_0lab direction for APS
fprintf('\n...setting S_0lab direction...');
S_0lab_dir = [0; 0; 1];

fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function LS_shape_DCS = DCS_to_LS(S)
fprintf('DCS_to_LS\n...transforming DCS shape to LS shape...');
%% Declaring variables from shape structure
% detector conjugated space shape
DCS_shape_DRS = S.DCS_shape_DRS;

% detector parameters
D = S.D;
d = S.d;

% sample pixel size
p_sam = S.p_sam;

% wavelength
lambda = S.lambda;

% orthogonal grid coordinates
N1grid = S.N1grid;
N2grid = S.N2grid;
N3grid = S.N3grid;

% rotation matrices
R_dqp_12 = S.R_dqp_12;
R_dqp_3 = S.R_dqp_3;

% scattering vector
Q_lab = S.Q_lab;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making basis vectors
% finding size of the array and creating lab space dimensions
[N1, N2, N3] = size(DCS_shape_DRS);

% making x, y, z lab space vectors
x = p_sam*[1; 0; 0];
y = p_sam*[0; 1; 0];
z = p_sam*[0; 0; 1];

% making q_1', q_2', q_3' detector reciprocal space vectors
q_1p = 2*pi/lambda*d/D*R_dqp_12*[1; 0; 0];
q_2p = 2*pi/lambda*d/D*R_dqp_12*[0; 1; 0];
q_3p = R_dqp_3*Q_lab-Q_lab;

% make x', y', and z' detector conjugated space basis vectors, adapted from Berenguer et al. PRB 88, 144101 (2013).
V_DRS = dot(cross(q_3p, q_2p), q_1p)*N1*N2*N3;
xp = 2*pi*cross(N2*q_2p, N3*q_3p)./V_DRS;
yp = 2*pi*cross(N3*q_3p, N1*q_1p)./V_DRS;
zp = 2*pi*cross(N1*q_1p, N2*q_2p)./V_DRS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map from detector conjugated space to lab space:
% make the T_DCS_to_LS matrix
T_DCS_to_LS = [xp, yp, zp]\[x, y, z]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map non-orthogonal coordinates to orthogonal coordinates 
N1gridp = T_DCS_to_LS(1,1)*N1grid + T_DCS_to_LS(1,2)*N2grid + T_DCS_to_LS(1,3)*N3grid;
N2gridp = T_DCS_to_LS(2,1)*N1grid + T_DCS_to_LS(2,2)*N2grid + T_DCS_to_LS(2,3)*N3grid;
N3gridp = T_DCS_to_LS(3,1)*N1grid + T_DCS_to_LS(3,2)*N2grid + T_DCS_to_LS(3,3)*N3grid;

% interpolate data in the detector conjugated frame to the lab frame
LS_shape_DCS = interp3(N1grid, N2grid, N3grid, DCS_shape_DRS, N1gridp, N2gridp, N3gridp, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
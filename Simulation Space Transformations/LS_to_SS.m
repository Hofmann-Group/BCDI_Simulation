function SS_shape_LS = LS_to_SS(S, LS_shape)
% maps a LS shape to SS
fprintf('\nLS_to_SS\n...transforming LS shape to SS shape...');
%% Making basis vectors
% making x, y, z lab space vectors
x = S.p_sam*[1; 0; 0];
y = S.p_sam*[0; 1; 0];
z = S.p_sam*[0; 0; 1];

% making x_sam , y_sam , z_sam sample space vectors
x_sam = S.R_xyz*x;
y_sam = S.R_xyz*y;
z_sam = S.R_xyz*z;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map from sample space to lab space
% make the T_LS_to_SS matrix
T_LS_to_SS = [x, y, z]\[x_sam, y_sam, z_sam]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% mapping lab space coordinates to sample space coordinates
samXgrid = T_LS_to_SS(1,1)*S.N1grid+ T_LS_to_SS(1,2)*S.N2grid+ T_LS_to_SS(1,3)*S.N3grid;
samYgrid = T_LS_to_SS(2,1)*S.N1grid+ T_LS_to_SS(2,2)*S.N2grid+ T_LS_to_SS(2,3)*S.N3grid;
samZgrid = T_LS_to_SS(3,1)*S.N1grid+ T_LS_to_SS(3,2)*S.N2grid+ T_LS_to_SS(3,3)*S.N3grid;

% interpolate data in the lab frame to the sample frame
SS_shape_LS = interp3(S.N1grid, S.N2grid, S.N3grid, LS_shape, samXgrid, samYgrid, samZgrid, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
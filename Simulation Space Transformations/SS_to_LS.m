function LS_shape_SS = SS_to_LS(S)
fprintf('SS_to_LS\n...transforming SS shape to LS shape...');
%% Declaring variables from shape structure
% sample space shape
shape = S.SS_shape;

% rotation matrx
R_xyz= S.R_xyz;

% sample pixel size
p_sam = S.p_sam;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making basis vectors
% making x, y, z lab space vectors
x = p_sam*[1; 0; 0];
y = p_sam*[0; 1; 0];
z = p_sam*[0; 0; 1];

% making x_sam , y_sam , z_sam sample space vectors
x_sam = R_xyz*x;
y_sam = R_xyz*y;
z_sam = R_xyz*z;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map from sample space to lab space
% make the T_SS_to_LS matrix
T_SS_to_LS = [x_sam, y_sam, z_sam]\[x, y, z]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% mapping sample space coordinates to lab space coordinates
labXgrid = T_SS_to_LS(1,1)*S.realXgrid + T_SS_to_LS(1,2)*S.realYgrid + T_SS_to_LS(1,3)*S.realZgrid;
labYgrid = T_SS_to_LS(2,1)*S.realXgrid + T_SS_to_LS(2,2)*S.realYgrid + T_SS_to_LS(2,3)*S.realZgrid;
labZgrid = T_SS_to_LS(3,1)*S.realXgrid + T_SS_to_LS(3,2)*S.realYgrid + T_SS_to_LS(3,3)*S.realZgrid;

% interpolate data in the sample frame to the lab frame
LS_shape_SS = interp3(S.realXgrid, S.realYgrid, S.realZgrid, shape, labXgrid, labYgrid, labZgrid, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
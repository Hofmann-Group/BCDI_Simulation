function LS_shape_SS = SS_to_LS(S)
% maps a SS shape to LS
fprintf('SS_to_LS\n...transforming SS shape to LS shape...');
%% Map from sample space to lab space
% make the T_SS_to_LS matrix
T_SS_to_LS = [S.x_sam, S.y_sam, S.z_sam]\[S.x, S.y, S.z]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% mapping sample space coordinates to lab space coordinates
N1gridlab = T_SS_to_LS(1,1)*S.N1grid + T_SS_to_LS(1,2)*S.N2grid + T_SS_to_LS(1,3)*S.N3grid;
N2gridlab = T_SS_to_LS(2,1)*S.N1grid + T_SS_to_LS(2,2)*S.N2grid + T_SS_to_LS(2,3)*S.N3grid;
N3gridlab = T_SS_to_LS(3,1)*S.N1grid + T_SS_to_LS(3,2)*S.N2grid + T_SS_to_LS(3,3)*S.N3grid;

% interpolate data in the sample frame to the lab frame
LS_shape_SS = interp3(S.N1grid, S.N2grid, S.N3grid, S.SS_shape, N1gridlab, N2gridlab, N3gridlab, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
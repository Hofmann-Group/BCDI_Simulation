function DCS_shape_LS = LS_to_DCS(S)
% maps a LS shape to DCS
fprintf('LS_to_DCS\n...transforming LS shape to DCS shape...');
%% Map from lab space to detector conjugated space:
% make the T_LS_to_DCS matrix
T_LS_to_DCS = [S.x, S.y, S.z]\[S.xp, S.yp, S.zp]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map orthogonal coordinates to non-orthogonal coordinates
N1gridp = T_LS_to_DCS(1,1)*S.N1grid + T_LS_to_DCS(1,2)*S.N2grid + T_LS_to_DCS(1,3)*S.N3grid;
N2gridp = T_LS_to_DCS(2,1)*S.N1grid + T_LS_to_DCS(2,2)*S.N2grid + T_LS_to_DCS(2,3)*S.N3grid;
N3gridp = T_LS_to_DCS(3,1)*S.N1grid + T_LS_to_DCS(3,2)*S.N2grid + T_LS_to_DCS(3,3)*S.N3grid;

% interpolate data in lab frame to detector conjugated frame
DCS_shape_LS = interp3(S.N1grid, S.N2grid, S.N3grid, S.LS_shape_SS, N1gridp, N2gridp, N3gridp, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function LS_shape_DCS = DCS_to_LS(S)
% maps a DCS shape to LS
fprintf('DCS_to_LS\n...transforming DCS shape to LS shape...');
%% Map from detector conjugated space to lab space:
% make the T_DCS_to_LS matrix
T_DCS_to_LS = [S.xp, S.yp, S.zp]\[S.x, S.y, S.z]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map non-orthogonal coordinates to orthogonal coordinates 
N1gridp = T_DCS_to_LS(1,1)*S.N1grid + T_DCS_to_LS(1,2)*S.N2grid + T_DCS_to_LS(1,3)*S.N3grid;
N2gridp = T_DCS_to_LS(2,1)*S.N1grid + T_DCS_to_LS(2,2)*S.N2grid + T_DCS_to_LS(2,3)*S.N3grid;
N3gridp = T_DCS_to_LS(3,1)*S.N1grid + T_DCS_to_LS(3,2)*S.N2grid + T_DCS_to_LS(3,3)*S.N3grid;

% interpolate data in the detector conjugated frame to the lab frame
LS_shape_DCS = interp3(S.N1grid, S.N2grid, S.N3grid, S.DCS_shape_DRS, N1gridp, N2gridp, N3gridp, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function DRS_shape_RLS = RLS_to_DRS(S)
fprintf('RLS_to_DRS\n...transforming RLS shape to DRS shape...');
%% Declaring variables from shape structure
% reciprocal lab space shape
I_RLS_shape_LS = S.I_RLS_shape_LS;

% orthogonal grid coordinates
N1grid = S.N1grid;
N2grid = S.N2grid;
N3grid = S.N3grid;

% non-orthogonal grid coordinates
N1gridp = S.N1gridp;
N2gridp = S.N2gridp;
N3gridp = S.N3gridp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate data in the reciprocal lab frame to the detector reciprocal frame
DRS_shape_RLS = interp3(N1grid, N2grid, N3grid, I_RLS_shape_LS, N1gridp, N2gridp, N3gridp, 'linear', 0);
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function DRS_shape_RLS = RLS_to_DRS(S)
% maps a RLS shape to DRS
fprintf('RLS_to_DRS\n...transforming RLS shape to DRS shape...');
%% Interpolate data in the reciprocal lab frame to the detector reciprocal frame
DRS_shape_RLS = interp3(S.N1grid, S.N2grid, S.N3grid, S.I_RLS_shape_LS, S.N1gridp, S.N2gridp, S.N3gridp, 'linear', 0);
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


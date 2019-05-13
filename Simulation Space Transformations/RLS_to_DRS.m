function DRS_shape_RLS = RLS_to_DRS(S)
fprintf('RLS_to_DRS\n...transforming RLS shape to DRS shape...');
%% Declaring variables from shape structure
% reciprocal lab space shape
I_RLS_shape_LS = S.I_RLS_shape_LS;

% reciprocal space grid
reciprocalXgrid = S.reciprocalXgrid;
reciprocalYgrid = S.reciprocalYgrid;
reciprocalZgrid = S.reciprocalZgrid;

% detector reciprocal space grid
det_reciprocalXgrid = S.det_reciprocalXgrid;
det_reciprocalYgrid = S.det_reciprocalYgrid;
det_reciprocalZgrid = S.det_reciprocalZgrid;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate data in the reciprocal lab frame to the detector reciprocal frame
DRS_shape_RLS = interp3(reciprocalXgrid, reciprocalYgrid, reciprocalZgrid, I_RLS_shape_LS, det_reciprocalXgrid, det_reciprocalYgrid, det_reciprocalZgrid, 'linear', 0);
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


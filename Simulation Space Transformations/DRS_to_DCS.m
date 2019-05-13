function DCS_shape_DRS = DRS_to_DCS(S)
fprintf('DRS_to_DCS\n...transforming DRS shape to DCS shape...');
%% Declaring variables from shape structure
% reciprocal lab space shape
RLS_shape_LS = S.RLS_shape_LS;

% reciprocal space grid
reciprocalXgrid = S.reciprocalXgrid;
reciprocalYgrid = S.reciprocalYgrid;
reciprocalZgrid = S.reciprocalZgrid;

% detector reciprocal space grid
det_reciprocalXgrid = S.det_reciprocalXgrid;
det_reciprocalYgrid = S.det_reciprocalYgrid;
det_reciprocalZgrid = S.det_reciprocalZgrid;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating detector conjugated space shape using inverse FT of RLS shape
DCS_shape_DRS= interp3(reciprocalXgrid, reciprocalYgrid, reciprocalZgrid, RLS_shape_LS, det_reciprocalXgrid, det_reciprocalYgrid, det_reciprocalZgrid, 'linear', 0); % using FT of shape rather than intensity
DCS_shape_DRS = fftshift(ifftn(ifftshift(DCS_shape_DRS)));
DCS_shape_DRS = abs(DCS_shape_DRS);

% take the twin
F = ifftshift(fftn(fftshift(DCS_shape_DRS)));
DCS_shape_DRS = fftshift(ifftn(ifftshift(conj(F)))); 
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
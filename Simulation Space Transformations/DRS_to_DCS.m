function DCS_shape_DRS = DRS_to_DCS(S)
fprintf('DRS_to_DCS\n...transforming DRS shape to DCS shape...');
%% Declaring variables from shape structure
% reciprocal lab space shape
RLS_shape_LS = S.RLS_shape_LS;

% orthogonal grid coordinates
N1grid = S.N1grid;
N2grid = S.N2grid;
N3grid = S.N3grid;

% non-orthogonal grid coordinates
N1gridp = S.N1gridp;
N2gridp = S.N2gridp;
N3gridp = S.N3gridp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating detector conjugated space shape using inverse FT of RLS shape
DCS_shape_DRS= interp3(N1grid, N2grid, N3grid, RLS_shape_LS, N1gridp, N2gridp, N3gridp, 'linear', 0); % using FT of shape rather than intensity
DCS_shape_DRS = fftshift(ifftn(ifftshift(DCS_shape_DRS)));
DCS_shape_DRS = abs(DCS_shape_DRS);

% take the twin
F = ifftshift(fftn(fftshift(DCS_shape_DRS)));
DCS_shape_DRS = fftshift(ifftn(ifftshift(conj(F)))); 
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
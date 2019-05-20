function DCS_shape_DRS = DRS_to_DCS(S)
% maps a DRS shape to DCS
fprintf('DRS_to_DCS\n...transforming DRS shape to DCS shape...');
%% Calculating detector conjugated space shape using inverse FT of RLS shape
DCS_shape_DRS= interp3(S.N1grid, S.N2grid, S.N3grid, S.RLS_shape_LS, S.N1gridp, S.N2gridp, S.N3gridp, 'linear', 0); % using FT of shape rather than intensity
DCS_shape_DRS = fftshift(ifftn(ifftshift(DCS_shape_DRS)));
DCS_shape_DRS = abs(DCS_shape_DRS);

% take the twin
F = ifftshift(fftn(fftshift(DCS_shape_DRS)));
DCS_shape_DRS = fftshift(ifftn(ifftshift(conj(F)))); 
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
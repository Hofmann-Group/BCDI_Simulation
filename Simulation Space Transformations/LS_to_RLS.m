function [RLS_shape_LS,I_RLS_shape_LS] = LS_to_RLS(S)
% maps a LS shape to RLS
fprintf('LS_to_RLS\n...transforming LS shape to RLS shape...');
%% Calculating reciprocal lab space and intensity
% taking fourier transform of lab shape
RLS_shape_LS = fftshift(fftn(ifftshift(S.LS_shape_SS)));

% calculating intensity of diffraction pattern
I_RLS_shape_LS = RLS_shape_LS.*conj(RLS_shape_LS);

fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function [RLS_shape_LS,I_RLS_shape_LS] = LS_to_RLS(S)
fprintf('LS_to_RLS\n...transforming LS shape to RLS shape...');
%% Declaring variables from shape structure
% lab space shape
shape = S.LS_shape_SS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating reciprocal lab space and intensity
% taking fourier transform of lab shape
RLS_shape_LS = fftshift(fftn(ifftshift(shape)));

% calculating intensity of diffraction pattern
I_RLS_shape_LS = RLS_shape_LS.*conj(RLS_shape_LS);

fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
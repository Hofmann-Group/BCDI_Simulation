function save_arrays_and_tiff(S, dir, binary_threshold)
fprintf('\n\n<<<SAVING FILES>>>\n\n');
fprintf('save_arrays_and_tiff\n...making new folder...');
%% Choosing dtheta or dphi in file name
if strcmp(S.rocking_angle, 'dtheta')
    rocking_angle = round(S.dtheta, 5);
elseif strcmp(S.rocking_angle, 'dphi')
    rocking_angle = round(S.dphi, 5);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making a new folder to save the arrays
% choosing the folder AND file name, feel free to change
folder_name = [S.name,'_(', num2str(S.hkl(1)), num2str(S.hkl(2)), num2str(S.hkl(3)),')_', num2str(-S.gamma_spec),'_gamma_',num2str(S.delta_spec),'_delta_',num2str(rocking_angle),'_',S.rocking_angle];

% making the array folder
mkdir(dir, folder_name);

% creating array folder directory
folder_dir = [dir, '/', folder_name];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving -SAM, -LAB, -AMP, -PH, and -SUP arrays in detector conjugated space
fprintf(['\n...saving -SAM, -LAB, -AMP, -PH, and -SUP (binarized DCS shape) arrays in "',folder_dir, '"...']);
% sample space shape
array = S.SS_shape;
save([folder_dir, '/', folder_name, '-SAM.mat'],'array');

% lab space shape
array = S.LS_shape_SS;
save([folder_dir, '/', folder_name, '-LAB.mat'],'array');

% amplitude (not intensity)
array = abs(S.DCS_shape_DRS);
save([folder_dir, '/', folder_name, '-AMP.mat'],'array');

% phase
array = angle(S.DCS_shape_DRS);
save([folder_dir, '/', folder_name, '-PH.mat'],'array');

% binary shape
array = abs(S.DCS_shape_DRS) > binary_threshold;
save([folder_dir, '/', folder_name, '-SUP.mat'],'array');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving intensity TIFF file
temp = [dir, '/', folder_name];
fprintf(['\n...saving intensity TIFF file in "', dir, '"...']);
if exist(sprintf('%s',temp,'.tif'), 'file')==2
    delete(sprintf('%s',temp,'.tif'));
end
% saveastiff(single(S.DRS_shape_RLS), sprintf('%s',temp,'.tif'));
saveastiff(single(flip(S.DRS_shape_RLS,1)), sprintf('%s',temp,'.tif')); 
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
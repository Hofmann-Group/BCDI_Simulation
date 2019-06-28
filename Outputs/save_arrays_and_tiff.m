function save_arrays_and_tiff(S, dir, binary_threshold)
% saves -SAM, -LAB, -AMP, -PH, -BIN and TIFF files
fprintf('\n\n<<<SAVING FILES>>>\n\n');
fprintf('save_arrays_and_tiff\n...making new folder...');
%% Choosing dtheta or dphi in file name
if strcmp(S.rocking_axis, 'dtheta')
    rocking_angle = round(S.dtheta, 5);
elseif strcmp(S.rocking_axis, 'dphi')
    rocking_angle = round(S.dphi, 5);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making a new folder to save the arrays
% choosing the folder AND file name, feel free to change
folder_name = [S.name,'_(', num2str(S.hkl(1)), num2str(S.hkl(2)), num2str(S.hkl(3)),')_', num2str(-S.gamma_bl),'_gamma_',num2str(S.delta_bl),'_delta_',num2str(rocking_angle),'_',S.rocking_axis];

% making the array folder
mkdir(dir, folder_name);

% creating array folder directory
folder_dir = [dir, '/', folder_name];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving -SAM, -LAB, -AMP, -PH, and -SUP arrays in detector conjugated space
fprintf(['\n...saving -SAM, -LAB, -AMP, -PH, and -BIN arrays in "',folder_dir, '"...']);
% sample space shape
array = S.SS_shape;
save([folder_dir, '/', folder_name, '-SAM.mat'],'array');

% lab space shape
array = S.LS_shape_SS;
save([folder_dir, '/', folder_name, '-LAB.mat'],'array');

% amplitude of detector conjugated shape
array = abs(S.DCS_shape_DRS);
save([folder_dir, '/', folder_name, '-AMP.mat'],'array');

% phase of detector conjugated shape
array = angle(S.DCS_shape_DRS);
save([folder_dir, '/', folder_name, '-PH.mat'],'array');

% binarized detector conjugated shape
array = abs(S.DCS_shape_DRS) > binary_threshold;
save([folder_dir, '/', folder_name, '-BIN.mat'],'array');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving intensity TIFF file
temp = [dir, '/', folder_name];
fprintf(['\n...saving intensity TIFF file in "', dir, '"...']);
if exist(sprintf('%s',temp,'.tif'), 'file')==2
    delete(sprintf('%s',temp,'.tif'));
end

saveastiff(single(flip(S.DRS_shape_RLS, 1)), sprintf('%s', temp, '.tif')); % need to flip along the first dimension because the detector "sees" the incoming beam as if it were looking down on the sample
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
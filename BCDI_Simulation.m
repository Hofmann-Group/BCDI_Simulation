%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCDI Simulation 
% Version 1.0
% May 2019
% Written By: David Yang
% University of Oxford, Dept. of Engineering Science
% 
% PURPOSE: To simulate a BCDI experiment for a given shape and reflection
% 
% USER-DEFINED SECTIONS:
% 1. Sample Space Shape Details
% - collects size and shape details about the sample space (SS) object
% 
% 2. BCDI Measurement Variables
% - user specifies the beamline measurement details and reflection
% 
% 3. Saving file options
% - save SAM, LAB, AMP, PH, SUP (binarized DCS shape) arrays and TIFF file
% 
% 4. Plot Options
% - user can set plot parameters
% 
% 5. BCDI Space Plots
% - plots the different BCDI spaces characteristic of a BCDI measurement
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
fprintf('                       BCDI Simulation\n');
fprintf('                          Version 1.0\n');
fprintf('                            May 2019\n');
fprintf('                     Written By: David Yang\n');
fprintf('      University of Oxford, Dept. of Engineering Science\n');
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Sample Space Shape Details
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 2^8; % square detector length in number of pixels, should ideally be a power of 2
length = 2^6; % number of pixels that the object occupies, should ideally be <= N/4
a_0 = 3.16522e-10; % for tungsten in m
size = 1000*10^-9; % total length of object in sample space along the x, y or z axis in m
p_sam = size/length; % pixel size of object in sample space in m
name = 'Cylinder'; % name of the shape (it will be on all files generated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x_pix, y_pix, z_pix] = meshgrid(-(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2, -(N-1)/2:(N-1)/2); % creating a meshgrid for SS coordinates
SS_shape = zeros(N, N, N); % creating empty shape array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECT SAMPLE SHAPE %%%%%%%%%%%%%%%%%%%%%%%%%%%
SS_shape((x_pix.^2 + z_pix.^2) < (length/2)^2 & abs(y_pix) < length/2) = complex(1,0); % cylinder
% SS_shape(abs(x_pix) < length & abs(y_pix) < length & abs(z_pix) < length) = complex(1,0); % cube
% SS_shape((abs(x_pix) + abs(y_pix) + abs(z_pix)) < length) = complex(1,0); % octahedral
% SS_shape(((x_pix/2).^2 + (y_pix/1).^2 + (z_pix/3).^2) < (length/2)^2) = complex(1,0); % elippsoid
% SS_shape((x_pix.^2 + z_pix.^2) < (length/2)^2 & abs(y_pix) < length/2 & x_pix < 0) = complex(1,0); % semi-cylinder
% SS_shape(y_pix <= 0 & y_pix >= -length/2 & x_pix <= 0 & x_pix >= -length & z_pix <= 0 & z_pix >= -length & x_pix.^2 + z_pix.^2 <= (length/2)^2 |...
%     y_pix <= 0 & y_pix >= -length/2 & x_pix <= 0 & x_pix >= -length & z_pix >= 0 & z_pix <= length & abs(x_pix) < (length/2) & abs(z_pix) < (length/2) | ...
%     y_pix <= 0 & y_pix >= -length/2 & x_pix >= 0 & x_pix <= length & z_pix <= 0 & z_pix >= -length & x_pix.^2 + z_pix.^2 <= (length/2)^2 |...
%     y_pix <= 0 & y_pix >= -length/2 & x_pix >= 0 & x_pix <= length & z_pix >= 0 & z_pix <= length & abs(x_pix) < (length/2) & abs(z_pix) < (length/2) |...
%     y_pix > 0 & y_pix <= length/2 & x_pix >= 0 & x_pix <= length & z_pix <= 0 & z_pix >= -length & abs(x_pix) < (length/2) & abs(z_pix) < (length/2) |...
%     y_pix > 0 & y_pix <= length/2 & x_pix >= 0 & x_pix <= length & z_pix >= 0 & z_pix <= length & x_pix.^2 + z_pix.^2 <= (length/2)^2 |...
%     (y_pix > 0 & (abs(x_pix.^2) + abs(z_pix).^2 + abs(y_pix).^2 <= (length/2)^2))) ...
%     = complex(1,0); % test shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. BCDI Measurement Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%% CRYSTAL ORIENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
hkl = [-1; 2; 0]; % pick a reflection
U_0 = eye(3); % orientation matrix, otherwise set to eye(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 34-ID-C PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 55e-6; % detector pixel size in m
lambda = (12.398/10.0)/10*10^-9; % wavelength in m
D = p_sam*N*d/lambda; % detector distance in m, put absolute value if known
rocking_angle = 'dtheta'; % choose 'dtheta' or 'dphi' rocking 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3. Saving file options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_files = 0; % 1 to save -SAM, -LAB, -AMP, -PH, -SUP and TIFF files
save_files_dir = 'Examples'; % the folder to put the TIFF file and array folder
binary_threshold = 0.5; % binary threshold applied to abs(DCS_shape_DRS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 4. Plot Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT VIEWPOINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
viewpoint = [-180, -90]; % viewpoint = [az, el], x-y plane
% viewpoint = [-180, 0]; % viewpoint = [az, el], x-z plane
% viewpoint = [-90, 0]; % viewpoint = [az, el], y-z plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION PLOT OPTIONS %%%%%%%%%%%%%%%%%%%%%%%
dir = 'Examples'; % directory to reconstructions
folder_name = 'Rec-Cylinder_(-120)_dtheta-00274-ERlrHIOlr2000-NM-SW'; % reconstruction folder (which should be the root of the individual file names)
twin = 1; % 0 or 1 to flip LAB shape, toggle to flip reconstruction so that it matches simulation
p_sam_REC = p_sam; % a Command Window Output from reconstruction, in m (not nm). Comment/leave null if unknown/not needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THRESHOLDING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_threshold = 0.3; % for creating binary masks from thresholded amplitude to calcuate centre of mass
plot_threshold = 0.3; % amplitude threshold for figure plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 5. BCDI Space Plots
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SAMPLE SPACE (SS) %%%%%%%%%%%%%%%%%%%%%%%%%%
plot_SS_figures = 1; % 1 to show sample space figure
plot_SS_axes = 1; % 1  to show axes in the figure
plot_SS_beams = 1; % 1  to show S_0lab, S_lab and Q_lab in the figure
plot_SS_grid = 1; % 1  to gridlines in the figure
SS_shape_plot = 1; % 1 to show SS shape in the figure
SS_shape_LS_plot = 1; % 1 to show SS shape calculated from LS in the figure
SS_shape_REC_plot = 1; % 1 to show SS shape from reconstruction in the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT LAB SPACE (LS) %%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_LS_figures = 1; % 1 to show lab space figure
plot_LS_axes = 1; % 1  to show axes in the figure
plot_LS_beams = 1; % 1  to show S_0lab, S_lab and Q_lab in the figure
plot_LS_grid = 1; % 1  to gridlines in the figure
LS_shape_SS_plot = 1; % 1 to show LS shape calculated from SS in the figure
LS_shape_DCS_plot = 1; % 1 to show LS shape calculated from DCS in the figure
LS_shape_REC_plot = 1; % 1 to show LS shape from reconstruction in the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% PLOT RECIPROCAL LAB SPACE (RLS) %%%%%%%%%%%%%%%%%%%%%
plot_RLS_figures = 1; % 1 to show reciprocal lab space figure
plot_RLS_axes = 1; % 1  to show axes in the figure
plot_RLS_beams = 1; % 1  to show S_0lab, S_lab and Q_lab in the figure
plot_RLS_grid = 1; % 1  to gridlines in the figure
plot_RLS_log10scale = 1; % 1 for log10 intensity, 0 for absolute intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% PLOT DETECTOR RECIPROCAL SPACE (DRS) %%%%%%%%%%%%%%%%%%
plot_DRS_figures = 1; % 1 to show detector reciprocal space figure
plot_DRS_axes = 1; % 1  to show axes in the figure
plot_DRS_beams = 1; % 1  to show S_0lab, S_lab and Q_lab in the figure
plot_DRS_grid = 1; % 1  to gridlines in the figure
plot_DRS_log10scale = 1; % 1 for log10 intensity, 0 for absolute intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% PLOT DETECTOR CONJUGATED SPACE (DCS) %%%%%%%%%%%%%%%%%%
plot_DCS_figures = 1; % 1 to show detector conjugated space figure
plot_DCS_axes = 1; % 1  to show axes in the figure
plot_DCS_beams = 1; % 1  to show S_0lab, S_lab and Q_lab in the figure
plot_DCS_grid = 1; % 1  to gridlines in the figure
DCS_shape_DRS_plot = 1; % 1 to show DCS shape calculated from DRS in the figure
DCS_shape_LS_plot = 1; % 1 to show DCS shape calculated from LS in the figure
DCS_shape_REC_plot = 1; % 1 to show LS shape from reconstruction in the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating paths to scripts in folders
addpath(genpath('APS angles'));
addpath(genpath('Coordinates'));
addpath(genpath('Examples'));
addpath(genpath('Plots'));
addpath(genpath('Simulation Space Transformations'));
addpath(genpath('Tiff'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Placing variables in structure array
fprintf('\n\n<<<SETTING UP SIMULATION EXPERIMENT>>>\n\n');

% making the structure that contains all the experiment information
S = struct;

% sample parameters
fprintf('collecting shape parameters...\n');
S.N = N;
S.length = length;
S.a_0 = a_0;
S.size = size;
S.p_sam = p_sam;
S.name = name;
fprintf('...done\n');

% select sample shape
fprintf('creating shape...\n');
S.SS_shape = SS_shape;
fprintf('...done\n');

% crystal orientation
fprintf('gathering instrument information...\n');
S.hkl = hkl;
S.U_0 = U_0;

% 34-ID-C parameters
S.d = d;
S.lambda = lambda;
S.D = D;
S.rocking_angle = rocking_angle;
S.p_sam_REC = p_sam_REC; % a Command Window Output from reconstruction, in m (not nm). Comment/leave null if unknown/not needed 
fprintf('...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating diffraction beams and angles for APS geometry
fprintf('\n\n<<<PERFORMING BEAMLINE-SPECIFIC CALCULATIONS>>>\n\n');

% calculating diffraction angles based on chosen reflection
[S.theta_spec, S.chi_spec, S.phi_spec, S.delta_spec, S.gamma_spec, S.Q_lab, S.S_lab, S.S_0lab] = motor_angles_APS(S); % returns as SPEC angles: gamma_spec as -gamma and chi_spec as 90-chi

% calculating rocking increment angle and creating matrices
if strcmp(S.rocking_angle, 'dtheta')
    S.dtheta = rocking_increment_APS(S);
    [S.R_dqp_12, S.R_dqp_3, S.R_xyz, S.S_0lab_dir] = plugin_APS_34IDC(S.theta_spec, S.chi_spec, S.phi_spec, S.delta_spec, S.gamma_spec, S.dtheta, S.rocking_angle);
elseif strcmp(S.rocking_angle, 'dphi')
    S.dphi = rocking_increment_APS(S);
    [S.R_dqp_12, S.R_dqp_3, S.R_xyz, S.S_0lab_dir] = plugin_APS_34IDC(S.theta_spec, S.chi_spec, S.phi_spec, S.delta_spec, S.gamma_spec, S.dphi, S.rocking_angle);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating coordinates in lab, reciprocal lab and detector reciprocal space
fprintf('\n\n<<<MAKING COORDINATES>>>\n\n');
[S.N1grid, S.N2grid, S.N3grid] = orthogonal_coordinates(S.SS_shape);
[S.N1gridp, S.N2gridp, S.N3gridp] = non_orthogonal_coordinates(S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining shape in various spaces
fprintf('\n\n<<<PERFORMING SPACE TRANSFORMATIONS>>>\n\n');

% calculating the lab space shape from sample space
S.LS_shape_SS = SS_to_LS(S);

% calculating the reciprocal lab space shape from lab space
[S.RLS_shape_LS, S.I_RLS_shape_LS] = LS_to_RLS(S);

% calculating the detector reciprocal space shape from reciprocal lab space
S.DRS_shape_RLS = RLS_to_DRS(S);

% calculating the detector conjugated shape from detector reciprocal space 
S.DCS_shape_DRS = DRS_to_DCS(S);

% calculating the real space shape from detector conjugated shape
S.LS_shape_DCS = DCS_to_LS(S);

% calculating the detector conjugated shape from real space shape
S.DCS_shape_LS = LS_to_DCS(S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving -SAM, -LAB, -AMP, -PH, -SUP (binarized DCS shape) and TIFF files
if save_files == 1
    save_arrays_and_tiff(S, save_files_dir, binary_threshold);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
fprintf('\n\n<<<PLOTS>>>\n\n');
% plotting sample space (SS)
if plot_SS_figures == 1
    if exist('dir', 'var') && exist('folder_name', 'var')
        temp_dir_LAB = [dir, '/', folder_name, '/', folder_name,'-LAB.mat'];
    elseif exist('dir', 'var')
        folder_name = '';
        temp_dir_LAB = [dir, '/', folder_name, '/', folder_name,'-LAB.mat'];
    elseif exist('folder_name', 'var')
        temp_dir_LAB = [folder_name, '/', folder_name,'-LAB.mat'];
    else
        temp_dir_LAB = '';
    end
    plot_SS(S, plot_SS_axes, plot_SS_beams, plot_SS_grid, viewpoint, temp_dir_LAB, mask_threshold, plot_threshold, twin, SS_shape_plot, SS_shape_LS_plot, SS_shape_REC_plot);
end

% plotting lab space (LS)
if plot_LS_figures == 1
    if exist('dir', 'var') && exist('folder_name', 'var')
        temp_dir_LAB = [dir, '/', folder_name, '/', folder_name,'-LAB.mat'];
    elseif exist('dir', 'var')
        folder_name = '';
        temp_dir_LAB = [dir, '/', folder_name, '/', folder_name,'-LAB.mat'];
    elseif exist('folder_name', 'var')
        temp_dir_LAB = [folder_name, '/', folder_name,'-LAB.mat'];
    else
        temp_dir_LAB = '';
    end
    plot_LS(S, plot_LS_axes, plot_LS_beams, plot_LS_grid, viewpoint, temp_dir_LAB, mask_threshold, plot_threshold, twin, LS_shape_SS_plot, LS_shape_DCS_plot, LS_shape_REC_plot);
end

% plotting reciprocal lab space (RLS)
if plot_RLS_figures == 1
    plot_RLS(S, plot_RLS_axes, plot_RLS_beams, plot_RLS_grid, viewpoint, plot_RLS_log10scale);
end

% plotting detector reciprocal space (DRS)
if plot_DRS_figures == 1
	plot_DRS(S, plot_DRS_axes, plot_DRS_beams, plot_DRS_grid, viewpoint, plot_DRS_log10scale);
end

% plotting detector conjugated space (DCS)
if plot_DCS_figures == 1
    if exist('dir', 'var') && exist('folder_name', 'var')
        temp_dir_AMP = [dir, '/', folder_name, '/', folder_name,'-AMP.mat'];
        temp_dir_PH = [dir, '/', folder_name, '/', folder_name,'-PH.mat'];
    elseif exist('dir', 'var')
        folder_name = '';
        temp_dir_AMP = [dir, '/', folder_name, '/', folder_name,'-AMP.mat'];
        temp_dir_PH = [dir, '/', folder_name, '/', folder_name,'-PH.mat'];    
    elseif exist('folder_name', 'var')
        temp_dir_AMP = [folder_name, '/', folder_name,'-AMP.mat'];
        temp_dir_PH = [folder_name, '/', folder_name,'-PH.mat'];  
    else
        temp_dir_AMP = '';
        temp_dir_PH = '';
    end
    plot_DCS(S, plot_DCS_axes, plot_DCS_beams, plot_DCS_grid, viewpoint, temp_dir_AMP, temp_dir_PH, mask_threshold, plot_threshold, twin, DCS_shape_DRS_plot, DCS_shape_LS_plot, DCS_shape_REC_plot);
end

fprintf('\n\n<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
S
fprintf('\n                     FINISHED SIMULATION\n');
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
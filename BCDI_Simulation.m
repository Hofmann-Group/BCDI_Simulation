%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCDI Simulation
% Version 1.0
% May 2019
% Written By: David Yang
% University of Oxford, Dept. of Engineering Science, Hofmann Group
% 
% PURPOSE: To simulate a BCDI experiment for a given shape and reflection
% 
% USER-DEFINED SECTIONS:
% 1. Sample Space Shape Details
% - collects size and shape details about the sample space (SS) object
% 
% 2. BCDI Measurement Variables
% - user specifies the reflection and beamline measurement details 
% 
% 3. File Saving options
% - SAM.mat: sample space shape
% - LAB.mat: lab space shape
% - AMP.mat: amplitude of detector conjugated shape
% - PH.mat: phase of detector conjugated shape
% - BIN.mat: binarized detector conjugated shape
% - TIFF: slices through the bragg peak
% 
% 4. Plot Options
% - viewpoint
% - shapes to plot
% - log10
% - axes
% - grid
% - diffraction beams
% - mask and plot thresholds
% 
% 5. BCDI Space Plots
% - plot_SS: plots sample space
% - plot_LS: plots lab space
% - plot_RLS: plots reciprocal lab space
% - plot_DRS: plots detector reciprocal space
% - plot_DCS: plots detector conjugated space
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
fprintf('                         BCDI Simulation\n');
fprintf('                            Version 1.0\n');
fprintf('                             May 2019\n');
fprintf('                      Written By: David Yang\n');
fprintf(' University of Oxford, Dept. of Engineering Science, Hofmann Group\n');
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Sample Space Shape Details
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAMPLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.N = 2^8; % square detector length in number of pixels, should ideally be a power of 2
S.length = 2^6; % number of pixels that the object occupies, should ideally be <= N/4
S.a_0 = 3.16522e-10; % lattice constant for tungsten in m (arbitrarily chosen)
S.size = 1000*10^-9; % length of the cube that surrounds the object in sample space in m
S.p_sam = S.size/S.length; % pixel size of object in sample space in m
S.name = 'Cylinder'; % name of the shape (it will be on all files generated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x_pix, y_pix, z_pix] = meshgrid(-(S.N-1)/2:(S.N-1)/2, -(S.N-1)/2:(S.N-1)/2, -(S.N-1)/2:(S.N-1)/2); % creating a meshgrid for SS coordinates
S.SS_shape = zeros(S.N, S.N, S.N); % creating empty shape array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% SELECT SAMPLE SHAPE %%%%%%%%%%%%%%%%%%%%%%%%%%%
% user required to select shape by commenting/uncommenting the predefined shapes
S.SS_shape((x_pix.^2 + z_pix.^2) < (S.length/2)^2 & abs(y_pix) < S.length/2) = complex(1,0); % cylinder
% S.SS_shape(abs(x_pix) < S.length & abs(y_pix) < S.length & abs(z_pix) < S.length) = complex(1,0); % cube
% S.SS_shape((abs(x_pix) + abs(y_pix) + abs(z_pix)) < S.length) = complex(1,0); % octahedral
% S.SS_shape(((x_pix/2).^2 + (y_pix/1).^2 + (z_pix/3).^2) < (S.length/2)^2) = complex(1,0); % elippsoid
% S.SS_shape((x_pix.^2 + z_pix.^2) < (S.length/2)^2 & abs(y_pix) < S.length/2 & x_pix < 0) = complex(1,0); % semi-cylinder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2. BCDI Measurement Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%% CRYSTAL ORIENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
S.hkl = [-1; 2; 0]; % pick a reflection
S.U_0 = eye(3); % 3x3 orientation matrix, otherwise set to eye(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% 34-ID-C PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
S.d = 55e-6; % detector pixel size in m
S.lambda = (12.398/10.0)/10*10^-9; % wavelength in m
S.D = S.p_sam*S.N*S.d/S.lambda; % detector distance in m, put absolute value if known
S.rocking_axis = 'dtheta'; % choose 'dtheta' or 'dphi' rocking
S.rocking_increment = 0; % set to 0 for automatic calculation based setting magnitude of q'_3 to equal the magnitude of q'_1 or 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3. Saving file options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_files = 1; % 1 to save -SAM, -LAB, -AMP, -PH, -BIN and TIFF files
save_files_dir = 'Simulated Data'; % the folder to put the saved files
binary_threshold = 0.5; % binary threshold applied to abs(DCS_shape_DRS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 4. Plot Configurations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT VIEWPOINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options.viewpoint = [-180, -90]; % viewpoint = [az, el], x-y plane
% Options.viewpoint = [-180, 0]; % viewpoint = [az, el], x-z plane
% Options.viewpoint = [-90, 0]; % viewpoint = [az, el], y-z plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options.axes = 1; % 1 to show axes in the figure
Options.beams = 1; % 1 to show S_0lab, S_lab and Q_lab in the figure

Options.grid = 1; % 1 to show gridlines in the figure
Options.log10scale = 1; % 1 for log10 intensity, 0 for absolute intensity

Options.plot_SIM_fwd = 1; % 1 to show simulated shape via forward calculation
Options.plot_SIM_rev = 1; % 1 to show simulated shape via reverse calculation 

Options.mask_threshold = 0.3; % for creating binary masks from thresholded amplitude to calcuate centre of mass
Options.plot_threshold = 0.3; % amplitude threshold for figure plots

Options.plot_REC = 1; % 1 to show shape calculated from reconstruction
Options.twin = 0; % 0 or 1 to flip LAB shape, toggle to take the conjugate reflection of the reconstructed shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% RECONSTRUCTION PLOT OPTIONS %%%%%%%%%%%%%%%%%%%%%%%
Options.dir = 'Reconstruction Examples';

Options.folder_name = 'Rec-Cylinder_(-120)_dtheta-00274-ERlrHIOlr2000-NM-SW'; % reconstruction folder (which should be the root of the individual file names)

S.p_sam_REC = S.p_sam; % a Command Window Output from reconstruction, in m (not nm). Comment/leave null if unknown/not needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 5. BCDI Space Plots
%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SAMPLE SPACE (SS) %%%%%%%%%%%%%%%%%%%%%%%%%%
plot_SS_figures = 1; % 1 to show sample space figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT LAB SPACE (LS) %%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_LS_figures = 1; % 1 to show lab space figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%% PLOT RECIPROCAL LAB SPACE (RLS) %%%%%%%%%%%%%%%%%%%%%
plot_RLS_figures = 1; % 1 to show reciprocal lab space figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% PLOT DETECTOR RECIPROCAL SPACE (DRS) %%%%%%%%%%%%%%%%%%
plot_DRS_figures = 1; % 1 to show detector reciprocal space figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% PLOT DETECTOR CONJUGATED SPACE (DCS) %%%%%%%%%%%%%%%%%%
plot_DCS_figures = 1; % 1 to show detector conjugated space figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating paths to scripts in folders
addpath(genpath('APS 34-ID-C angles'));
addpath(genpath('Coordinates'));
addpath(genpath('Basis Vectors'));
addpath(genpath('Outputs'));
addpath(genpath('Reconstruction Examples'));
addpath(genpath('Simulation Space Transformations'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating diffraction beams and angles for APS geometry
fprintf('\n\n<<<PERFORMING BEAMLINE-SPECIFIC CALCULATIONS>>>\n\n');

% calculating diffraction angles based on chosen reflection
[S.theta_bl, S.chi_bl, S.phi_bl, S.delta_bl, S.gamma_bl, S.Q_lab, S.S_lab, S.S_0lab] = motor_angles_APS_34IDC(S, 0); % returns as beamline angles. gamma_bl as -gamma and chi_bl as 90-chi

% calculating rocking increment (such that |q'_1|=|q'_2|=|q'_3|) and creating fundamental matrices
if strcmp(S.rocking_axis, 'dtheta')
    if S.rocking_increment == 0
        S.dtheta = rocking_increment_APS_34IDC(S);
    else
        S.dtheta = S.rocking_increment;
    end
    [S.R_dqp_12, S.R_dqp_3, S.R_xyz, S.S_0lab_dir] = plugin_APS_34IDC(S.theta_bl, S.chi_bl, S.phi_bl, S.delta_bl, S.gamma_bl, S.dtheta, S.rocking_axis);
elseif strcmp(S.rocking_axis, 'dphi')
    if S.rocking_increment == 0
        S.dphi = rocking_increment_APS_34IDC(S);
    else
        S.dphi = S.rocking_increment;
    end
    [S.R_dqp_12, S.R_dqp_3, S.R_xyz, S.S_0lab_dir] = plugin_APS_34IDC(S.theta_bl, S.chi_bl, S.phi_bl, S.delta_bl, S.gamma_bl, S.dphi, S.rocking_axis);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating orthogonal and non-orthogonal coordinates
fprintf('\n\n<<<MAKING COORDINATES>>>\n\n');
[S.N1grid, S.N2grid, S.N3grid] = orthogonal_coordinates(S.SS_shape);
[S.N1gridp, S.N2gridp, S.N3gridp] = non_orthogonal_coordinates(S);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining basis vectors in various spaces
fprintf('\n\n<<<MAKING BASIS VECTORS>>>\n\n');
[S.x_sam, S.y_sam, S.z_sam] = SS_basis_vectors(S);
[S.x, S.y, S.z] = LS_basis_vectors(S);
[S.q_x, S.q_y, S.q_z] = RLS_basis_vectors(S);
[S.q_1p, S.q_2p, S.q_3p] = DRS_basis_vectors(S);
[S.xp, S.yp, S.zp] = DCS_basis_vectors(S);


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
%% Saving -SAM, -LAB, -AMP, -PH, -BIN and TIFF files
if save_files == 1
    save_arrays_and_tiff(S, save_files_dir, binary_threshold);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots
fprintf('\n\n<<<PLOTS>>>\n\n');
% plotting sample space (SS)
if plot_SS_figures == 1
    plot_SS(S, Options);
end

% plotting lab space (LS)
if plot_LS_figures == 1
    plot_LS(S, Options);
end

% plotting reciprocal lab space (RLS)
if plot_RLS_figures == 1
    plot_RLS(S, Options);
end

% plotting detector reciprocal space (DRS)
if plot_DRS_figures == 1
	plot_DRS(S, Options);
end

% plotting detector conjugated space (DCS)
if plot_DCS_figures == 1
    plot_DCS(S, Options);
end


fprintf('\n\n<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
S
fprintf('\n                        FINISHED SIMULATION\n');
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
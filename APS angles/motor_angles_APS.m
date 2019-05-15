function [theta_spec, chi_spec, phi_spec, delta_spec, gamma_spec, Q_lab, S_lab, S_0lab] = motor_angles_APS(S, varargin)
fprintf('motor_angles_APS');
%% Lattice spacings
% real space lattice spacing
a_0 = S.a_0; % Sample lattice parameter

% reciprocal space lattice vector matrix
B_0 = [2*pi/a_0 0 0; 0 2*pi/a_0 0; 0 0 2*pi/a_0];

% orientation matrix
U_0 = S.U_0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sample space 
% Q vector in sample coordinates
hkl = S.hkl;
Q_sam = U_0*B_0*hkl; 
Q_sam_mag = norm(Q_sam);
Q_sam_n = Q_sam/Q_sam_mag;

% normal vector to substrate in sample coordinates
n_sam = [0; 1; 0]; % 34-ID-C (APS) convention
% n_sam = [0; 0; 0]; % ID01 (ESRF) convention

% angle between Q_sam and S_0sam:
lambda = S.lambda;
ang_Q_sam_S_0sam = real(acosd(Q_sam_mag/(2*2*pi/lambda))); % why is this 2*(2*pi/lambda)? ans: use Bragg's law and an sine angle sum identity

% unit normal to Q_sam (Q vector) and n_sam_n (sample surface normal)
ip_sam = cross(Q_sam, n_sam); 

% check in case the Q vector is aligned with the sample normal
if ip_sam == 0
    ip_sam = [1; 0; 0];
end
ip_sam_n = ip_sam/norm(ip_sam);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the first set of incident and exit beams
fprintf('\n...calculating first set of incident and exit beams...');
% Sample coordinates
% compute S_0 in sample coordiantes using Q_sam_n and ip_sam_n:
S_0sam_temp_1 = -(tand(ang_Q_sam_S_0sam)*ip_sam_n + Q_sam_n);
S_0sam_n_1 = S_0sam_temp_1/norm(S_0sam_temp_1); % unit vector along S_0 direction
S_0sam_1 = S_0sam_n_1*(2*pi/lambda); % S_0 vector in sample coordinates

% now compute s in sample coordiantes using Q_sam and S_0sam
S_sam_1 = S_0sam_1 + Q_sam; % this is the exit beam vector
S_sam_n_1 = S_sam_1/norm(S_sam_1); % unit vector along s direction


% Lab coordinates
% compute phi, chi and theta, based on projection of S onto xz plane
phi_1 = 2*atand(S_0sam_n_1(2)/S_0sam_n_1(3)); % solved rotation matrices to get phi and chi based on the constraint that S_0sam projected onto the xz plane should be parallel to [a; 0; b]
chi_1 = atand((S_0sam_n_1(3)*sind(phi_1)-S_0sam_n_1(2)*cosd(phi_1))/S_0sam_n_1(1));
R_chi_phi_1 = rotz(chi_1)*rotx(phi_1);
S_0sam_1_proj_temp = R_chi_phi_1*S_0sam_1; % S_0sam projected onto the xz plane, parallel to [a; 0; b]
theta_1 = atand(-S_0sam_1_proj_temp(1)/S_0sam_1_proj_temp(3));

R_xyz_1 = roty(theta_1)*rotz(chi_1)*rotx(phi_1); % sample rotation matrix
S_0lab_1 = R_xyz_1*S_0sam_1; % should be parallel to [0; 0; 1];
S_lab_1 = R_xyz_1*S_sam_1;
Q_lab_1 = R_xyz_1*Q_sam;

% Calculate detector angles
S_0lab_n_1 = R_xyz_1*S_sam_n_1; % exit beams for option 1

gamma_1 = asind(-S_0lab_n_1(2)); % rotate about x
delta_1 = asind(S_0lab_n_1(1)/cosd(gamma_1)); % rotate about y

det_1 = roty(delta_1)*rotx(gamma_1)*[0; 0; 1]; % should be the same direction as R_xyz*S_sam_n_1 = S_lab

% fixing the dectector vector orientation and delta angle
if sign(det_1(3))~=sign(S_0lab_n_1(3))
    delta_1 = 180-delta_1;
    det_1 = roty(delta_1)*rotx(gamma_1)*[0; 0; 1];
%     disp('changed delta_1');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the second set of incident and exit beams - reciprocal of set 1.
fprintf('\n...calculating second set of incident and exit beams...');
% Sample Coordinates
% compute the second set of S_0sam and S in sample coordinates using the reciprocal of set 1:
S_0sam_2 = -S_sam_1; % this is the incident beam vector
S_sam_2 = -S_0sam_1; % this is the exit beam vector

% now compute unit vectors
S_0sam_n_2 = S_0sam_2/norm(S_0sam_2); % unit vector along S_0sam direction
S_sam_n_2 = S_sam_2/norm(S_sam_2); % unit vector along s direction


% Lab coordinates
% compute phi, chi and theta, based on projection of s onto xz plane
phi_2 = 2*atand(S_0sam_n_2(2)/S_0sam_n_2(3)); % solved rotation matrices to get phi and chi based on the constraint that S_0sam projected onto the xz plane should be parallel to [a; 0; b]
chi_2 = atand((S_0sam_n_2(3)*sind(phi_2)-S_0sam_n_2(2)*cosd(phi_2))/S_0sam_n_2(1));
R_chi_phi_2 = rotz(chi_2)*rotx(phi_2); 
S_0sam_2_proj_temp = R_chi_phi_2*S_0sam_2; % S_0sam projected onto the xz plane, parallel to [a; 0; b]
theta_2 = atand(-S_0sam_2_proj_temp(1)/S_0sam_2_proj_temp(3));

R_xyz_2 = roty(theta_2)*rotz(chi_2)*rotx(phi_2); % sample rotation matrix
S_0lab_2 = R_xyz_2*S_0sam_2; % should be parallel to [0; 0; 1];
S_lab_2 = R_xyz_2*S_sam_2;
Q_lab_2 = R_xyz_2*Q_sam;

% Calculate the detector angles
S_0lab_n_2 = R_xyz_2*S_sam_n_2; % exit beams for option 2: incident and exit beams from option 1 swapped - checked, this is true in sam

gamma_2 = asind(-S_0lab_n_2(2)); % rotate about x
delta_2 = asind(S_0lab_n_2(1)/cosd(gamma_2)); % rotate about y

det_2 = roty(delta_2)*rotx(gamma_2)*[0; 0; 2*pi/lambda]; % should be the same direction as R_xyz*S_sam_n_2 = S_lab

% fixing the dectector vector orientation and delta angle
if sign(det_2(3))~=sign(S_0lab_n_2(3))
    delta_2 = 180-delta_2;
    det_2 = roty(delta_2)*rotx(gamma_2)*[0; 0; 1];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checks from Felix's Code
% work out the other combination for delta and gamma that could be used if
% abs(delta)>90 and abs(gamma)>90

gamma_1_opt2 = nan;
delta_1_opt2 = nan;
gamma_2_opt2 = nan;
delta_2_opt2 = nan;

if abs(delta_1)>90
    if sign(delta_1)==1
        delta_1_opt2 = delta_1 - 180;
    elseif sign(delta_1)==-1
        delta_1_opt2 = delta_1 + 180;
    end

    if sign(gamma_1)==1
        gamma_1_opt2 = 180 - gamma_1;
    elseif sign(gamma_1)==-1
        gamma_1_opt2 = -180 - gamma_1;
    end
end

if abs(delta_2)>90
    if sign(delta_2)==1
        delta_2_opt2 = delta_2 - 180;
    elseif sign(delta_2)==-1
        delta_2_opt2 = delta_2 + 180;
    end

    if sign(gamma_2)==1
        gamma_2_opt2 = 180 - gamma_2;
    elseif sign(gamma_2)==-1
        gamma_2_opt2 = -180 - gamma_2;
    end
end

% can display these if you want
% delta_1_opt2
% gamma_1_opt2
% 
% delta_2_opt2
% gamma_2_opt2
        
% check that the rotation matrices are correct...
% R_dqp_12_1 = roty(delta_1)*rotx(gamma_1);
% det_01_n_1 = R_dqp_12_1'*S_0lab_n_1;  % takes diffracted beam in lab coordinates and rotates back to 0 detector angles to make sure this lines up with incindent beam
% 
% R_dqp_12_1_opt2 = roty(delta_1_opt2)*rotx(gamma_1_opt2);
% det_01_n_1_opt2 = R_dqp_12_1_opt2'*S_0lab_n_1;
% 
% R_dqp_12_2 = roty(delta_2)*rotx(gamma_2);
% det_01_n_2 = R_dqp_12_2'*S_0lab_n_2;
% 
% R_dqp_12_2_opt2 = roty(delta_2_opt2)*rotx(gamma_2_opt2);
% det_01_n_2_opt2 = R_dqp_12_2_opt2'*S_0lab_n_2;


% det_01_n_1 % should all be [0; 0; 1]
% det_01_n_1_opt2
% det_01_n_2
% det_01_n_2_opt2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choosing which set to use
if sign(S_0lab_1(3))==1
    theta_spec = theta_1;
    chi_spec = 90-chi_1; % spec convention
    phi_spec = phi_1;
    delta_spec = delta_1;
    gamma_spec = -gamma_1; % spec convention
    Q_lab = Q_lab_1;
    S_lab = S_lab_1;
    S_0lab = S_0lab_1;
    fprintf('\n...used the first set...');
else
    theta_spec = theta_2;
    chi_spec = 90-chi_2; % spec convention
    phi_spec = phi_2;
    delta_spec = delta_2;
    gamma_spec = -gamma_2; % spec convention
    Q_lab = Q_lab_2;
    S_lab = S_lab_2;
    S_0lab = S_0lab_2;
    fprintf('\n...used the second set...');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Geometry Vector Plots
if nargin == 2 && varargin{1} == 1
    fprintf('\n...plotting geometry vector plots...');
    % First set of incident and exit beams (sample coordinates)
    figure('Name','Sample and Lab Coordinates (1st Set)');
    axes_1_1 = subplot(2,1,1);
    hold on;
    
    % diffraction beams
    beam_S_0lab = quiver3(-0.9*S_0sam_1(1),-0.9*S_0sam_1(2),-0.9*S_0sam_1(3),0.9*S_0sam_1(1),0.9*S_0sam_1(2),0.9*S_0sam_1(3));
    set(beam_S_0lab,'Color','blue','Linewidth', 2, 'AutoScale','off');
    text(-S_0sam_1(1),-S_0sam_1(2),-S_0sam_1(3),'S_{0lab}','Color','blue','FontSize',14);
    beam_S_lab = quiver3(0,0,0,0.9*S_sam_1(1),0.9*S_sam_1(2),0.9*S_sam_1(3));
    set(beam_S_lab,'Color','red','Linewidth', 2, 'AutoScale','off');
    text(S_sam_1(1),S_sam_1(2),S_sam_1(3),'S_{lab}','Color','red','FontSize',14);
    beam_Q_lab = quiver3(0,0,0,0.9*Q_sam(1),0.9*Q_sam(2),0.9*Q_sam(3));
    set(beam_Q_lab,'Color','green','Linewidth', 2, 'AutoScale','off');
    text(Q_sam(1),Q_sam(2),Q_sam(3),'Q_{lab}','Color','green','FontSize',14);
    
    % geometry vectors
    beam_n_sam = quiver3(0,0,0,0.9*n_sam(1)*pi/lambda/2,0.9*n_sam(2)*pi/lambda/2,0.9*n_sam(3)*pi/lambda/2);
    set(beam_n_sam,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(n_sam(1)*pi/lambda/2,n_sam(2)*pi/lambda/2,n_sam(3)*pi/lambda/2,'n_{sam}','Color','black','FontSize',14);
    beam_hkl = quiver3(0,0,0,0.9*hkl(1)*pi/lambda/2,0.9*hkl(2)*pi/lambda/2,0.9*hkl(3)*pi/lambda/2);
    set(beam_hkl,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(hkl(1)*pi/lambda/2,hkl(2)*pi/lambda/2,hkl(3)*pi/lambda/2,'hkl','Color','black','FontSize',14);
    
    % x, y, z axes
    x_axis = quiver3(0,0,0,0.9*2*pi/lambda,0,0);
    set(x_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(2*pi/lambda,0,0,'x','Color','black','FontSize',14);
    y_axis = quiver3(0,0,0,0,0.9*2*pi/lambda,0);
    set(y_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,2*pi/lambda,0,'y','Color','black','FontSize',14);
    z_axis = quiver3(0,0,0,0,0,0.9*2*pi/lambda);
    set(z_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,0,2*pi/lambda,'z','Color','black','FontSize',14);
    
    % plot parameters
    grid on;
    axis equal;
    view(-180, -90);
    title('Sample Coordinates (1st Set)');
    xlim([-6e10, 6e10]); ylim([-6e10, 6e10]); zlim([-6e10, 6e10]);
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % First set of incident and exit beams (lab coordinates)
    axes_1_2 = subplot(2,1,2);
    hold on;
    
    % diffraction beams
    beam_S_0lab = quiver3(-0.9*S_0lab_1(1),-0.9*S_0lab_1(2),-0.9*S_0lab_1(3),0.9*S_0lab_1(1),0.9*S_0lab_1(2),0.9*S_0lab_1(3));
    set(beam_S_0lab,'Color','blue','Linewidth', 2, 'AutoScale','off');
    text((-S_0lab_1(1)),(-S_0lab_1(2)),(-S_0lab_1(3)),'S_{0lab}','Color','blue','FontSize',14);
    beam_S_lab = quiver3(0,0,0,0.9*S_lab_1(1),0.9*S_lab_1(2),0.9*S_lab_1(3));
    set(beam_S_lab,'Color','red','Linewidth', 2, 'AutoScale','off');
    text((S_lab_1(1)),(S_lab_1(2)),(S_lab_1(3)),'S_{lab}','Color','red','FontSize',14);
    beam_Q_lab = quiver3(0,0,0,0.9*Q_lab_1(1),0.9*Q_lab_1(2),0.9*Q_lab_1(3));
    set(beam_Q_lab,'Color','green','Linewidth', 2, 'AutoScale','off');
    text((Q_lab_1(1)),(Q_lab_1(2)),(Q_lab_1(3)),'Q_{lab}','Color','green','FontSize',14);
    
    % geometry vectors
    beam_hkl = quiver3(0,0,0,0.9*hkl(1)*pi/lambda/2,0.9*hkl(2)*pi/lambda/2,0.9*hkl(3)*pi/lambda/2);
    set(beam_hkl,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((hkl(1)*pi/lambda/2),(hkl(2)*pi/lambda/2),(hkl(3)*pi/lambda/2),'hkl','Color','black','FontSize',14);
    
    % x, y, z axes
    x_axis = quiver3(0,0,0,0.9*2*pi/lambda,0,0);
    set(x_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(2*pi/lambda,0,0,'x','Color','black','FontSize',14);
    y_axis = quiver3(0,0,0,0,0.9*2*pi/lambda,0);
    set(y_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,2*pi/lambda,0,'y','Color','black','FontSize',14);
    z_axis = quiver3(0,0,0,0,0,0.9*2*pi/lambda);
    set(z_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,0,2*pi/lambda,'z','Color','black','FontSize',14);
    
    % plot parameters
    grid on;
    axis equal;
    view(-180, -90);
    title('Lab Coordinates (1st Set)');
    xlim([-6e10, 6e10]); ylim([-6e10, 6e10]); zlim([-6e10, 6e10]);
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % detector vector (should overlap S_lab)
    beam_det_1 = quiver3(0,0,0,0.9*det_1(1)*2*pi/lambda,0.9*det_1(2)*2*pi/lambda,0.9*det_1(3)*2*pi/lambda);
    set(beam_det_1,'Color','yellow','Linewidth', 2, 'AutoScale','off');
    text(det_1(1)*2*pi/lambda,det_1(2)*2*pi/lambda,det_1(3)*2*pi/lambda,'det_1','Color','yellow','FontSize',14);
    
    % linking plot axes for first set of beams
    axislink_1 = linkprop([axes_1_1,axes_1_2],{'CameraPosition','CameraUpVector'});
    setappdata(gcf, 'StoreTheLink', axislink_1);
    
    
    % Second set of incident and exit beams (sample coordinates)
    figure('Name','Sample and Lab Coordinates (2nd Set)');
    axes_2_1 = subplot(2,1,1);
    hold on;
    
    % diffraction beams
    beam_S_0lab = quiver3(-0.9*S_0sam_2(1),-0.9*S_0sam_2(2),-0.9*S_0sam_2(3),0.9*S_0sam_2(1),0.9*S_0sam_2(2),0.9*S_0sam_2(3));
    set(beam_S_0lab,'Color','blue','Linewidth', 2, 'AutoScale','off');
    text(-S_0sam_2(1),-S_0sam_2(2),-S_0sam_2(3),'S_{0lab}','Color','blue','FontSize',14);
    beam_S_lab = quiver3(0,0,0,0.9*S_sam_2(1),0.9*S_sam_2(2),0.9*S_sam_2(3));
    set(beam_S_lab,'Color','red','Linewidth', 2, 'AutoScale','off');
    text(S_sam_2(1),S_sam_2(2),S_sam_2(3),'S_{lab}','Color','red','FontSize',14);
    beam_Q_lab = quiver3(0,0,0,0.9*Q_sam(1),0.9*Q_sam(2),0.9*Q_sam(3));
    set(beam_Q_lab,'Color','green','Linewidth', 2, 'AutoScale','off');
    text((Q_sam(1)),(Q_sam(2)),(Q_sam(3)),'Q_{lab}','Color','green','FontSize',14);
    
    % geometry vectors
    beam_n_sam = quiver3(0,0,0,0.9*n_sam(1)*pi/lambda/2,0.9*n_sam(2)*pi/lambda/2,0.9*n_sam(3)*pi/lambda/2);
    set(beam_n_sam,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((n_sam(1)*pi/lambda/2),(n_sam(2)*pi/lambda/2),(n_sam(3)*pi/lambda/2),'n_s_c','Color','black','FontSize',14);
    beam_hkl = quiver3(0,0,0,0.9*hkl(1)*pi/lambda/2,0.9*hkl(2)*pi/lambda/2,0.9*hkl(3)*pi/lambda/2);
    set(beam_hkl,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((hkl(1)*pi/lambda/2),(hkl(2)*pi/lambda/2),(hkl(3)*pi/lambda/2),'hkl','Color','black','FontSize',14);
    
    % x, y, z axes
    x_axis = quiver3(0,0,0,0.9*2*pi/lambda,0,0);
    set(x_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(2*pi/lambda,0,0,'x','Color','black','FontSize',14);
    y_axis = quiver3(0,0,0,0,0.9*2*pi/lambda,0);
    set(y_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,2*pi/lambda,0,'y','Color','black','FontSize',14);
    z_axis = quiver3(0,0,0,0,0,0.9*2*pi/lambda);
    set(z_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,0,2*pi/lambda,'z','Color','black','FontSize',14);
    % % Projection of S_0sam_n_1 onto the xz plane
    % S_0sam_2_proj = S_0sam_2-dot(S_0sam_2,[0 1 0]')*[0 1 0]';
    % S_0sam_n_2_proj = S_0sam_2_proj/norm(S_0sam_2_proj);
    % beam_proj = quiver3((-S_0sam_2_proj(1)),(-S_0sam_2_proj(2)),(-S_0sam_2_proj(3)),(S_0sam_2_proj(1)),(S_0sam_2_proj(2)),(S_0sam_2_proj(3))); % last three indices are the S vector direction
    % set(beam_proj,'Color','red','Linewidth', 2, 'AutoScale','off');
    % text((-S_0sam_2_proj(1)),(-S_0sam_2_proj(2)),(-S_0sam_2_proj(3)),'S_{0proj}','Color','red','FontSize',14);
    
    % plot parameters
    grid on;
    axis equal;
    view(-180, -90);
    title('Sample Coordinates (2nd Set)');
    xlim([-6e10, 6e10]); ylim([-6e10, 6e10]); zlim([-6e10, 6e10]);
    xlabel('x'); ylabel('y'); zlabel('z');

    % Second set of incident and exit beams (lab coordinates)
    axes_2_2 = subplot(2,1,2);
    hold on;
    
    % diffraction beams
    beam_S_0lab = quiver3(-0.9*S_0lab_2(1),-0.9*S_0lab_2(2),-0.9*S_0lab_2(3),0.9*S_0lab_2(1),0.9*S_0lab_2(2),0.9*S_0lab_2(3));
    set(beam_S_0lab,'Color','blue','Linewidth', 2, 'AutoScale','off');
    text((-S_0lab_2(1)),(-S_0lab_2(2)),(-S_0lab_2(3)),'S_{0lab}','Color','blue','FontSize',14);
    beam_S_lab = quiver3(0,0,0,0.9*S_lab_2(1),0.9*S_lab_2(2),0.9*S_lab_2(3));
    set(beam_S_lab,'Color','red','Linewidth', 2, 'AutoScale','off');
    text((S_lab_2(1)),(S_lab_2(2)),(S_lab_2(3)),'S_{lab}','Color','red','FontSize',14);
    beam_Q_lab = quiver3(0,0,0,0.9*Q_lab_2(1),0.9*Q_lab_2(2),0.9*Q_lab_2(3));
    set(beam_Q_lab,'Color','green','Linewidth', 2, 'AutoScale','off');
    text(Q_lab_2(1),Q_lab_2(2),Q_lab_2(3),'Q_{lab}','Color','green','FontSize',14);
    
    % geometry vectors
    beam_hkl = quiver3(0,0,0,0.9*hkl(1)*pi/lambda/2,0.9*hkl(2)*pi/lambda/2,0.9*hkl(3)*pi/lambda/2);
    set(beam_hkl,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((hkl(1)*pi/lambda/2),(hkl(2)*pi/lambda/2),(hkl(3)*pi/lambda/2),'hkl','Color','black','FontSize',14);
    
    % x, y, z axes
    x_axis = quiver3(0,0,0,0.9*2*pi/lambda,0,0);
    set(x_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(2*pi/lambda,0,0,'x','Color','black','FontSize',14);
    y_axis = quiver3(0,0,0,0,0.9*2*pi/lambda,0);
    set(y_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,2*pi/lambda,0,'y','Color','black','FontSize',14);
    z_axis = quiver3(0,0,0,0,0,0.9*2*pi/lambda);
    set(z_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,0,2*pi/lambda,'z','Color','black','FontSize',14);
    
    % plot parameters
    grid on;
    axis equal;
    view(-180, -90);
    title('Lab Coordinates (2nd Set)');
    xlim([-6e10, 6e10]); ylim([-6e10, 6e10]); zlim([-6e10, 6e10]);
    xlabel('x'); ylabel('y'); zlabel('z');
    
    % detector vector (should overlap S_lab)
    beam_det_2 = quiver3(0,0,0,0.9*det_2(1)*2*pi/lambda,0.9*det_2(2)*2*pi/lambda,0.9*det_2(3)*2*pi/lambda); % last three indices belong to the S vector direction
    set(beam_det_2,'Color','yellow','Linewidth', 2, 'AutoScale','off');
    text(det_2(1)*2*pi/lambda,det_2(2)*2*pi/lambda,det_2(3)*2*pi/lambda,'det_2','Color','yellow','FontSize',14);
    
    % linking plot axes for first set of beams
    axislink_2 = linkprop([axes_2_1,axes_2_2],{'CameraPosition','CameraUpVector'});
    setappdata(gcf, 'StoreTheLink', axislink_2);
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
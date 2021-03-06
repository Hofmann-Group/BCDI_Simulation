function plot_RLS(S, Options)
% plots the reciprocal lab space
fprintf('plot_RLS');
%% Plotting reciprocal lab space shape from lab space
fprintf('\n...plotting RLS intensity from LS...');
% applying log10 scale if necessary
if Options.log10scale == 1
    S.I_RLS_shape_LS = log10(S.I_RLS_shape_LS);
end

% creating plot
figure('Name','Reciprocal Lab Space (RLS)');
hold on;

% plot grids
N1grid = S.N1grid*2*pi/S.lambda*S.d/S.D;
N2grid = S.N2grid*2*pi/S.lambda*S.d/S.D;
N3grid = S.N3grid*2*pi/S.lambda*S.d/S.D;

% intensity in lab reciprocal space
plot_I_RLS_shape_LS = patch(isosurface(N1grid,N2grid,N3grid,S.I_RLS_shape_LS,6));
set(plot_I_RLS_shape_LS, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha',1);

% scaling
FT_scale = 1/S.N;

% % diffraction beams
% if Options.beams == 1
%     fprintf('\n...plotting diffraction beams...');
%     beam_S_0lab = quiver3(0,0,-0.9*norm(S.S_0lab)*FT_scale,0.9*(S.S_0lab(1,1)*FT_scale),0.9*(S.S_0lab(2,1)*FT_scale),0.9*(S.S_0lab(3,1)*FT_scale));
%     set(beam_S_0lab,'Color','blue','Linewidth',2,'MaxHeadSize',0.5, 'AutoScale','off');
%     text(0,0,-norm(S.S_0lab)*FT_scale,'S_{0lab}','Color','blue','FontSize',14);
%     beam_S_lab = quiver3(0,0,0,0.9*(S.S_lab(1,1)*FT_scale),0.9*(S.S_lab(2,1)*FT_scale),0.9*(S.S_lab(3,1)*FT_scale)); % last three indices are the s vector direction
%     set(beam_S_lab,'Color','red','Linewidth',2,'MaxHeadSize',0.5,'AutoScale','off');
%     text((S.S_lab(1,1)*FT_scale),(S.S_lab(2,1)*FT_scale),(S.S_lab(3,1)*FT_scale),'S_{lab}','Color','red','FontSize',14);
%     beam_Q_lab = quiver3(0,0,0,0.9*(S.Q_lab(1,1)*FT_scale),0.9*(S.Q_lab(2,1)*FT_scale),0.9*(S.Q_lab(3,1)*FT_scale)); % last three indices are the s vector direction
%     set(beam_Q_lab,'Color','green','Linewidth',2,'MaxHeadSize',0.5, 'AutoScale','off');
%     text((S.Q_lab(1,1)*FT_scale),(S.Q_lab(2,1)*FT_scale),(S.Q_lab(3,1)*FT_scale),'Q_{lab}','Color','green','FontSize',14);
% end
    
% q_x, q_y and q_z axes
if Options.axes == 1
    fprintf('\n...plotting axes...');
    q_x_axis = quiver3(0,0,0,(0.9*norm(S.S_0lab)*FT_scale),0,0); % last three indices are the s vector direction
    set(q_x_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((1*norm(S.S_0lab)*FT_scale),0,0,'q_x','Color','black','FontSize',14);
    q_y_axis = quiver3(0,0,0,0,(0.9*norm(S.S_0lab)*FT_scale),0); % last three indices are the s vector direction
    set(q_y_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,(1*norm(S.S_0lab)*FT_scale),0,'q_y','Color','black','FontSize',14);
    q_z_axis = quiver3(0,0,0,0,0,(0.9*norm(S.S_0lab)*FT_scale)); % last three indices are the s vector direction
    set(q_z_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text(0,0,(1*norm(S.S_0lab)*FT_scale),'q_z','Color','black','FontSize',14);
end

% plot parameters
title(['Reciprocal Lab Space Shape with detector at \gamma = ', num2str(-S.gamma_bl), char(176), ' and \delta = ', num2str(S.delta_bl), char(176)]); xlabel('q_x (m^-^1)'); ylabel('q_y (m^-^1)'); zlabel('q_z (m^-^1)');
set(gca,'XDir','normal');
set(gca,'YDir','normal');
daspect([1,1,1]);
axis equal;
axis vis3d xy;
view(Options.viewpoint(1), Options.viewpoint(2));
lighting gouraud;
camlight('headlight');
xlim([min(min(min(N1grid))) max(max(max(N1grid)))]); ylim([min(min(min(N2grid))) max(max(max(N2grid)))]); zlim([min(min(min(N3grid))) max(max(max(N3grid)))]);
if Options.grid == 1
    grid on;
end
legend(plot_I_RLS_shape_LS,'RLS shape from LS');
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function plot_DRS(S, Options)
% plots the detector reciprocal space
fprintf('plot_DRS');
%% Plotting detector reciprocal space from reciprocal lab space
fprintf('\n...plotting DRS intensity from RLS...');
% applying log10 scale if necessary
if Options.log10scale == 1
    S.DRS_shape_RLS = log10(S.DRS_shape_RLS);
end

% creating plot
figure('Name','Detector Reciprocal Space (DRS)');
hold on;

% plot grids
N1grid = S.N1grid*2*pi/S.lambda*S.d/S.D;
N2grid = S.N2grid*2*pi/S.lambda*S.d/S.D;
N3grid = S.N3grid*2*pi/S.lambda*S.d/S.D;

% intensity in detector reciprocal space
plot_DRS_shape_RLS = patch(isosurface(N1grid,N2grid,N3grid,S.DRS_shape_RLS,6));
set(plot_DRS_shape_RLS, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha',1);

% scaling
FT_scale = 1/S.N;

% diffraction beams
if Options.beams == 1
    fprintf('\n...plotting diffraction beams...');
    beam_S_0lab = quiver3(0,0,-0.9*norm(S.S_0lab)*FT_scale,0.9*(S.S_0lab(1,1)*FT_scale),0.9*(S.S_0lab(2,1)*FT_scale),0.9*(S.S_0lab(3,1)*FT_scale));
    set(beam_S_0lab,'Color','blue','Linewidth',2,'MaxHeadSize',0.5, 'AutoScale','off');
    text(0,0,-norm(S.S_0lab)*FT_scale,'S_{0lab}','Color','blue','FontSize',14);
    beam_S_lab = quiver3(0,0,0,0.9*(S.S_lab(1,1)*FT_scale),0.9*(S.S_lab(2,1)*FT_scale),0.9*(S.S_lab(3,1)*FT_scale)); % last three indices are the s vector direction
    set(beam_S_lab,'Color','red','Linewidth',2,'MaxHeadSize',0.5,'AutoScale','off');
    text((S.S_lab(1,1)*FT_scale),(S.S_lab(2,1)*FT_scale),(S.S_lab(3,1)*FT_scale),'S_{lab}','Color','red','FontSize',14);
    beam_Q_lab = quiver3(0,0,0,0.9*(S.Q_lab(1,1)*FT_scale),0.9*(S.Q_lab(2,1)*FT_scale),0.9*(S.Q_lab(3,1)*FT_scale)); % last three indices are the s vector direction
    set(beam_Q_lab,'Color','green','Linewidth',2,'MaxHeadSize',0.5, 'AutoScale','off');
    text((S.Q_lab(1,1)*FT_scale),(S.Q_lab(2,1)*FT_scale),(S.Q_lab(3,1)*FT_scale),'Q_{lab}','Color','green','FontSize',14);
end

% plot parameters
title(['Detector Reciprocal Space Shape with detector at \gamma = ', num2str(-S.gamma_bl), char(176), ' and \delta = ', num2str(S.delta_bl), char(176)]); xlabel('q''_1 (m^-^1)'); ylabel('q''_2 (m^-^1)'); zlabel('q''_3 (m^-^1)');
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
legend(plot_DRS_shape_RLS,'DRS shape from RLS')
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

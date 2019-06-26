function plot_LS(S, Options)
% plots the lab space
fprintf('plot_LS');
%% Plotting lab space shape from sample space 
fprintf('\n...plotting LS shape from SS...');
figure('Name','Lab Space (LS)');
hold on;

% shifting LS shape from SS to the centre of mass
fprintf('\n...shifting LS shape from SS to the centre of mass...');
LS_shape_SS_MASK = single(abs(S.LS_shape_SS) > Options.mask_threshold);
structure_element = strel('sphere', 3);
LS_shape_SS_MASK = imerode(imdilate(LS_shape_SS_MASK, structure_element),structure_element); % takes care of dislocation cores
LS_shape_SS_COM = ceil(centerOfMass(LS_shape_SS_MASK));
S.LS_shape_SS = circshift(S.LS_shape_SS, size(S.LS_shape_SS)/2-LS_shape_SS_COM);

% plot grids
N1grid = S.N1grid*S.p_sam;
N2grid = S.N2grid*S.p_sam;
N3grid = S.N3grid*S.p_sam;

% plotting LS shape from SS
if Options.plot_SIM_fwd == 1
    plot_LS_shape_SS = patch(isosurface(N1grid, N2grid, N3grid, S.LS_shape_SS));
    set(plot_LS_shape_SS, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    legend('LS shape from SS');
end

% plot parameters
title(['Lab Space Shape with detector at \gamma = ', num2str(-S.gamma_bl), char(176), ' and \delta = ', num2str(S.delta_bl), char(176)]); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); 
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

% scaling
scale = norm(S.S_lab)/S.N*2/(S.p_sam);

% diffraction beams
if Options.beams == 1
    fprintf('\n...plotting diffraction beams...');
    beam_S_0lab = quiver3(0,0,-0.9*norm(S.S_0lab)/scale,0.9*(S.S_0lab(1,1)/scale),0.9*(S.S_0lab(2,1)/scale),0.9*(S.S_0lab(3,1)/scale));
    set(beam_S_0lab,'Color','blue','Linewidth',2,'MaxHeadSize',0.5,'AutoScale','off');
    text(0,0,-norm(S.S_0lab)/scale,'S_{0lab}','Color','blue','FontSize',14);
    beam_S_lab = quiver3(0,0,0,0.9*(S.S_lab(1,1)/scale),0.9*(S.S_lab(2,1)/scale),0.9*(S.S_lab(3,1)/scale));
    set(beam_S_lab,'Color','red','Linewidth',2,'MaxHeadSize',0.5,'AutoScale','off');
    text((S.S_lab(1,1)/scale),(S.S_lab(2,1)/scale),(S.S_lab(3,1)/scale),'S_{lab}','Color','red','FontSize',14);
    beam_Q_lab = quiver3(0,0,0,0.9*(S.Q_lab(1,1)/scale),0.9*(S.Q_lab(2,1)/scale),0.9*(S.Q_lab(3,1)/scale));
    set(beam_Q_lab,'Color','green','Linewidth',2,'MaxHeadSize',0.5,'AutoScale','off');
    text((S.Q_lab(1,1)/scale),(S.Q_lab(2,1)/scale),(S.Q_lab(3,1)/scale),'Q_{lab}','Color','green','FontSize',14);
end
    
% x, y and z axes
if Options.axes == 1
    fprintf('\n...plotting axes...');
    beam_x = quiver3(0,0,0,(0.9*norm(S.S_0lab)/scale),(0),(0));
    set(beam_x,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((1*norm(S.S_0lab)/scale),(0),(0),'x','Color','black','FontSize',14);
    beam_y = quiver3(0,0,0,(0),(0.9*norm(S.S_0lab)/scale),(0));
    set(beam_y,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((0),(1*norm(S.S_0lab)/scale),(0),'y','Color','black','FontSize',14);
    beam_z = quiver3(0,0,0,(0),(0),(0.9*norm(S.S_0lab)/scale));
    set(beam_z,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((0),(0),(1*norm(S.S_0lab)/scale),'z','Color','black','FontSize',14);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting lab space shape from detector conjugated space
if isfield(S,'LS_shape_DCS') && Options.plot_SIM_rev == 1    
    fprintf('\n...plotting LS shape from DCS...');
    % shifting LS shape from DCS to the centre of mass
    fprintf('\n...shifting LS shape from DCS to the centre of mass...');
    LS_shape_DCS = S.LS_shape_DCS;
    LS_shape_DCS_MASK = single(abs(LS_shape_DCS) > Options.mask_threshold);
    structure_element = strel('sphere', 3);
    LS_shape_DCS_MASK = imerode(imdilate(LS_shape_DCS_MASK, structure_element),structure_element); % takes care of dislocation cores
    LS_shape_DCS_COM = ceil(centerOfMass(LS_shape_DCS_MASK));
    LS_shape_DCS = circshift(LS_shape_DCS, size(LS_shape_DCS)/2-LS_shape_DCS_COM);

    % plotting LS shape from DCS
    plot_LS_shape_DCS = patch(isosurface(N1grid, N2grid, N3grid, LS_shape_DCS));
    set(plot_LS_shape_DCS, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    if Options.plot_SIM_fwd == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS], 'LS shape from SS', 'LS shape from DCS');
    else
    	legend('LS shape from DCS');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting lab space shape from reconstruction
% locating the reconstruction file
if isfield(Options, 'dir') && isfield(Options, 'folder_name')
    Options.temp_dir_LAB = [Options.dir, '/', Options.folder_name, '/', Options.folder_name,'-LAB.mat'];
elseif isfield(Options, 'dir')
    Options.folder_name = '';
    Options.temp_dir_LAB = [Options.dir, '/', Options.folder_name, '/', Options.folder_name,'-LAB.mat'];
elseif isfield(Options, 'folder_name')
    Options.temp_dir_LAB = [Options.folder_name, '/', Options.folder_name,'-LAB.mat'];
else
    Options.temp_dir_LAB = '';
end
    
if Options.plot_REC == 1 && exist(Options.temp_dir_LAB, 'file') == 2
    fprintf('\n...plotting LS shape from REC...');
    % loading LS shape from reconstruction
    LS_shape_REC = cell2mat(struct2cell(load(Options.temp_dir_LAB)));

    % taking conjugate reflection if necessary
    if Options.twin == 1
        F = ifftshift(fftn(fftshift(LS_shape_REC)));
        LS_shape_REC = fftshift(ifftn(ifftshift(conj(F))));
    end

    % shifting LS shape from REC to the centre of mass
    fprintf('\n...shifting LS shape from REC to the centre of mass...');
    LS_shape_REC_MASK = single(abs(LS_shape_REC) > Options.mask_threshold);
    structure_element = strel('sphere', 3);
    LS_shape_REC_MASK = imerode(imdilate(LS_shape_REC_MASK, structure_element),structure_element); % takes care of dislocation cores
    LS_shape_REC_COM = ceil(centerOfMass(LS_shape_REC_MASK));
    LS_shape_REC = circshift(LS_shape_REC, size(LS_shape_REC)/2-LS_shape_REC_COM);

    % plotting LS shape from reconstruction
    if isfield(S,'p_sam_REC') && S.p_sam_REC > 0
        fprintf('\n...interpolating simulation final pixel size and reconstruction final pixel size...')
        p_sam_REC = S.p_sam_REC; % Final Sample Pixel Size from reconstruction Command Window Output
        [RECX, RECY, RECZ] = meshgrid(-(size(LS_shape_REC,1)-1)/2*p_sam_REC:p_sam_REC:(size(LS_shape_REC,1)-1)/2*p_sam_REC, -(size(LS_shape_REC,2)-1)/2*p_sam_REC:p_sam_REC:(size(LS_shape_REC,2)-1)/2*p_sam_REC, -(size(LS_shape_REC,3)-1)/2*p_sam_REC:p_sam_REC:(size(LS_shape_REC,3)-1)/2*p_sam_REC);
        LS_shape_REC = interp3(RECX, RECY, RECZ, LS_shape_REC, N1grid, N2grid, N3grid);
    end
    LS_shape_REC_AMP  = abs(LS_shape_REC);
    LS_shape_REC_PH = angle(LS_shape_REC);
    [faces,verts,colors] = isosurface(N1grid, N2grid, N3grid, LS_shape_REC_AMP, Options.plot_threshold, LS_shape_REC_PH);
    plot_LS_shape_REC = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
    c = colorbar;
    ylabel(c, 'Phase');

    % putting legend in figure
    if isfield(S,'LS_shape_DCS') && Options.plot_SIM_fwd == 1 && Options.plot_SIM_rev == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS, plot_LS_shape_REC], 'LS shape from SS', 'LS shape from DCS', 'LS shape from REC')
    elseif isfield(S,'LS_shape_DCS') && Options.plot_SIM_rev == 1
        legend([plot_LS_shape_DCS, plot_LS_shape_REC], 'LS shape from DCS', 'LS shape from REC')
    elseif Options.plot_SIM_fwd == 1 && Options.plot_SIM_rev == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS], 'LS shape from SS', 'LS shape from DCS');
    elseif Options.plot_SIM_fwd == 1
        legend([plot_LS_shape_SS, plot_LS_shape_REC], 'LS shape from SS', 'LS shape from REC');
    else
        legend(plot_LS_shape_REC, 'LS shape from REC')
    end

    % overlap textbox
    fprintf('\n...calculating overlap...');
    % centring masks for overlap calculation
    LS_shape_SS_MASK = circshift(LS_shape_SS_MASK, size(LS_shape_SS_MASK)/2-LS_shape_SS_COM);
    LS_shape_REC_MASK = circshift(LS_shape_REC_MASK, size(LS_shape_REC_MASK)/2-LS_shape_REC_COM);
    % calculating overlap
    LS_shape_SS_REC_overlap = round(abs((1-abs(sum(sum(sum(LS_shape_REC_MASK - LS_shape_SS_MASK))))/sum(sum(sum(LS_shape_SS_MASK))))*100), 2);
    annotation('textbox',[0.17, 0.1, .3, .3], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','String',['Overlap: ',num2str(LS_shape_SS_REC_overlap),'%'], 'BackgroundColor', 'white','FitBoxToText','on');
else
    fprintf('\n...cannot find the reconstruction folder and/or file...');
    if Options.plot_SIM_fwd == 1 && Options.plot_SIM_rev == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS], 'LS shape from SS', 'LS shape from DCS');
    elseif Options.plot_SIM_fwd == 1
        legend(plot_LS_shape_SS, 'LS shape from SS');
    end
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
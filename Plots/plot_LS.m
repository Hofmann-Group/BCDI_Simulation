function plot_LS(S, axes_plot, beams_plot, grid_plot, viewpoint, temp_dir_LAB, mask_threshold, plot_threshold, twin, LS_shape_SS_plot, LS_shape_DCS_plot, LS_shape_REC_plot)
fprintf('plot_LS');
%% Plotting lab space shape from sample space 
fprintf('\n...plotting LS shape from SS...');
figure('Name','Lab Space (LS)');
hold on;

% shifting LS shape from SS to the centre of mass
fprintf('\n...shifting LS shape from SS to the centre of mass...');
LS_shape_SS = S.LS_shape_SS;
LS_shape_SS_MASK = single(abs(LS_shape_SS) > mask_threshold);
structure_element = strel('sphere', 3);
LS_shape_SS_MASK = imerode(imdilate(LS_shape_SS_MASK, structure_element),structure_element); % takes care of dislocation cores
LS_shape_SS_COM = ceil(centerOfMass(LS_shape_SS_MASK));
LS_shape_SS = circshift(LS_shape_SS, size(LS_shape_SS)/2-LS_shape_SS_COM);

% plotting LS shape from SS
if LS_shape_SS_plot == 1
    plot_LS_shape_SS = patch(isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, LS_shape_SS));
    set(plot_LS_shape_SS, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 1);
    legend([plot_LS_shape_SS], 'LS shape from SS');
end

% plot parameters
title(['Lab Space Shape with detector at \gamma = ', num2str(-S.gamma_spec), char(176), ' and \delta = ', num2str(S.delta_spec), char(176)]); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); 
set(gca,'XDir','normal');
set(gca,'YDir','normal');
daspect([1,1,1]);
axis equal;
axis vis3d xy;
view(viewpoint(1), viewpoint(2));
lighting gouraud;
camlight('headlight');
xlim([min(min(min(S.N1grid*S.p_sam))) max(max(max(S.N1grid*S.p_sam)))]); ylim([min(min(min(S.N2grid*S.p_sam))) max(max(max(S.N2grid*S.p_sam)))]); zlim([min(min(min(S.N3grid*S.p_sam))) max(max(max(S.N3grid*S.p_sam)))]);
if grid_plot == 1
    grid on;
end

% diffraction beams
if beams_plot == 1
    fprintf('\n...plotting diffraction beams...');
    scale = norm(S.S_lab)/S.N*2/(S.p_sam);
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
if axes_plot == 1
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
if isfield(S,'LS_shape_DCS') && LS_shape_DCS_plot == 1    
    fprintf('\n...plotting LS shape from DCS...');
    % shifting LS shape from DCS to the centre of mass
    fprintf('\n...shifting LS shape from DCS to the centre of mass...');
    LS_shape_DCS = S.LS_shape_DCS;
    LS_shape_DCS_MASK = single(abs(LS_shape_DCS) > mask_threshold);
    structure_element = strel('sphere', 3);
    LS_shape_DCS_MASK = imerode(imdilate(LS_shape_DCS_MASK, structure_element),structure_element); % takes care of dislocation cores
    LS_shape_DCS_COM = ceil(centerOfMass(LS_shape_DCS_MASK));
    LS_shape_DCS = circshift(LS_shape_DCS, size(LS_shape_DCS)/2-LS_shape_DCS_COM);

    % plotting LS shape from DCS
    plot_LS_shape_DCS = patch(isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, LS_shape_DCS));
    set(plot_LS_shape_DCS, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    if LS_shape_SS_plot == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS], 'LS shape from SS', 'LS shape from DCS');
    else
    	legend('LS shape from DCS');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting lab space shape from reconstruction
if LS_shape_REC_plot == 1 && exist(sprintf('%s', temp_dir_LAB), 'file') == 2
    fprintf('\n...plotting LS shape from REC...');
    % loading LS shape from reconstruction
    LS_shape_REC = cell2mat(struct2cell(load(temp_dir_LAB)));
%         LS_shape_REC = padcroparray(LS shape from DCS, S.nnc);

    % taking conjugate reflection if necessary
    if twin == 1
        F = ifftshift(fftn(fftshift(LS_shape_REC)));
        LS_shape_REC = fftshift(ifftn(ifftshift(conj(F))));
    end

    % shifting LS shape from REC to the centre of mass
    fprintf('\n...shifting LS shape from REC to the centre of mass...');
    LS_shape_REC_MASK = single(abs(LS_shape_REC) > mask_threshold);
    structure_element = strel('sphere', 3);
    LS_shape_REC_MASK = imerode(imdilate(LS_shape_REC_MASK, structure_element),structure_element); % takes care of dislocation cores
    LS_shape_REC_COM = ceil(centerOfMass(LS_shape_REC_MASK));
    LS_shape_REC = circshift(LS_shape_REC, size(LS_shape_REC)/2-LS_shape_REC_COM);

    % plotting LS shape from reconstruction
    if isfield(S,'p_sam_REC') && S.p_sam_REC > 0
        fprintf('\n...interpolating simulation final pixel size and reconstruction final pixel size...')
        p_sam_REC = S.p_sam_REC; % Final Sample Pixel Size from reconstruction Command Window Output
        [RECX, RECY, RECZ] = meshgrid(-(size(LS_shape_REC,1)-1)/2*p_sam_REC:p_sam_REC:(size(LS_shape_REC,1)-1)/2*p_sam_REC, -(size(LS_shape_REC,2)-1)/2*p_sam_REC:p_sam_REC:(size(LS_shape_REC,2)-1)/2*p_sam_REC, -(size(LS_shape_REC,3)-1)/2*p_sam_REC:p_sam_REC:(size(LS_shape_REC,3)-1)/2*p_sam_REC);
        LS_shape_REC = interp3(RECX, RECY, RECZ, LS_shape_REC, S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam);
    end
    LS_shape_REC_AMP  = abs(LS_shape_REC);
    LS_shape_REC_PH = angle(LS_shape_REC);
%         plot_LS_shape_REC = patch(isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, LS_shape_REC_AMP));
%         set(plot_LS_shape_REC, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    [faces,verts,colors] = isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, LS_shape_REC_AMP, plot_threshold, LS_shape_REC_PH);
    plot_LS_shape_REC = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
    c = colorbar;
    ylabel(c, 'Phase');

    % putting legend in figure
    if isfield(S,'LS_shape_DCS') && LS_shape_SS_plot == 1 && LS_shape_DCS_plot == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS, plot_LS_shape_REC], 'LS shape from SS', 'LS shape from DCS', 'LS shape from REC')
    elseif isfield(S,'LS_shape_DCS') && LS_shape_DCS_plot == 1
        legend([plot_LS_shape_DCS, plot_LS_shape_REC], 'LS shape from DCS', 'LS shape from REC')
    elseif LS_shape_SS_plot == 1 && LS_shape_DCS_plot == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS], 'LS shape from SS', 'LS shape from DCS');
    elseif LS_shape_SS_plot == 1
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
    if LS_shape_SS_plot == 1 && LS_shape_DCS_plot == 1
        legend([plot_LS_shape_SS, plot_LS_shape_DCS], 'LS shape from SS', 'LS shape from DCS');
    elseif LS_shape_SS_plot == 1
        legend(plot_LS_shape_SS, 'LS shape from SS');
    end
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function plot_SS(S, axes_plot, beams_plot, grid_plot, viewpoint, temp_dir_LAB, mask_threshold, plot_threshold, twin, SS_shape_plot, SS_shape_LS_plot, SS_shape_REC_plot)
fprintf('plot_SS');
%% Plotting sample space shape from user
fprintf('\n...plotting SS shape...');
figure('Name','Sample Space (SS)');
hold on;

% shifting original object to the centre of mass
fprintf('\n...shifting SS shape to the centre of mass...');
SS_shape = S.SS_shape;
SS_shape_MASK = single(abs(SS_shape) > mask_threshold);
structure_element = strel('sphere', 3);
SS_shape_MASK = imerode(imdilate(SS_shape_MASK, structure_element),structure_element); % takes care of dislocation cores
SS_shape_COM = ceil(centerOfMass(SS_shape_MASK));
SS_shape = circshift(SS_shape, size(SS_shape)/2-SS_shape_COM);

% original SS shape
if SS_shape_plot == 1
    plot_SS_shape = patch(isosurface(S.realXgrid, S.realYgrid, S.realZgrid, SS_shape));
    set(plot_SS_shape, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 1);
end

% plot parameters
title(['Sample Space Shape with detector at \gamma = ', num2str(-S.gamma_spec), char(176), ' and \delta = ', num2str(S.delta_spec), char(176)]); xlabel('x_{sam} (m)'); ylabel('y_{sam} (m)'); zlabel('z_{sam} (m)'); 
set(gca,'XDir','normal');
set(gca,'YDir','normal');
daspect([1,1,1]);
axis equal;
axis vis3d xy;
view(viewpoint(1), viewpoint(2));
lighting gouraud;
camlight('left');
xlim([min(min(min(S.realXgrid))) max(max(max(S.realXgrid)))]); ylim([min(min(min(S.realYgrid))) max(max(max(S.realYgrid)))]); zlim([min(min(min(S.realZgrid))) max(max(max(S.realZgrid)))]);
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
    
% x_sam, y_sam and z_sam axes
if axes_plot == 1
    fprintf('\n...plotting axes...');
    x_sam_axis = quiver3(0,0,0,(0.9*norm(S.S_0lab)/scale),(0),(0));
    set(x_sam_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((1*norm(S.S_0lab)/scale),(0),(0),'x_{sam}','Color','black','FontSize',14);
    y_sam_axis = quiver3(0,0,0,(0),(0.9*norm(S.S_0lab)/scale),(0));
    set(y_sam_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((0),(1*norm(S.S_0lab)/scale),(0),'y_{sam}','Color','black','FontSize',14);
    z_sam_axis = quiver3(0,0,0,(0),(0),(0.9*norm(S.S_0lab)/scale));
    set(z_sam_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((0),(0),(1*norm(S.S_0lab)/scale),'z_{sam}','Color','black','FontSize',14);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting sample space shape from lab space
if SS_shape_LS_plot == 1   
    fprintf('\n...plotting SS shape from LS...');
    % shifting SS shape from LS to the centre of mass
    fprintf('\n...shifting SS shape from LS to the centre of mass...');
    SS_shape_LS = LS_to_SS(S, S.LS_shape_DCS);
    SS_shape_LS_MASK = single(abs(SS_shape_LS) > mask_threshold);
    structure_element = strel('sphere', 3);
    SS_shape_LS_MASK = imerode(imdilate(SS_shape_LS_MASK, structure_element),structure_element); % takes care of dislocation cores
    SS_shape_LS_COM = ceil(centerOfMass(SS_shape_LS_MASK));
    SS_shape_LS = circshift(SS_shape_LS, size(SS_shape_LS)/2-SS_shape_LS_COM);

    % plotting LS shape from DCS
    plot_SS_shape_LS = patch(isosurface(S.realXgrid, S.realYgrid, S.realZgrid, SS_shape_LS));
    set(plot_SS_shape_LS, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    
    if SS_shape_plot == 1
        legend([plot_SS_shape, plot_SS_shape_LS], 'SS shape', 'SS shape from LS');
    else
    	legend('SS shape from LS');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting sample space shape from reconstruction
if SS_shape_REC_plot == 1 && exist(sprintf('%s', temp_dir_LAB), 'file') == 2
    fprintf('\n...plotting SS shape from REC...');
    % loading SS shape from reconstruction and mapping it to SS
    LS_shape_REC = cell2mat(struct2cell(load(temp_dir_LAB)));    
    SS_shape_REC = LS_to_SS(S, LS_shape_REC);
    
    % taking conjugate reflection if necessary
    if twin == 1
        F = ifftshift(fftn(fftshift(SS_shape_REC)));
        SS_shape_REC = fftshift(ifftn(ifftshift(conj(F))));
    end

    % shifting SS shape from REC to the centre of mass
    fprintf('\n...shifting SS shape from REC to the centre of mass...');
    SS_shape_REC_MASK = single(abs(SS_shape_REC) > mask_threshold);
    structure_element = strel('sphere', 3);
    SS_shape_REC_MASK = imerode(imdilate(SS_shape_REC_MASK, structure_element),structure_element); % takes care of dislocation cores
    SS_shape_REC_COM = ceil(centerOfMass(SS_shape_REC_MASK));
    SS_shape_REC = circshift(SS_shape_REC, size(SS_shape_REC)/2-SS_shape_REC_COM);

    % plotting SS shape from reconstruction
    if isfield(S,'p_sam_REC') && S.p_sam_REC > 0
        fprintf('\n...interpolating simulation final pixel size and reconstruction final pixel size...')
        p_sam_REC = S.p_sam_REC; % Final Sample Pixel Size from reconstruction Command Window Output
        [RECX, RECY, RECZ] = meshgrid(-(size(SS_shape_REC,1)-1)/2*p_sam_REC:p_sam_REC:(size(SS_shape_REC,1)-1)/2*p_sam_REC, -(size(SS_shape_REC,2)-1)/2*p_sam_REC:p_sam_REC:(size(SS_shape_REC,2)-1)/2*p_sam_REC, -(size(SS_shape_REC,3)-1)/2*p_sam_REC:p_sam_REC:(size(SS_shape_REC,3)-1)/2*p_sam_REC);
        SS_shape_REC = interp3(RECX, RECY, RECZ, SS_shape_REC, S.realXgrid, S.realYgrid, S.realZgrid);
    end
    SS_shape_REC_AMP  = abs(SS_shape_REC);
    SS_shape_REC_PH = angle(SS_shape_REC);
%         plot_SS_shape_REC = patch(isosurface(S.realXgrid, S.realYgrid, S.realZgrid, SS_shape_REC_AMP));
%         set(plot_SS_shape_REC, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    [faces,verts,colors] = isosurface(S.realXgrid, S.realYgrid, S.realZgrid, SS_shape_REC_AMP, plot_threshold, SS_shape_REC_PH);
    plot_SS_shape_REC = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
    c = colorbar;
    ylabel(c, 'Phase');

    % putting legend in figure
    if SS_shape_plot == 1 && SS_shape_LS_plot == 1
        legend([plot_SS_shape, plot_SS_shape_LS, plot_SS_shape_REC], 'SS shape from SS', 'SS shape from LS', 'SS shape from REC')
    elseif SS_shape_LS_plot == 1
        legend([plot_SS_shape_LS, plot_SS_shape_REC], 'SS shape from LS', 'SS shape from REC')
    elseif SS_shape_plot == 1 && SS_shape_LS_plot == 1
        legend([plot_SS_shape, plot_SS_shape_LS], 'SS shape from SS', 'SS shape from LS');
    elseif SS_shape_plot == 1
        legend([plot_SS_shape, plot_SS_shape_REC], 'SS shape from SS', 'SS shape from REC');
    else
        legend(plot_SS_shape_REC, 'SS shape from REC')
    end

    % overlap textbox
    fprintf('\n...calculating overlap...');
    % centring masks for overlap calculation
    SS_shape_MASK = circshift(SS_shape_MASK, size(SS_shape_MASK)/2-SS_shape_COM);
    SS_shape_REC_MASK = circshift(SS_shape_REC_MASK, size(SS_shape_REC_MASK)/2-SS_shape_REC_COM);
    % calculating overlap
    SS_shape_REC_overlap = round(abs((1-abs(sum(sum(sum(SS_shape_REC_MASK - SS_shape_MASK))))/sum(sum(sum(SS_shape_MASK))))*100), 2);
    annotation('textbox',[0.17, 0.1, .3, .3], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','String',['Overlap: ',num2str(SS_shape_REC_overlap),'%'], 'BackgroundColor', 'white','FitBoxToText','on');
else
    fprintf('\n...cannot find the reconstruction folder and/or file...');
    if SS_shape_plot == 1 && SS_shape_LS_plot == 1
        legend([plot_SS_shape, plot_SS_shape_LS], 'SS shape', 'SS shape from LS');
    elseif SS_shape_plot == 1
        legend(plot_SS_shape, 'SS shape');
    end
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
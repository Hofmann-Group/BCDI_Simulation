function plot_DCS(S, axes_plot, beams_plot, grid_plot, viewpoint, temp_dir_AMP, temp_dir_PH, mask_threshold, plot_threshold, twin, DCS_shape_DRS_plot, DCS_shape_LS_plot, DCS_shape_REC_plot)
fprintf('plot_DCS');
%% Plotting detector conjugated space shape from detector reciprocal space
fprintf('\n...plotting DCS shape from DRS...');
figure('Name','Detector Conjugated Space (DCS)');
hold on;

% shifting DCS shape from DRS to the centre of mass
fprintf('\n...shifting DCS shape from DRS to the centre of mass...');
DCS_shape_DRS = S.DCS_shape_DRS;
DCS_shape_DRS_MASK = single(abs(DCS_shape_DRS) > 0.3);
structure_element = strel('sphere', 3);
DCS_shape_DRS_MASK = imerode(imdilate(DCS_shape_DRS_MASK, structure_element),structure_element); % takes care of dislocation cores
DCS_shape_DRS_COM = ceil(centerOfMass(DCS_shape_DRS_MASK));
DCS_shape_DRS = circshift(DCS_shape_DRS, size(DCS_shape_DRS)/2-DCS_shape_DRS_COM);

% plotting DCS shape from DRS
if DCS_shape_DRS_plot == 1
    plot_DCS_shape_DRS = patch(isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, DCS_shape_DRS));
    set(plot_DCS_shape_DRS, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    legend('DCS shape from DRS')
end

% plot parameters
title(['Detector Conjugated Space Shape with detector at \gamma = ', num2str(-S.gamma_spec), char(176), ' and \delta = ', num2str(S.delta_spec), char(176)]); xlabel('x'' (m)'); ylabel('y'' (m)'); zlabel('z'' (m)');
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
    x_1p_axis = quiver3(0,0,0,(0.9*norm(S.S_0lab)/scale),(0),(0));
    set(x_1p_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((1*norm(S.S_0lab)/scale),(0),(0),'x''','Color','black','FontSize',14);
    y_1p_axis = quiver3(0,0,0,(0),(0.9*norm(S.S_0lab)/scale),(0));
    set(y_1p_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((0),(1*norm(S.S_0lab)/scale),(0),'y''','Color','black','FontSize',14);
    z_1p_axis = quiver3(0,0,0,(0),(0),(0.9*norm(S.S_0lab)/scale));
    set(z_1p_axis,'Color','black','Linewidth', 2, 'AutoScale','off');
    text((0),(0),(1*norm(S.S_0lab)/scale),'z''','Color','black','FontSize',14);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting detector conjugated space shape from lab space
if isfield(S,'DCS_shape_LS') && DCS_shape_LS_plot == 1    
    fprintf('\n...plotting DCS shape from LS...');
    % shifting DCS shape from LS to the centre of mass
    fprintf('\n...shifting DCS shape from LS to the centre of mass...');
    DCS_shape_LS = S.DCS_shape_LS;
    DCS_shape_LS_MASK = single(abs(DCS_shape_LS) > mask_threshold);
    structure_element = strel('sphere', 3);
    DCS_shape_LS_MASK = imerode(imdilate(DCS_shape_LS_MASK, structure_element),structure_element); % takes care of dislocation cores
    DCS_shape_LS_COM = ceil(centerOfMass(DCS_shape_LS_MASK));
    DCS_shape_LS = circshift(DCS_shape_LS, size(DCS_shape_LS)/2-DCS_shape_LS_COM);

    % detector conjugated space shape from lab space
    plot_DCS_shape_LS = patch(isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, DCS_shape_LS));
    set(plot_DCS_shape_LS, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha',0.3);
    
    if DCS_shape_DRS_plot == 1
        legend([plot_DCS_shape_DRS, plot_DCS_shape_LS], 'DCS shape from DRS', 'DCS shape from LS')
    else
    	legend('DCS shape from LS')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting detector conjugated shape from reconstruction
try
    if DCS_shape_REC_plot == 1 && exist(sprintf('%s', temp_dir_AMP), 'file') == 2 && exist(sprintf('%s', temp_dir_PH), 'file') == 2
        fprintf('\n...plotting DCS shape from REC...');
        % loading reconstructed detector conjugated shape
        DCS_shape_REC_AMP = cell2mat(struct2cell(load(temp_dir_AMP, '-mat')));
        DCS_shape_REC_PH = cell2mat(struct2cell(load(temp_dir_PH, '-mat'))); 
        DCS_shape_REC = DCS_shape_REC_AMP.*exp(1i.*DCS_shape_REC_PH);  

        % taking conjugate reflection if necessary
        if twin ~= 1 %~=1 instead of == 1 to match it with the detector conjugated shape
            F = ifftshift(fftn(fftshift(DCS_shape_REC)));
            DCS_shape_REC = fftshift(ifftn(ifftshift(conj(F))));
        end

        % shifting DCS shape from REC to the centre of mass
        fprintf('\n...shifting DCS shape from REC to the centre of mass...');
        DCS_shape_REC_MASK = single(abs(DCS_shape_REC) > mask_threshold);
        structure_element = strel('sphere', 3);
        DCS_shape_REC_MASK = imerode(imdilate(DCS_shape_REC_MASK, structure_element),structure_element); % takes care of dislocation cores
        DCS_shape_REC_COM = ceil(centerOfMass(DCS_shape_REC_MASK));
        DCS_shape_REC = circshift(DCS_shape_REC, size(DCS_shape_REC)/2-DCS_shape_REC_COM);

        % plotting reconstructed detector conjugated shape
        if isfield(S,'p_sam_REC') && S.p_sam_REC > 0
            fprintf('\n...interpolating simulation final pixel size and reconstruction final pixel size...')
            p_sam_REC = S.p_sam_REC; % Final Sample Pixel Size from reconstruction Command Window Output
            [RECX, RECY, RECZ] = meshgrid(-(size(DCS_shape_REC,1)-1)/2*p_sam_REC:p_sam_REC:(size(DCS_shape_REC,1)-1)/2*p_sam_REC, -(size(DCS_shape_REC,2)-1)/2*p_sam_REC:p_sam_REC:(size(DCS_shape_REC,2)-1)/2*p_sam_REC, -(size(DCS_shape_REC,3)-1)/2*p_sam_REC:p_sam_REC:(size(DCS_shape_REC,3)-1)/2*p_sam_REC);
            DCS_shape_REC = interp3(RECX, RECY, RECZ, DCS_shape_REC, S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam);
        end
        DCS_shape_REC_AMP  = abs(DCS_shape_REC);
        DCS_shape_REC_PH = angle(DCS_shape_REC);
        [faces,verts,colors] = isosurface(S.N1grid*S.p_sam, S.N2grid*S.p_sam, S.N3grid*S.p_sam, DCS_shape_REC_AMP, plot_threshold, DCS_shape_REC_PH);
        plot_DCS_shape_REC = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
        c = colorbar;
        ylabel(c, 'Phase');

        % putting legend in figure
        if isfield(S,'DCS_shape_LS')  && DCS_shape_DRS_plot == 1 && DCS_shape_LS_plot == 1
            legend([plot_DCS_shape_DRS, plot_DCS_shape_LS, plot_DCS_shape_REC], 'DCS shape from DRS', 'DCS shape from LS', 'DCS shape from REC')
        elseif isfield(S,'DCS_shape_LS') && DCS_shape_LS_plot == 1
            legend([plot_DCS_shape_LS, plot_DCS_shape_REC], 'DCS shape from LS', 'DCS shape from REC')
        elseif DCS_shape_DRS_plot == 1 && DCS_shape_LS_plot == 1
            legend([plot_DCS_shape_DRS, plot_DCS_shape_LS], 'DCS shape from DRS', 'DCS shape from LS')
        elseif DCS_shape_DRS_plot == 1
            legend([plot_DCS_shape_DRS, plot_DCS_shape_REC], 'DCS shape from DRS', 'DCS shape from REC');
        else
            legend([plot_DCS_shape_REC], 'DCS shape from REC')
        end

        % overlap textbox
        fprintf('\n...calculating overlap...');
        % centring masks for overlap calculation
        DCS_shape_DRS_MASK = circshift(DCS_shape_DRS_MASK, size(DCS_shape_DRS_MASK)/2-DCS_shape_DRS_COM);
        DCS_shape_REC_MASK = circshift(DCS_shape_REC_MASK, size(DCS_shape_REC_MASK)/2-DCS_shape_REC_COM);
        % calculating overlap
        DCS_shape_DRS_REC_overlap = round(abs((1-abs(sum(sum(sum(DCS_shape_REC_MASK - DCS_shape_DRS_MASK))))/sum(sum(sum(DCS_shape_DRS_MASK))))*100), 2);
        annotation('textbox',[0.17, 0.1, .3, .3], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','String',['Overlap: ',num2str(DCS_shape_DRS_REC_overlap),'%'], 'BackgroundColor', 'white','FitBoxToText','on');
    end
catch
    fprintf('\n...cannot find the reconstruction folder and/or file...');
    if DCS_shape_DRS_plot == 1 && DCS_shape_LS_plot == 1
        legend([plot_DCS_shape_DRS, plot_DCS_shape_LS], 'DCS shape from DRS', 'DCS shape from LS')
    elseif DCS_shape_DRS_plot == 1
        legend(plot_DCS_shape_DRS, 'DCS shape from DRS')
    end
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
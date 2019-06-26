function plot_DCS(S, Options)
% plots the detector conjugated space
fprintf('plot_DCS');
%% Plotting detector conjugated space shape from detector reciprocal space
fprintf('\n...plotting DCS shape from DRS...');
figure('Name','Detector Conjugated Space (DCS)');
hold on;

% shifting DCS shape from DRS to the centre of mass
fprintf('\n...shifting DCS shape from DRS to the centre of mass...');
DCS_shape_DRS_MASK = single(abs(S.DCS_shape_DRS) > 0.3);
structure_element = strel('sphere', 3);
DCS_shape_DRS_MASK = imerode(imdilate(DCS_shape_DRS_MASK, structure_element),structure_element); % takes care of dislocation cores
DCS_shape_DRS_COM = ceil(centerOfMass(DCS_shape_DRS_MASK));
S.DCS_shape_DRS = circshift(S.DCS_shape_DRS, size(S.DCS_shape_DRS)/2-DCS_shape_DRS_COM);

% plot grids
N1grid = S.N1grid*S.p_sam;
N2grid = S.N2grid*S.p_sam;
N3grid = S.N3grid*S.p_sam;

% plotting DCS shape from DRS
if Options.plot_SIM_fwd == 1
    plot_DCS_shape_DRS = patch(isosurface(N1grid, N2grid, N3grid, S.DCS_shape_DRS));
    set(plot_DCS_shape_DRS, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    legend('DCS shape from DRS')
end

% plot parameters
title(['Detector Conjugated Space Shape with detector at \gamma = ', num2str(-S.gamma_bl), char(176), ' and \delta = ', num2str(S.delta_bl), char(176)]); xlabel('x'' (m)'); ylabel('y'' (m)'); zlabel('z'' (m)');
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

% x, y and z axes
if Options.axes == 1
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
if isfield(S,'DCS_shape_LS') && Options.plot_SIM_rev == 1    
    fprintf('\n...plotting DCS shape from LS...');
    % shifting DCS shape from LS to the centre of mass
    fprintf('\n...shifting DCS shape from LS to the centre of mass...');
    DCS_shape_LS = S.DCS_shape_LS;
    DCS_shape_LS_MASK = single(abs(DCS_shape_LS) > Options.mask_threshold);
    structure_element = strel('sphere', 3);
    DCS_shape_LS_MASK = imerode(imdilate(DCS_shape_LS_MASK, structure_element),structure_element); % takes care of dislocation cores
    DCS_shape_LS_COM = ceil(centerOfMass(DCS_shape_LS_MASK));
    DCS_shape_LS = circshift(DCS_shape_LS, size(DCS_shape_LS)/2-DCS_shape_LS_COM);

    % detector conjugated space shape from lab space
    plot_DCS_shape_LS = patch(isosurface(N1grid, N2grid, N3grid, DCS_shape_LS));
    set(plot_DCS_shape_LS, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha',0.3);
    
    if Options.plot_SIM_fwd == 1
        legend([plot_DCS_shape_DRS, plot_DCS_shape_LS], 'DCS shape from DRS', 'DCS shape from LS')
    else
    	legend('DCS shape from LS')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting detector conjugated shape from reconstruction
% locating the reconstruction file
if isfield(Options, 'dir') && isfield(Options, 'folder_name')
    Options.temp_dir_AMP = [Options.dir, '/', Options.folder_name, '/', Options.folder_name,'-AMP.mat'];
    Options.temp_dir_PH = [Options.dir, '/', Options.folder_name, '/', Options.folder_name,'-PH.mat'];
elseif isfield(Options, 'dir')
    Options.folder_name = '';
    Options.temp_dir_AMP = [Options.dir, '/', Options.folder_name, '/', Options.folder_name,'-AMP.mat'];
    Options.temp_dir_PH = [Options.dir, '/', Options.folder_name, '/', Options.folder_name,'-PH.mat'];    
elseif isfield(Options, 'folder_name')
    Options.temp_dir_AMP = [Options.folder_name, '/', Options.folder_name,'-AMP.mat'];
    Options.temp_dir_PH = [Options.folder_name, '/', Options.folder_name,'-PH.mat'];  
else
    Options.temp_dir_AMP = '';
    Options.temp_dir_PH = '';
end
    
if Options.plot_REC == 1 && exist(Options.temp_dir_AMP, 'file') == 2 && exist(Options.temp_dir_PH, 'file') == 2
    fprintf('\n...plotting DCS shape from REC...');   
    % loading reconstructed detector conjugated shape
    DCS_shape_REC_AMP = cell2mat(struct2cell(load(Options.temp_dir_AMP, '-mat')));
    DCS_shape_REC_PH = cell2mat(struct2cell(load(Options.temp_dir_PH, '-mat'))); 
    DCS_shape_REC = DCS_shape_REC_AMP.*exp(1i.*DCS_shape_REC_PH);  

    % taking conjugate reflection if necessary
    if Options.twin ~= 1 %~=1 instead of == 1 to match it with the detector conjugated shape
        F = ifftshift(fftn(fftshift(DCS_shape_REC)));
        DCS_shape_REC = fftshift(ifftn(ifftshift(conj(F))));
    end

    % shifting DCS shape from REC to the centre of mass
    fprintf('\n...shifting DCS shape from REC to the centre of mass...');
    DCS_shape_REC_MASK = single(abs(DCS_shape_REC) > Options.mask_threshold);
    structure_element = strel('sphere', 3);
    DCS_shape_REC_MASK = imerode(imdilate(DCS_shape_REC_MASK, structure_element),structure_element); % takes care of dislocation cores
    DCS_shape_REC_COM = ceil(centerOfMass(DCS_shape_REC_MASK));
    DCS_shape_REC = circshift(DCS_shape_REC, size(DCS_shape_REC)/2-DCS_shape_REC_COM);

    % plotting reconstructed detector conjugated shape
    if isfield(S,'p_sam_REC') && S.p_sam_REC > 0
        fprintf('\n...interpolating simulation final pixel size and reconstruction final pixel size...')
        p_sam_REC = S.p_sam_REC; % Final Sample Pixel Size from reconstruction Command Window Output
        [RECX, RECY, RECZ] = meshgrid(-(size(DCS_shape_REC,1)-1)/2*p_sam_REC:p_sam_REC:(size(DCS_shape_REC,1)-1)/2*p_sam_REC, -(size(DCS_shape_REC,2)-1)/2*p_sam_REC:p_sam_REC:(size(DCS_shape_REC,2)-1)/2*p_sam_REC, -(size(DCS_shape_REC,3)-1)/2*p_sam_REC:p_sam_REC:(size(DCS_shape_REC,3)-1)/2*p_sam_REC);
        DCS_shape_REC = interp3(RECX, RECY, RECZ, DCS_shape_REC, N1grid, N2grid, N3grid);
    end
    DCS_shape_REC_AMP  = abs(DCS_shape_REC);
    DCS_shape_REC_PH = angle(DCS_shape_REC);
    [faces,verts,colors] = isosurface(N1grid, N2grid, N3grid, DCS_shape_REC_AMP, Options.plot_threshold, DCS_shape_REC_PH);
    plot_DCS_shape_REC = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
    c = colorbar;
    ylabel(c, 'Phase');

    % putting legend in figure
    if isfield(S,'DCS_shape_LS') && Options.plot_SIM_fwd == 1 && Options.plot_SIM_rev == 1
        legend([plot_DCS_shape_DRS, plot_DCS_shape_LS, plot_DCS_shape_REC], 'DCS shape from DRS', 'DCS shape from LS', 'DCS shape from REC')
    elseif isfield(S,'DCS_shape_LS') && Options.plot_SIM_rev == 1
        legend([plot_DCS_shape_LS, plot_DCS_shape_REC], 'DCS shape from LS', 'DCS shape from REC')
    elseif Options.plot_SIM_fwd == 1 && Options.plot_SIM_rev == 1
        legend([plot_DCS_shape_DRS, plot_DCS_shape_LS], 'DCS shape from DRS', 'DCS shape from LS')
    elseif Options.plot_SIM_fwd == 1
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

else
    fprintf('\n...cannot find the reconstruction folder and/or file...');
    if Options.plot_SIM_fwd == 1 && Options.plot_SIM_rev == 1
        legend([plot_DCS_shape_DRS, plot_DCS_shape_LS], 'DCS shape from DRS', 'DCS shape from LS')
    elseif Options.plot_SIM_fwd == 1
        legend(plot_DCS_shape_DRS, 'DCS shape from DRS')
    end
end
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% PostProcessingScript_3D_Geometry
% File path & prefix of output results
% Function: Plot 3D solid configuration of topology optimization results
% Description: Visualize original 3D geometric shape & material distribution
%              No deformation, displacement or buckling mode is superimposed

% ---------------- Plot 3D geometry for specified iteration steps ---------
% File path & prefix of output results
root = 'Mesh18x96x24_20260612091023';
prefix = 'Mesh18x96x24_BCV';
filetitle = [root '\' prefix];

% Plot settings
it = [300]; % 
threshold = 0.5;
isSave = 0;

% Load basic model parameters and mesh information
load([filetitle '_Input.mat']);
load([filetitle '_MeshInfo.mat']);

% Recover geometric dimensions in Y and Z directions
[Ly,Lz] = deal(nely/nelx*Lx,nelz/nelx*Lx);

% Maximum display displacement (reserved parameter)
maxDU = Lx/5;

% Loop over target iteration steps
for i = it
    % Load data of current iteration
    load([filetitle '_xPhysMat_It' num2str(i) '.mat']); % Physical density field 'xPhysMat' 
    load([filetitle '_U_It' num2str(i) '.mat']); % Global displacement vector 'U'
    load([filetitle '_phi_It' num2str(i) '.mat']); % Buckling mode shape vector 'phi'
    
    % Select valid solid elements above density threshold
    eIndex = xPhysMat(:)>threshold;
    faces = reshape(elements(eIndex,efaces')',4,[])';
    
     % Smooth density field, extract isosurfaces and end caps for 3D reconstruction
    isovals = smooth3(xPhysMat,'box',1);
    [isof1,isov1] = isosurface(centerX,centerY,centerZ,isovals,0.5);
    [isof2,isov2] = isocaps(centerX,centerY,centerZ,isovals,0.5);
    
    % Create figure window
    figName = sprintf('It_%d_Design', i);
    fig = figure('Name', figName);
    
    % Draw 3D solid geometry: main surface + end caps
    patch('Faces',isof1,'Vertices',isov1,'FaceColor','b','EdgeColor','none');
    patch('Faces',isof2,'Vertices',isov2,'FaceColor','#D95319','EdgeColor','none');
    
    % Light & Material Settings (Enhance 3D Sense) 
    lighting gouraud;                % Smooth lighting mode (better for curved surfaces)
    material dull;                   % Matte material, conform to engineering drawing style
    % Main light: Focus on XOZ plane (Y direction offset small, highlight XOZ surface)
    light('Position', [Lx, 0.1*Ly, Lz], 'Style', 'infinite');
    % Fill light: Soften shadows, prevent overly dark areas
    light('Position', [-0.5*Lx, 0.1*Ly, 0.8*Lz], 'Style', 'infinite', 'Color', [0.5 0.5 0.5]);
    % Increase ambient light to brighten overall scene
    set(gca, 'AmbientLightColor', [0.8 0.8 0.8]);
    
    
    % Figure layout and display settings
    title(['It.' num2str(i)]); 
    view([135,25]);
    axis equal;grid on;
    xlabel('X');ylabel('Y');zlabel('Z');
    
    % Save
    if isSave == 1
        % Define save path and filename
        saveName = fullfile(root, sprintf('%s_It%d_3D_Geometry.tif', prefix, i));
        
        % Check if file exists, pop up confirm dialog
        if exist(saveName, 'file')
            choice = questdlg('File already exists, overwrite it?', ...
                'File Overwrite Warning', 'Yes', 'No', 'No');
            if strcmp(choice, 'No')
                fprintf('Skip saving: %s\n', saveName);
                continue;
            end
        end
        
        % Export high-resolution image
        exportgraphics(gcf, saveName, 'Resolution', 600);
        fprintf('Image saved successfully: %s\n', saveName);
    end
end
%% PostProcessingScript_2D_Waveform_GIF
% Function: Generate GIF animation of 2D cross-section waveform
% Description: Create animated GIF for middle Z-section material distribution
%              Unified naming rule, overwrite prompt and standard code style

% ---------------- Plot 2D waveform animation for iteration steps ---------
% File path & prefix of output results
root = 'Mesh18x96x24_20260611173945';
prefix = 'Mesh18x96x24_BCV';
filetitle = [root '\' prefix];

% Animation & plot settings
totalLoop = 800;             % Total iteration steps for animation
delay = 0.02;                % GIF frame delay time (unit: second)
figWidth = 16;               % Figure width (unit: cm)
figHeight = 3;               % Figure height (unit: cm)
isSave = 0;                  % 1 = Enable GIF saving, 0 = Disable

% Define full save name with unified prefix
saveName = fullfile(root, sprintf('%s_It1-To-It%d_Waveform_Anim.gif', prefix, totalLoop));

% Load basic model parameters
load([filetitle '_Input.mat']);

% Create figure window and set layout
fig = figure();
set(gcf,'Units','centimeters','Position',[10,10,figWidth,figHeight]);
set(gca,'Units','centimeters','Position',[0,0,figWidth,figHeight]);

% Check file existence and pop up overwrite prompt
if isSave == 1 && exist(saveName, 'file')
    choice = questdlg('File already exists, overwrite it?', ...
        'File Overwrite Warning', 'Yes', 'No', 'No');
    if strcmp(choice, 'No')
        fprintf('Skip generating GIF: %s\n', saveName);
        return;
    end
end

% Loop over each iteration to generate animation frames
for idx = 1 : totalLoop
    % Load physical density field of current iteration
    load([filetitle '_xPhysMat_It' num2str(idx) '.mat']);
    
    % Update figure name
    fig.Name = sprintf('%s_Waveform_Frame_%d', prefix, idx);
    
    % Extract middle Z cross-section and plot grayscale waveform
    midSlice = floor(nelz / 2);
    imagesc(1 - xPhysMat(:,:,midSlice));
    colormap('gray');
    axis equal tight off;
    
    % Add iteration number text label
    text(nely/2, nelx-5, ['It. ' num2str(idx,'%3d')], ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    
    drawnow;  % Refresh figure immediately
    
    % Capture frame and write to GIF file
    if isSave == 1
        frame = getframe(fig);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
        if idx == 1
            % Create new GIF file for first frame
            imwrite(A, map, saveName, "gif", "LoopCount", Inf, "DelayTime", delay);
        else
            % Append subsequent frames to existing GIF
            imwrite(A, map, saveName, "gif", "WriteMode", "append", "DelayTime", delay);
        end
    end
end

if isSave == 1
    fprintf('GIF animation saved successfully: %s\n', saveName);
end
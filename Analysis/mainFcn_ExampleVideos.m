%% Information integration for decision-making in desert locusts
% Locust swarms can extend over several hundred kilometers, and starvation
% compels this ancient pest to devour everything in its path. Theory
% suggests that gregarious behavior benefits foraging efficiency, yet the
% role of social cohesion in locust foraging decisions remains elusive. To
% this end, we collected high-resolution tracking data of individual and
% grouped gregarious desert locusts in a 2-choice behavioral assay with
% animals deciding between patches of either similar or different quality.
% Carefully maintaining the animals' identities allowed us to monitor what
% each individual has experienced and to estimate the leaky accumulation
% process of personally acquired and, when available, socially derived
% evidence. We fitted these data to a model based on Bayesian estimation
% to gain insight into the locust social decision-making system for patch
% selection. By disentangling the relative contribution of each information
% class, our study suggests that locusts balance incongruent evidence but
% reinforce congruent ones. We provide insight into the collective foraging
% decisions of social (but non-eusocial) insects and present locusts as a
% powerful empirical system to study individual choices and their
% consequent collective dynamics.
%
% This is script creates an example video showing how a group of 15 animals
% were foraging under the unequal condition
%
% Version: 15-Jan-2022 (MATLAB R2022a)

% Tidy up
clc; clear all; close all
% Add path to external export function
% (https://github.com/altmany/export_fig)
addpath(genpath('altmany-export_fig'))
% Keep things down
warning('off')

%% Settings
SET.ExampleGrp = '15';
SET.ExampleCond = 'UE';
SET.ExampleTrialName = '15UE20191206';
SET.CurrDir = pwd;
SET.VideoName = [SET.CurrDir,'\ExampleVideo\',SET.ExampleTrialName,'.mp4'];
VidObj = VideoReader(SET.VideoName);
SET.TrackingName = [SET.CurrDir,'\ExampleVideo\',SET.ExampleTrialName,'_tracked.csv'];
Tracking = readtable(SET.TrackingName);
SET.AnnotationName = [SET.CurrDir,'\ExampleVideo\',SET.ExampleTrialName,'_annotation.mat'];
Annotation = load(SET.AnnotationName);
SET.InteractionRange = 7; %cm
SET.ArenaDiameter= 90; %cm
SET.FrameRate = 25;

% Max Density
bins = zeros(VidObj.Height, VidObj.Width);
bins(floor(VidObj.Height/2), floor(VidObj.Width/2)) = 1;
for iOrder = 1:2
    bins = imgaussfilt(bins, SET.InteractionRange/(SET.ArenaDiameter/(Annotation.Annotation.ROI.Par(3)*2)));
end%iOrder
SET.MaxDensity = max(max(bins));

%% Video
figure('Color', 'w', 'Units', 'normalized', 'Position', [0.125 0.125 0.5 0.5])
cnt_Frame = 1;
for iFrame = 50:50:45000

    % Get frame
    frame = read(VidObj, iFrame);
    % Processing
    frame = imadjust(rgb2gray(frame));
    frame = cat(3,frame,frame,frame);

    % Get tracking
    idx = find(Tracking.frame == iFrame);
    xPos = Tracking.pos_x(idx);
    yPos = Tracking.pos_y(idx);

    subplot(2,4,[1 2 5 6])
    hold on
    imagesc(frame)
    hCol = colorbar;
    hCol.Visible = 'off';
    % Depict Centroids
    mkr = scatter(Tracking.pos_x(idx), Tracking.pos_y(idx),...
        10,[189 042 132]/255,'filled'); mkr.MarkerEdgeColor = 'none'; clear mkr;
    % Indicate patches and arena
    th = 0:pi/50:2*pi;
    plot(Annotation.Annotation.ROI.Par(3) * cos(th) + Annotation.Annotation.ROI.Par(1),  Annotation.Annotation.ROI.Par(3) * sin(th) + Annotation.Annotation.ROI.Par(2), 'Color', [189 042 132]/255)
    plot(Annotation.Annotation.Masks.Circular(1,3) * cos(th) + Annotation.Annotation.Masks.Circular(1,1),  Annotation.Annotation.Masks.Circular(1,3) * sin(th) + Annotation.Annotation.Masks.Circular(1,2), 'Color', [189 042 132]/255)
    plot(Annotation.Annotation.Masks.Circular(2,3) * cos(th) + Annotation.Annotation.Masks.Circular(2,1),  Annotation.Annotation.Masks.Circular(2,3) * sin(th) + Annotation.Annotation.Masks.Circular(2,2), 'Color', [189 042 132]/255)
    % Iterate over all animals and indicate their IDs
    ID_maxDist = mean(size(rgb2gray(frame)))*0.05;
    for iAni = 1:length(idx)

        % Get center of frame
        ID_center = [size(frame,1)/2; size(frame,2)/2];
        % Get angle to ID center
        ID_dir = [Tracking.pos_x(idx(iAni)); Tracking.pos_y(idx(iAni))]- ID_center;
        if ID_dir(1)==0 && ID_dir(2)==0
            ID_dir = [1 1];
        end
        ID_len = ID_maxDist;
        % Get position of anootation
        ID_dir = (ID_dir/norm(ID_dir))*ID_len;
        ID_pos_text = ID_dir + 10 + [Tracking.pos_x(idx(iAni)); Tracking.pos_y(idx(iAni))];
        ID_pos_line = ID_dir*0.85 + [Tracking.pos_x(idx(iAni)); Tracking.pos_y(idx(iAni))];
        if atan2d(ID_dir(2), ID_dir(1)) >= -45 && atan2d(ID_dir(2), ID_dir(1)) <= 45
            ID_align = 'left';
        elseif atan2d(ID_dir(2), ID_dir(1)) <= -135 || atan2d(ID_dir(2), ID_dir(1)) >= 135
            ID_align = 'right';
        else
            ID_align = 'center';
        end
        % Depict everything
        text(ID_pos_text(1), ID_pos_text(2), sprintf('%02d', iAni), 'Color', [189 042 132]/255, 'Interpreter', 'none', 'HorizontalAlignment', ID_align)
        plot([Tracking.pos_x(idx(iAni)), ID_pos_line(1)],[Tracking.pos_y(idx(iAni)), ID_pos_line(2)], 'Color', [189 042 132]/255, 'LineWidth', 0.5)
    end
    axis equal
    axis off
    xlim([-250, size(frame,1)+250])
    ylim([-250, size(frame,2)+250])
    title(['time: ', sprintf('%03d',iFrame/25), '.0s'])

    % ---------------------------------------------------------------------

    % Get frame
    frame = read(VidObj, iFrame);

    % Get tracking
    idx = find(Tracking.frame == iFrame);
    xPos = floor(Tracking.pos_x(idx));
    yPos = floor(Tracking.pos_y(idx));

    % Preallocation for density
    density = zeros(size(frame,1), size(frame,2));

    % Fill in locations
    for iAni = 1:length(idx)
        density(xPos(iAni),yPos(iAni)) = 1;
    end

    % Smooth
    density_smooth = density;
    for iOrder = 1:2
        density_smooth = imgaussfilt(density_smooth, SET.InteractionRange/(SET.ArenaDiameter/(Annotation.Annotation.ROI.Par(3)*2)));
    end

    % Normalize
    density_smooth = density_smooth/SET.MaxDensity;

    subplot(2,4,[3 4 7 8])
    hold on
    imagesc(frame)
    hDens = imagesc(density_smooth');
    hDens.AlphaData = 0.75;
    colormap(SubFcn.ColMapInferno)
    caxis([0 12])
    hCol = colorbar;
    hCol.Label.String = 'local density';
    % Indicate patches and arena
    th = 0:pi/50:2*pi;
    plot(Annotation.Annotation.ROI.Par(3) * cos(th) + Annotation.Annotation.ROI.Par(1),  Annotation.Annotation.ROI.Par(3) * sin(th) + Annotation.Annotation.ROI.Par(2), 'Color', [0.75 0.75 0.75])
    plot(Annotation.Annotation.Masks.Circular(1,3) * cos(th) + Annotation.Annotation.Masks.Circular(1,1),  Annotation.Annotation.Masks.Circular(1,3) * sin(th) + Annotation.Annotation.Masks.Circular(1,2), 'Color', [0.75 0.75 0.75])
    plot(Annotation.Annotation.Masks.Circular(2,3) * cos(th) + Annotation.Annotation.Masks.Circular(2,1),  Annotation.Annotation.Masks.Circular(2,3) * sin(th) + Annotation.Annotation.Masks.Circular(2,2), 'Color', [0.75 0.75 0.75])
    % Depict Centroids
    mkr = scatter(Tracking.pos_x(idx), Tracking.pos_y(idx),...
        10,[1 1 1],'filled'); mkr.MarkerEdgeColor = 'none'; clear mkr;
    % Show tails
    for iAni = 1:15
        % Get tail
        idx_ani = find(strcmp(Tracking.id, ['A', sprintf('%02d', iAni)]));
        if iFrame<2*VidObj.FrameRate
            idx_tail = idx_ani(1:iFrame);
        else
            idx_tail  = idx_ani(1+iFrame-2*VidObj.FrameRate:iFrame);
        end
        plot(Tracking.pos_x(idx_tail), Tracking.pos_y(idx_tail), 'Color', [1 1 1], 'LineWidth', 0.5)
    end
    axis equal
    axis off
    xlim([-250, size(frame,1)+250])
    ylim([-250, size(frame,2)+250])
    title(['time: ', sprintf('%03d',iFrame/25), '.0s'])

    % Save snapshot
    if cnt_Frame == 1
        export_fig('FIG\raw\mov1', '-pdf')
    end

    % Record video
    F(cnt_Frame) = getframe(gcf) ;
    cnt_Frame = cnt_Frame+1;
    drawnow


    clf
end%iFrame


%% Create the video

% Writer with 50 fps (2x)
writerObj = VideoWriter('ExampleVideo_both.mp4', 'MPEG-4');
writerObj.FrameRate = 50;

% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
close all




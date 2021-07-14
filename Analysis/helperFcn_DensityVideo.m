clear all; close all; clc

%%
SET.VideoName = "D:\DATA\experiment_2_locust_feeding\Data\Videos\15UE20191206.mp4";
VidObj = VideoReader(SET.VideoName);

SET.TrackingName = "D:\DATA\experiment_2_locust_feeding\Data\Tracking_BM\15UE20191206_tracked.csv";
Tracking = readtable(SET.TrackingName);

SET.AnnotationName = "D:\DATA\experiment_2_locust_feeding\Data\Videos\15UE20191206_annotation.mat";
Annotation = load(SET.AnnotationName);

SET.InteractionRange = 7; %cm
SET.ArenaDiameter= 90; %cm

% Max Density
bins = zeros(VidObj.Height, VidObj.Width);
bins(floor(VidObj.Height/2), floor(VidObj.Width/2)) = 1;
for iOrder = 1:2
    bins = imgaussfilt(bins, SET.InteractionRange/(SET.ArenaDiameter/(Annotation.Annotation.ROI.Par(3)*2)));
end%iOrder
SET.MaxDensity = max(max(bins));

%%
figure('Color', 'w')
cnt_Frame = 1;
for iFrame = 1:8:30000
    
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
    
    hold on
    imagesc(frame)
    hDens = imagesc(density_smooth');
    hDens.AlphaData = 0.5;
    colormap(SubFcn.ColMapInferno)
    caxis([0 12])
    hCol = colorbar;
    hCol.Label.String = 'local density';
    axis equal
    axis off
    
    title([num2str(round(iFrame / 25,2)), ' s'])
    
    % Record video
    F(cnt_Frame) = getframe(gcf) ;
    cnt_Frame = cnt_Frame+1;
    drawnow
    
    clf
end%iFrame


% Create the video writer with 50 fps (16x)
writerObj = VideoWriter('Example_Density.avi');
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
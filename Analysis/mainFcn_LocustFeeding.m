%% Information integration for nutritional decision-making in desert locusts
% Swarms of the migratory desert locust can extend over several hundred 
% square kilometres, and starvation compels this ancient pest to devour 
% everything in its path. Theory suggests that gregarious behaviour 
% benefits foraging efficiency over a wide range of spatial food 
% distributions. However, despite the importance of identifying the 
% processes by which swarms locate and select feeding sites to predict 
% their progression, the role of social cohesion during foraging remains 
% elusive. We investigated the evidence accumulation and information 
% integration processes that underlie locusts' nutritional decision-making 
% by employing a Bayesian formalism on high-resolution tracking data from 
% foraging locusts. We tested individual gregarious animals and groups of 
% different sizes in a 2-choice behavioural assay in which food patch 
% qualities were either different or similar. We then predicted the 
% decisions of individual locusts based on personally acquired and socially 
% derived evidence by disentangling the relative contributions of each 
% information class. Our study suggests that locusts balance incongruent 
% evidence but reinforce congruent ones, resulting in more confident 
% assessments when evidence aligns. We provide new insights into the 
% interplay between personal experience and social context in locust 
% foraging decisions which constitute a powerful empirical system to study 
% local individual decisions and their consequent collective dynamics.
%
% This is the main analysis script.
%
% Version: 16-May-2022 (MATLAB R2022a)


% Tidy up
clc; clear all; close all
% Add path to external export function
% (https://github.com/altmany/export_fig)
addpath(genpath('altmany-export_fig'))
% Keep things down
warning('off')



%% Settings
% It is advisable to have a list of the different settings the data
% analysis is based on. This allows to quickly change metrics, as for
% example, calculating the mean instead of the median.


% ***** SCRIPT BEHAVIOUR *****
% Bring settings important for the behaviour of this script all the way to
% the top. This should help to keep an overview of what is happening
% Set whether to pool data again
SET.RecollectData = 1;
% Plot centered, roatated and normalized raw data. Also, if enabled, allow
% the user to adjust thresholds for gap filling (see next settings).
SET.DisplayRawData = 0;


% ***** FILES *****
% Set the path to the folder that contains the tracking results.
SET.Path2Data = 'C:\Users\Yannick\SynologyDrive\Drive\University\PhD\Experiments\experiment_2_locust_feeding\Data\Tracking\';
% Create a list of all data files
foldercontent = dir(SET.Path2Data);
SET.FileList = {foldercontent(~[foldercontent.isdir]).name}';
% Only keep files containing trajectories. Assume each trajectory file
% (file name ending with '_tracked.csv') has a corresponding annotation
% (file name ending with '_annotation.mat'). Excluding the annotations
% helps iterating over data
SET.FileList = SET.FileList(find(cellfun(@isempty, strfind(SET.FileList, '_annotation.mat'))));
% Help identifying data and annotation
SET.DataFileSuffix = '_tracked.csv';
SET.AnnotationFileSuffix = '_annotation.mat';
% Set patch to table containing before???after weights
SET.PatchWeights = 'C:\Users\Yannick\SynologyDrive\Drive\University\PhD\Experiments\experiment_2_locust_feeding\Data\PatchWeights.csv';
PatchWeightsTable = readtable(SET.PatchWeights);

% ***** INFO ON DATA COLLECTION *****
% Set frame rate of recordings
SET.FrameRate = 25; %[fps]
% Set data length (some trials are longer, some are shorter than 1h by a
% few frames)
SET.CutAfter = 30; %[min]
SET.CutAfter = SET.CutAfter*60*SET.FrameRate;
% Set arena diameter
SET.dArena = 90; %[cm]
% Angles of patch positions
SET.PatchAngle = [-90; -150; 150; 90; 30; -30];


% ***** DATA ANALYSIS *****
% Set additional distance to patch edge as tolerance
SET.ToleranceDistance = 4; %[4cm = 1BL]
% Set bins in which data will be divided for the heatmaps. Resulting number
% of bins will be SET.HeatMapGrid^2.
SET.HeatMapGrid = 1000;
% Bootsrapping
SET.BootSamples = 5000;
% Interaction range for density filter
SET.InteractionRange_SD = 7; %(cm)
% Std for imgausfilt()
SET.SmoothValue = (SET.HeatMapGrid/SET.dArena)*SET.InteractionRange_SD;
% The minimum inter-bout-interval and bout duration
SET.minInterBoutInterval = 5; %(s, animal-centric)
SET.minBoutDuration = 5; %(s)
SET.minInterBoutInterval = SET.minInterBoutInterval*SET.FrameRate;
SET.minBoutDuration = SET.minBoutDuration*SET.FrameRate;
% Determine whether an animal is walking or not
SET.MotionThreshold = 0.25;

% ***** DISPLAY OF RESULTS *****
% Colors
% --- UE, A (blue)
SET.Color.UE_2.N01 = [084 132 255]/255;
SET.Color.UE_2.N05 = [063 131 222]/255;
SET.Color.UE_2.N10 = [042 130 189]/255;
SET.Color.UE_2.N15 = [021 128 156]/255;
SET.Color.UE_2.N30 = [000 127 123]/255;
% --- UE, B (pink)
SET.Color.UE_1.N01 = [255 084 138]/255;
SET.Color.UE_1.N05 = [222 063 135]/255;
SET.Color.UE_1.N10 = [189 042 132]/255;
SET.Color.UE_1.N15 = [156 021 130]/255;
SET.Color.UE_1.N30 = [132 000 127]/255;
% --- EQ, A (orange)
SET.Color.EQ_2.N01 = [255 132 084]/255;
SET.Color.EQ_2.N05 = [223 099 064]/255;
SET.Color.EQ_2.N10 = [191 066 044]/255;
SET.Color.EQ_2.N15 = [158 033 024]/255;
SET.Color.EQ_2.N30 = [126 000 004]/255;
% --- EQ, B (green)
SET.Color.EQ_1.N01 = [084 255 199]/255;
SET.Color.EQ_1.N05 = [064 223 149]/255;
SET.Color.EQ_1.N10 = [044 191 100]/255;
SET.Color.EQ_1.N15 = [024 159 050]/255;
SET.Color.EQ_1.N30 = [004 127 000]/255;
% Colormap for heatmaps
SET.ColorHeat = SubFcn.ColMapInferno(1000);
% Set names of conditions
SET.ConditionNames.Patch = {'EQ', 'UE'};
SET.ConditionNames.Group = {'N01', 'N05', 'N10', 'N15', 'N30'};
SET.FigureAppendix = '_raw';


%% Collect daand pool data

% Only collect data if set by the user. If data has already been pooled and
% saved, simply loading it will save time
if SET.RecollectData == 1
    
    % Indicate progress
    h = waitbar(0, {'Please wait while we collect data ...'; ['file: 0/', num2str(length(SET.FileList))]});
    
    % Iterate over all data files
    for iFile = 1:length(SET.FileList)
        
        % Waitbar
        waitbar(iFile/length(SET.FileList), h, {'Please wait while we collect data ...'; ['file: ', num2str(iFile), '/', num2str(length(SET.FileList))]});
        
        % Keep information about the current data file in the
        % variale "currFile"
        currFile.FileName = SET.FileList{iFile}(1:end-length(SET.DataFileSuffix));
        currFile.Condition = currFile.FileName(3:4); % Uniform labeling of files allows us to extract the information we need.
        currFile.N_str = ['N', currFile.FileName(1:2)]; % Uniform labeling of files allows us to extract the information we need.
        currFile.N = str2double(currFile.FileName(1:2)); % Uniform labeling of files allows us to extract the information we need.
        
        % Load data
        currFile.Data = readtable([SET.Path2Data, [currFile.FileName, SET.DataFileSuffix]]);
        
        % Load annotation
        currFile.Annotation = load([SET.Path2Data, [currFile.FileName, SET.AnnotationFileSuffix]]);
        
        % Save annotation
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc = currFile.Annotation.PatchA.loc;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Radius = currFile.Annotation.PatchA.radius;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.DetectionRadius = currFile.Annotation.PatchA.radius + SET.ToleranceDistance/currFile.Annotation.Arena.radius_cm;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc = currFile.Annotation.PatchB.loc;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Radius = currFile.Annotation.PatchB.radius;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.DetectionRadius = currFile.Annotation.PatchB.radius + SET.ToleranceDistance/currFile.Annotation.Arena.radius_cm;
        
        % Get how much food was consumed
        idx = find(strcmp(PatchWeightsTable.id, currFile.FileName));
        % Check whether patches are correctly assigned to letters A & B
        switch currFile.Condition
            case 'UE'
                if (PatchWeightsTable.A_cond(idx) == 1) & (PatchWeightsTable.B_cond(idx) == -1)
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Consumption = ...
                        PatchWeightsTable.A_before(idx) - PatchWeightsTable.A_after(idx);
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Consumption = ...
                        PatchWeightsTable.B_before(idx) - PatchWeightsTable.B_after(idx);
                elseif (PatchWeightsTable.A_cond(idx) == -1) & (PatchWeightsTable.B_cond(idx) == 1)
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Consumption = ...
                        PatchWeightsTable.B_before(idx) - PatchWeightsTable.B_after(idx);
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Consumption = ...
                        PatchWeightsTable.A_before(idx) - PatchWeightsTable.A_after(idx);
                elseif isnan(PatchWeightsTable.A_cond(idx)) || isnan(PatchWeightsTable.B_cond(idx))
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Consumption = NaN;
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Consumption = NaN;
                else
                    error('Patch identities do not match labels')
                end
            case 'EQ'
                % Get orientation of patches
                phi_A = atan2d(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc(2), DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc(1)) + (currFile.Annotation.Transformation.rotation_deg*-1);
                phi_B = atan2d(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc(2), DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc(1)) + (currFile.Annotation.Transformation.rotation_deg*-1);
                if phi_A>180
                    phi_A = -360+phi_A;
                elseif phi_A<-180
                    phi_A = 360+phi_A;
                end
                if phi_B>180
                    phi_B = -360+phi_B;
                elseif phi_B<-180
                    phi_B = 360+phi_B;
                end
                [~, idx_A] = min(abs(SET.PatchAngle-phi_A));
                [~, idx_B] = min(abs(SET.PatchAngle-phi_B));
                
                if (PatchWeightsTable.A_loc(idx) == idx_A) & (PatchWeightsTable.B_loc(idx) == idx_B)
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Consumption = ...
                        PatchWeightsTable.A_before(idx) - PatchWeightsTable.A_after(idx);
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Consumption = ...
                        PatchWeightsTable.B_before(idx) - PatchWeightsTable.B_after(idx);
                elseif (PatchWeightsTable.A_loc(idx) == idx_B) & (PatchWeightsTable.B_loc(idx) == idx_A)
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Consumption = ...
                        PatchWeightsTable.B_before(idx) - PatchWeightsTable.B_after(idx);
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Consumption = ...
                        PatchWeightsTable.A_before(idx) - PatchWeightsTable.A_after(idx);
                elseif isnan(PatchWeightsTable.A_cond(idx)) || isnan(PatchWeightsTable.B_cond(idx))
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Consumption = NaN;
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Consumption = NaN;
                else
                    error('Patch identities do not match labels')
                end
        end
        clear idx*
        
        % Get number of unique animals IDs
        uniqueAnimals = unique(currFile.Data.id);
        
        % Check whether there is a mismatch with the filename
        if length(uniqueAnimals) ~= currFile.N
            str = ['Mismatch between number of animals specified in file name and number of unique IDs found. ', ...
                'File: ', currFile.FileName];
            warning(str)
            clear str
        end %if mismatch in number of animals
        
        % Preallocation
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x =                      nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_y =                      nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).rawSpeed =                   nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion =         nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion_walk =    nan(length(uniqueAnimals),1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion_stand =   nan(length(uniqueAnimals),1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A =               nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B =               nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A =                  nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B =                  nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingInterval =            nan(length(uniqueAnimals),1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).NonFeedingInterval =         nan(length(uniqueAnimals),1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FoundPatch_A =               nan(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FoundPatch_B =               nan(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_A =                  nan(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_B =                  nan(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_A =           nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_B =           nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_Both =        nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Switching_t1 =               zeros(2, 2);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Switching_t2 =               zeros(4, 2);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AffectionIndex =             nan(SET.CutAfter, length(uniqueAnimals));
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_A =      cell(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_B =      cell(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PreferenceIndex =            nan(SET.CutAfter, 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A =             [];
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B =             [];
        
        % Iterate over all animals
        for iAni = 1:length(uniqueAnimals)
            
            % Get the indices for the current animal
            idxAnimals = find(strcmp(currFile.Data.id, uniqueAnimals{iAni}));
            
            % Get the current animals positional data
            pos_x = currFile.Data.pos_x(idxAnimals);
            pos_y = currFile.Data.pos_y(idxAnimals);
            
            % Each animal gets its own column
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x(:, iAni) = pos_x;
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_y(:, iAni) = pos_y;
            
            % Also extract raw speed profiles for each animal. This could
            % later be used to detemine jumps
            helper = [0, 0; diff(pos_x), diff(pos_y)];
            speed = sqrt(sum(helper'.*helper'))';
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).rawSpeed(:, iAni) = speed*SET.FrameRate*(SET.dArena/2);
            clear speed
            
            % Determine whether animal is walking or standing
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion(:, iAni) = ...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).rawSpeed(:, iAni) > SET.MotionThreshold;
            
            % Determine whether an animal is at the food patch A or not
            % --- Calculate the distance to the patch A
            helper = [pos_x, pos_y] - DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc;
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A(:, iAni) = sqrt(sum(helper'.*helper'))';
            % --- Determine whether they are at the patch
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni) = ...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A(:, iAni) < DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.DetectionRadius;
            
            % Determine whether an animal is at the food patch B or not
            % --- Calculate the distance to the patch B
            helper = [pos_x, pos_y] - DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc;
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B(:, iAni) = sqrt(sum(helper'.*helper'))';
            % --- Determine whether they are at the patch
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni) = ...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B(:, iAni) < DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.DetectionRadius;
            
            % Correct visits. Consecutive visits that have an
            % inter-vist-interval shorter than SET.minInterBoutInterval
            % will be combined. Visits shorter than SET.minBoutDuration
            % will be kicked out
            % --- Patch A
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni) = SubFcn.CorrectVisits(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni), SET);
            % --- Patch B
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni) = SubFcn.CorrectVisits(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni), SET);
            
            % Determine visit history
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).VisitHistory(:, iAni) = SubFcn.VisitHistory(...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni),...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni));
            
            % Determine walking and standing intervals
            % --- walk
            helper = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion(:, iAni);
            helper(logical(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)+DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)))=0;
            temp = regionprops(bwlabel(helper));
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion_walk(iAni,1) = ...
                mean([temp.Area]/SET.FrameRate);
            % --- stand
            helper = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion(:, iAni)<1;
            helper(logical(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)+DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)))=0;
            temp = regionprops(bwlabel(helper));
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion_stand(iAni,1) = ...
                mean([temp.Area]/SET.FrameRate);
            
            % Determine feeding intervals
            temp1 = bwlabel((DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni) + DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni))>0);
            uniqueInterval = unique(temp1);
            if max(uniqueInterval)~=0
                temp2 = [];
                uniqueInterval(uniqueInterval==0)=[];
                for iInterval = 1:length(uniqueInterval)
                    idx = find(temp1 == uniqueInterval(iInterval));
                    if idx(1)~=1 && idx(end)~=SET.CutAfter
                        temp2 = [temp2; length(idx)/SET.FrameRate];
                    end%if
                end%iInterval
            else
                temp2 = NaN;
            end%if
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingInterval(iAni,1) = nanmean(temp2);
            
            % Determine non-feeding intervals
            temp1 = bwlabel((DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni) + DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni))<1);
            uniqueInterval = unique(temp1);
            if max(uniqueInterval)~=0
                temp2 = [];
                uniqueInterval(uniqueInterval==0)=[];
                for iInterval = 1:length(uniqueInterval)
                    idx = find(temp1 == uniqueInterval(iInterval));
                    if idx(1)~=1 && idx(end)~=SET.CutAfter
                        temp2 = [temp2; length(idx)/SET.FrameRate];
                    end%if
                end%iInterval
            else
                temp2 = NaN;
            end%if
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).NonFeedingInterval(iAni,1) = nanmean(temp2);
            
            % Determine how fast animals were finding the patches
            % --- Patch A
            if ~isempty(find(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni), 1))
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FoundPatch_A(iAni, 1) = find(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni), 1)/SET.FrameRate;
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FoundPatch_A(iAni, 1) = NaN;
            end
            % --- Patch B
            if ~isempty(find(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni), 1))
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FoundPatch_B(iAni, 1) = find(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni), 1)/SET.FrameRate;
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FoundPatch_B(iAni, 1) = NaN;
            end
            
            % Determine the number of visits to either patch
            % --- Patch A
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_A(iAni, 1) = length(find(diff(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)) == 1));
            % --- Patch B
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_B(iAni, 1) = length(find(diff(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)) == 1));
            
            % Determine the instantaneous visiting frequency, i.e. the
            % interval inbetween consecutive visits
            % --- Patch A
            if length(unique(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni))) > 1
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_A(:, iAni) = SubFcn.GetInstaFreq(find(diff(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)) == 1), SET.CutAfter, SET.FrameRate);
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_A(:, iAni) = zeros(SET.CutAfter, 1);
            end
            % --- Patch B
            if length(unique(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni))) > 1
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_B(:, iAni) = SubFcn.GetInstaFreq(find(diff(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)) == 1), SET.CutAfter, SET.FrameRate);
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_B(:, iAni) = zeros(SET.CutAfter, 1);
            end
            % --- Both
            temp = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni) + DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni);
            if length(unique(temp)) > 1
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_Both(:, iAni) = SubFcn.GetInstaFreq(find(diff(temp) == 1), SET.CutAfter, SET.FrameRate);
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).InstaVisitFreq_Both(:, iAni) = zeros(SET.CutAfter, 1);
            end
            clear temp
            
            % Determine switching probibilities
            % --- First, taking only the stay immedeatly prevous stay into
            % account
            helper_t1 = SubFcn.GetSwitching_t1(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni), DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni));
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Switching_t1 = ...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Switching_t1 + helper_t1;
            % --- Second, now take one more time step into account
            helper_t2 = SubFcn.GetSwitching_t2(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni), DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni));
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Switching_t2 = ...
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Switching_t2 + helper_t2;
            clear helper_*
            
            % Introduce an index reflecting the general interest for
            % the food patches. This index is based on on the distance
            % to the patches: -(A-B)/(A+B), for A being the patch labeled
            % with condition 1. The closer the index is to positive
            % one, the higher the affection for patch A. Similarly for
            % patch B, the closer the index is to negative one, the higher
            % the interest for this patch.
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AffectionIndex(:, iAni) = ...
                -(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A(:, iAni) - DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B(:, iAni)) ./ (DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A(:, iAni) + DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B(:, iAni));
            
            % Get how long on average an animal was feeding at patch A
            intervals = bwlabel(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni));
            intervals_unique = unique(intervals);
            intervals_unique(intervals_unique == 0) = [];
            interval_dur = regionprops(intervals);
            interval_dur = [interval_dur(:).Area];
            if ~isempty(interval_dur)
                % Get starting index, stopping index, and duration (in s) of feeding bout
                helper = [];
                for iStay = 1:length(intervals_unique)
                    helper = [...
                        helper;...
                        find(intervals == intervals_unique(iStay), 1, 'first'), find(intervals == intervals_unique(iStay), 1, 'last'), interval_dur(iStay)/SET.FrameRate];
                end%iStay
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_A{iAni, 1} = helper;
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_A{iAni, 1} = [];
            end
            clear intervals intervals_unique interval_dur helper
            
            % Get how long on average an animal was feeding at patch B
            intervals = bwlabel(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni));
            intervals_unique = unique(intervals);
            intervals_unique(intervals_unique == 0) = [];
            interval_dur = regionprops(intervals);
            interval_dur = [interval_dur(:).Area];
            if ~isempty(interval_dur)
                % Get starting index, stopping index, and duration (in s) of feeding bout
                helper = [];
                for iStay = 1:length(intervals_unique)
                    helper = [...
                        helper;...
                        find(intervals == intervals_unique(iStay), 1, 'first'), find(intervals == intervals_unique(iStay), 1, 'last'), interval_dur(iStay)/SET.FrameRate];
                end%iStay
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_B{iAni, 1} = helper;
            else
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_B{iAni, 1} = [];
            end
            clear intervals intervals_unique interval_dur helper
            
            % Survival Curves
            % --- Patch A
            area = struct2array(regionprops(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:,iAni) > 0, 'area'));
            for iStay = 1:length(area)
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A(end+1,1:area(iStay)) = 1;
            end
            clear area
            % --- Patch B
            area = struct2array(regionprops(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:,iAni) > 0, 'area'));
            for iStay = 1:length(area)
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B(end+1,1:area(iStay)) = 1;
            end
            clear area
            % --- Correct size
            if size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A, 2) < SET.CutAfter
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A = [...
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A,...
                    zeros(size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A,1),SET.CutAfter-size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_A,2))];
            end
            if size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B, 2) < SET.CutAfter
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B = [...
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B,...
                    zeros(size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B,1),SET.CutAfter-size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).SurvivalTime_B,2))];
            end
            
        end %iAni
        
        % Bin 2D data
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_2D = ...
            hist3([currFile.Data.pos_x, currFile.Data.pos_y], 'ctrs', {linspace(-1, 1, SET.HeatMapGrid) linspace(-1, 1, SET.HeatMapGrid)})/(size(currFile.Data,1)/currFile.N);
        
        % PreferenceIndex
        helper_A = sum(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x>0, 2);
        helper_B = sum(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x<0, 2);
        PreferenceIndex = (helper_A-helper_B)/(size(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x, 2));
        PreferenceIndex(isnan(PreferenceIndex)) = 0;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PreferenceIndex(:, 1) = PreferenceIndex;
        clear helper_* PreferenceIndex
        
        % PatchPreferenceIndex
        if currFile.N > 1
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchPreferenceIndex(:, 1) = ...
                (sum(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A, 2) - sum(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B, 2)) ...
                ./ (sum(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A, 2) + sum(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B, 2));
        else
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchPreferenceIndex(:, 1) = ...
                (DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A) - (DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B);
        end
        
        % Local density at food patches and animals
        % --- Patch A
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_A = sum(normpdf(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2)), 2);
        % --- Patch B
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_B = sum(normpdf(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2)), 2);
        % --- Animals
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Density = zeros(SET.CutAfter, currFile.N);
        if currFile.N>1
            for iFrame = 1:SET.CutAfter
                % Get data
                X = [DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x(iFrame, :)', ...
                    DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_y(iFrame, :)'];
                % Compute the Euclidean distance.
                D = pdist(X);
                % Arrange in squareform
                Z = squareform(D);
                % Get animal-centric density
                InterAniDensity = normpdf(Z, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Density(iFrame, :) = sum(InterAniDensity, 2);
                clear X D Z InterAniDensity
            end%iFrame
        end%if more than 1 animal
        
        
        % If set by the user, show trajectories
        if SET.DisplayRawData
            hFig = figure('Name', currFile.FileName, 'Color', 'w');
            % --- Traj ---
            subplot(1, 2, 1); hold on
            set(gca, 'Color', 'k')
            % Set colors
            col = SubFcn.ColMapInferno(currFile.N);
            if strcmp(currFile.Condition, 'EQ')
                col = col(:, [2, 1, 3]);
            end
            % Plot trajectory for each animal
            for iAni = 1:length(uniqueAnimals)
                plot(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x(:, iAni), DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_y(:, iAni), 'Color', col(iAni, :))
            end
            % Indicate food patches
            SubFcn.DrawCircle(... --- A ---
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc(1), ... x
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc(2), ... y
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Radius, ... r
                [.75 .75 .75]);% color
            SubFcn.DrawCircle(... --- B ---
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc(1), ... x
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc(2), ... y
                DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Radius, ... r
                [.75 .75 .75]);% color
            % Indicate arena
            SubFcn.DrawCircle(... --- B ---
                0, 0, 1, [.75 .75 .75]);
            axis equal
            xlim([-1 1])
            ylim([-1 1])
            xticks([])
            yticks([])
            box on
            % --- Heatmap ---
            subplot(1, 2, 2); hold on
            baseline_max = zeros(SET.HeatMapGrid,SET.HeatMapGrid);
            baseline_max(floor(SET.HeatMapGrid/2), floor(SET.HeatMapGrid/2))=1;
            baseline_max = imgaussfilt(baseline_max, SET.SmoothValue, 'FilterDomain', 'spatial', 'FilterSize', SET.HeatMapGrid-1);
            baseline_max = max(max(baseline_max));
            img = imgaussfilt(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_2D, SET.SmoothValue, 'FilterDomain', 'spatial')/baseline_max;
            imagesc(flipud(img))
            axis equal
            xlim([1 , SET.HeatMapGrid])
            ylim([1 , SET.HeatMapGrid])
            xticks([])
            yticks([])
            set(gca, 'view', [-90 90])
            box on
            colormap(SET.ColorHeat)
            % Adjust range in col map
            caxis([0, quantile(reshape(img, [1, numel(img)]), 0.99)]);
            colorbar
            % Save figure
            if exist('export_fig') > 0
                if ~isdir('FIG\raw\Trials'); mkdir('FIG\raw\Trials'); end
                export_fig(['FIG\raw\Trials\', currFile.FileName, '_rawData'], '-pdf')
            else
                warning('The external function export_fig has not been added to the path or does not exist!')
                warning('Figures will not be saved.')
            end
            close(hFig)
            clear hFig ax col img xvec yvec
        end%if display
        
    end %iFile
    
    
    
    
    %% ***** POOL DATA *****
    % Iterate over patch conditions: EQ and UE
    for iPatch = 1:length(SET.ConditionNames.Patch)
        
        % Iterate over group sizes: 1, 5, 10, 15, 30
        for iGroup = 1:length(SET.ConditionNames.Group)
            
            % Pre-specification
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Speed = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_walk = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_stand = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NonFeedingInterval = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PreferenceIndex = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AffectionIndex = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PatchPreferenceIndex = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Time2Feeding = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PropTimeFeeding = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).N_Feeding_TC = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NumOfVisits = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_Both = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t1 = zeros(2,2);
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t2 = zeros(4,2);
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AniCentric_FeedingInterval_GrpSize_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AniCentric_FeedingInterval_GrpSize_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).LeavingInterval_GrpSize_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).LeavingInterval_GrpSize_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).SurvivalCurves_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).SurvivalCurves_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBouts = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Splitting = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AmountFeeding_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AmountFeeding_B = [];
            
            % Get corresponding trials
            currTrialList = fieldnames(DATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}));
            
            % Iterate over all trials recorded under this combination of
            % conditions
            for iTrial = 1:length(currTrialList)
                
                % Shortcut
                currData = DATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).(currTrialList{iTrial});
                
                % --- Speed
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Speed = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Speed;...
                    nanmean(currData.rawSpeed,2)'];
                
                % --- intermittentMotion_walk
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_walk = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_walk;...
                    nanmean(currData.intermittentMotion_walk)];
                
                % --- intermittentMotion_stand
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_stand = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_stand;...
                    nanmean(currData.intermittentMotion_stand)];
                
                % --- FeedingInterval
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval;...
                    nanmean(currData.FeedingInterval)];
                
                % --- NonFeedingInterval
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NonFeedingInterval = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NonFeedingInterval;...
                    nanmean(currData.NonFeedingInterval)];
                
                % --- PreferenceIndex
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PreferenceIndex = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PreferenceIndex;...
                    nanmean(currData.PreferenceIndex,2)'];
                
                % --- AffectionIndex
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AffectionIndex = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AffectionIndex;...
                    nanmean(currData.AffectionIndex, 2)'];
                
                % --- PatchPreferenceIndex
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PatchPreferenceIndex = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PatchPreferenceIndex;...
                    nanmean(currData.PatchPreferenceIndex, 2)'];
                
                % --- Time2Feeding
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Time2Feeding = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Time2Feeding;...
                    nanmean(currData.FoundPatch_A), nanmean(currData.FoundPatch_B)];
                
                % --- PropTimeFeeding
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PropTimeFeeding = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PropTimeFeeding;...
                    sum(sum(currData.AtPatch_A))/sum(sum([currData.AtPatch_A;currData.AtPatch_B])), sum(sum(currData.AtPatch_B))/sum(sum([currData.AtPatch_A;currData.AtPatch_B]))];
                
                % --- N_Feeding_TC
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).N_Feeding_TC = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).N_Feeding_TC;...
                    (sum(currData.AtPatch_A, 2)+sum(currData.AtPatch_B, 2))'];
                
                % --- NumOfVisits
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NumOfVisits = [
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NumOfVisits;...
                    nanmean(currData.N_Visit_A), nanmean(currData.N_Visit_B)];
                
                % --- InstaVisitFreq_A
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_A = [
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_A;...
                    nanmean(currData.InstaVisitFreq_A,2)'];
                
                % --- InstaVisitFreq_B
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_B = [
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_B;...
                    nanmean(currData.InstaVisitFreq_B,2)'];
                
                % --- InstaVisitFreq_Both
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_Both = [
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).InstaVisitFreq_Both;...
                    nanmean(currData.InstaVisitFreq_Both,2)'];
                
                % --- Switching_t1
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t1 = ...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t1 + currData.Switching_t1;
                
                % --- Switching_t2
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t2 = ...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t2 + currData.Switching_t2;
                
                % --- FeedingInterval_GrpSize_A (col1: interval, col2: density at entry)
                helper_pop_A = sum(currData.AtPatch_A, 2);
                helper_pop_B = sum(currData.AtPatch_B, 2);
                helper_diff = diff(helper_pop_A);
                if range(helper_diff) ~= 0
                    idx = find(helper_diff>0); idx(1) = [];
                    N = size(currData.pos_x, 2); %0, number of animals in the arena
                    helper_interval = diff(idx)/SET.FrameRate; %1, interval
                    helper_pop_A = helper_pop_A(idx(2:end)); %2, grp size at entry
                    helper_pop_B = helper_pop_B(idx(2:end));
                    helper_density = currData.PatchDensity_A(idx(2:end)); %4, density at entry
                    % Correct density (subract the density imposed by joining
                    % animal)
                    helper_density = helper_density - normpdf(currData.PatchA.DetectionRadius, 0, SET.InteractionRange_SD/(SET.dArena/2)) / normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                    % Put everything together
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_A = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_A;...
                        helper_interval, helper_density];
                end
                clear helper_* idx N
                
                % --- FeedingInterval_GrpSize_B (col1: interval, col2: density at entry)
                helper_pop_A = sum(currData.AtPatch_A, 2);
                helper_pop_B = sum(currData.AtPatch_B, 2);
                helper_diff = diff(helper_pop_B);
                if range(helper_diff) ~= 0
                    idx = find(helper_diff>0); idx(1) = [];
                    N = size(currData.pos_x, 2); %0, number of animals in the arena
                    helper_interval = diff(idx)/SET.FrameRate; %1, interval
                    helper_pop_A = helper_pop_A(idx(2:end));
                    helper_pop_B = helper_pop_B(idx(2:end)); %2, grp size at entry
                    helper_density = currData.PatchDensity_B(idx(2:end)); %4, density at entry
                    % Correct density (subract the density imposed by joining
                    % animal)
                    helper_density = helper_density - normpdf(currData.PatchB.DetectionRadius, 0, SET.InteractionRange_SD/(SET.dArena/2)) / normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                    % Put everything together
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_B = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_B;...
                        helper_interval, helper_density];
                end
                clear helper_* idx N
                
                % --- AniCentric_FeedingInterval_GrpSize_A
                for iAni = 1:size(currData.pos_x,2)
                    clear helper*
                    helper_diff_A = diff(currData.AtPatch_A(:,iAni));
                    helper_diff = diff(currData.AtPatch_A(:,iAni)+currData.AtPatch_B(:,iAni));
                    if ~isempty(find(helper_diff_A==1))
                        idx = find(helper_diff_A>0);
                        idx_leave = find(helper_diff<0);
                        density = currData.PatchDensity_A(idx) - normpdf(currData.Dist2Patch_A(idx), 0, SET.InteractionRange_SD/(SET.dArena/2)) / normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                        for iJoin = 1:length(idx)
                            if idx(iJoin) == min(find(abs(helper_diff)))
                                helper_interval(iJoin,1) = (idx(iJoin))/SET.FrameRate;
                            else
                                helper = find(helper_diff==-1)-(idx(iJoin));
                                helper_interval(iJoin,1) = abs(max(helper(helper<0)))/SET.FrameRate;
                            end
                            helper_density(iJoin,1) = density(iJoin);
                            helper_history(iJoin,1) = currData.VisitHistory(idx(iJoin),iAni);
                        end%iJoin
                        % Put everything together
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AniCentric_FeedingInterval_GrpSize_A = [...
                            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AniCentric_FeedingInterval_GrpSize_A;...
                            helper_interval, helper_density, helper_history];
                    end%if
                end%iAni
                clear helper_* idx* density
                
                % --- AniCentric_FeedingInterval_GrpSize_B
                for iAni = 1:size(currData.pos_x,2)
                    clear helper*
                    helper_diff_B = diff(currData.AtPatch_B(:,iAni));
                    helper_diff = diff(currData.AtPatch_A(:,iAni)+currData.AtPatch_B(:,iAni));
                    if ~isempty(find(helper_diff_B==1))
                        idx = find(helper_diff_B>0);
                        idx_leave = find(helper_diff<0);
                        density = currData.PatchDensity_B(idx) - normpdf(currData.Dist2Patch_B(idx), 0, SET.InteractionRange_SD/(SET.dArena/2)) / normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                        for iJoin = 1:length(idx)
                            if idx(iJoin) == min(find(abs(helper_diff)))
                                helper_interval(iJoin,1) = (idx(iJoin))/SET.FrameRate;
                            else
                                helper = find(helper_diff==-1)-(idx(iJoin));
                                helper_interval(iJoin,1) = abs(max(helper(helper<0)))/SET.FrameRate;
                            end
                            helper_density(iJoin,1) = density(iJoin);
                            helper_history(iJoin,1) = currData.VisitHistory(idx(iJoin),iAni);
                        end%iJoin
                        % Put everything together
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AniCentric_FeedingInterval_GrpSize_B = [...
                            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AniCentric_FeedingInterval_GrpSize_B;...
                            helper_interval, helper_density, helper_history];
                    end%if
                end%iAni
                clear helper_* idx* density
                
                % --- LeavingInterval_GrpSize_A (col1: interval, col2: density at entry)
                helper_pop_A = sum(currData.AtPatch_A, 2);
                helper_diff = diff(helper_pop_A);
                if range(helper_diff) ~= 0
                    idx_join = find(helper_diff>0);
                    idx_leave = find(helper_diff<0);
                    helper_interval = nan(length(idx_join),1);
                    for iJoin = 1:length(idx_join)
                        helper_1 = idx_leave-idx_join(iJoin);
                        helper_2 = min(helper_1(helper_1>0));
                        if ~isnan(helper_2)
                            helper_interval(iJoin,1) = helper_2/SET.FrameRate;
                        end%if
                    end%iJoin
                    helper_density = currData.PatchDensity_A(idx_join); %2, density at entry
                    % Put everything together
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).LeavingInterval_GrpSize_A = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).LeavingInterval_GrpSize_A;...
                        helper_interval(~isnan(helper_interval)), helper_density(~isnan(helper_interval))];
                end
                clear helper_* idx_*
                
                % --- LeavingInterval_GrpSize_B (col1: interval, col2: density at entry)
                helper_pop_B = sum(currData.AtPatch_B, 2);
                helper_diff = diff(helper_pop_B);
                if range(helper_diff) ~= 0
                    idx_join = find(helper_diff>0);
                    idx_leave = find(helper_diff<0);
                    helper_interval = nan(length(idx_join),1);
                    for iJoin = 1:length(idx_join)
                        helper_1 = idx_leave-idx_join(iJoin);
                        helper_2 = min(helper_1(helper_1>0));
                        if ~isnan(helper_2)
                            helper_interval(iJoin,1) = helper_2/SET.FrameRate;
                        end%if
                    end%iJoin
                    helper_density = currData.PatchDensity_B(idx_join); %2, density at entry
                    % Put everything together
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).LeavingInterval_GrpSize_B = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).LeavingInterval_GrpSize_B;...
                        helper_interval(~isnan(helper_interval)), helper_density(~isnan(helper_interval))];
                end
                clear helper_* idx_*
                
                % SurvivalCurves_A
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).SurvivalCurves_A = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).SurvivalCurves_A;...
                    nanmean(currData.SurvivalTime_A,1)];
                
                % SurvivalCurves_B
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).SurvivalCurves_B = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).SurvivalCurves_B;...
                    nanmean(currData.SurvivalTime_B,1)];
                
                % --- FeedingBout_GrpSize_A
                helper_pop = sum(currData.AtPatch_A, 2);
                for iAni = 1:str2double(SET.ConditionNames.Group{iGroup}(2:end))
                    helper_stay = currData.FeedingBoutDuration_A{iAni, 1};
                    if ~isempty(helper_stay)
                        for iStay = 1:size(helper_stay, 1)
                            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_A = [...
                                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_A;...
                                helper_stay(iStay, 3), nanmean(helper_pop(helper_stay(iStay, 1):helper_stay(iStay, 2))), nanmean(currData.Density(helper_stay(iStay, 1):helper_stay(iStay, 2), iAni))];
                        end%iStay
                    end%if
                end%iAni
                clear helper_pop helper_stay
                
                % --- FeedingBout_GrpSize_B
                helper_pop = sum(currData.AtPatch_B, 2);
                for iAni = 1:str2double(SET.ConditionNames.Group{iGroup}(2:end))
                    helper_stay = currData.FeedingBoutDuration_B{iAni, 1};
                    if ~isempty(helper_stay)
                        for iStay = 1:size(helper_stay, 1)
                            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_B = [...
                                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_B;...
                                helper_stay(iStay, 3), nanmean(helper_pop(helper_stay(iStay, 1):helper_stay(iStay, 2))), nanmean(currData.Density(helper_stay(iStay, 1):helper_stay(iStay, 2), iAni))];
                        end%iStay
                    end%if
                end%iAni
                clear helper_pop helper_stay
                
                % --- FeedingBouts
                helper = currData.AtPatch_A' - currData.AtPatch_B';
                for iAni = 1:size(helper, 1)
                    Vq = interp1(1:size(helper, 2), helper(iAni, :), linspace(1, size(helper, 2), 1000), 'nearest');
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBouts = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBouts;...
                        Vq];
                end%iAni
                clear helpre Vq
                
                % --- Splitting
                helper = [...
                    (currData.PatchDensity_A+currData.PatchDensity_B), currData.PatchDensity_A./(currData.PatchDensity_A+currData.PatchDensity_B);...
                    (currData.PatchDensity_A+currData.PatchDensity_B), currData.PatchDensity_B./(currData.PatchDensity_A+currData.PatchDensity_B)];
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Splitting = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Splitting;...
                    helper];
                
                % --- AmountFeeding_A
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AmountFeeding_A = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AmountFeeding_A;...
                    sum(currData.PatchDensity_A), currData.PatchA.Consumption];
                
                % --- AmountFeeding_B
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AmountFeeding_B = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).AmountFeeding_B;...
                    sum(currData.PatchDensity_B), currData.PatchB.Consumption];
                
            end% iTrial
            
            % Sort FeedingBouts
            helper = PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBouts;
            idxAni = zeros(1, size(helper, 1));
            for iAni = 1:size(helper, 1)
                if isempty(find(helper(iAni, :) ~= 0, 1))
                    idxAni(iAni) = inf;
                else
                    idxAni(iAni) = find(helper(iAni, :) ~= 0, 1);
                end
            end%iAni
            [~, idxSort] = sort(idxAni);
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBouts = helper(idxSort, :);
            clear idx* helper
            
        end% iN
    end% iCod
    
    % Save pooled data
    save('PooledData.mat', 'DATA', 'PooledDATA', 'SET', '-v7.3')
    
    % Close waitbar
    close(h)
    
else
    
    % Load previously collected data
    load('PooledData.mat')
    
end% if recollect data

% Tidy up
clearvars -except DATA PooledDATA SET


%% Plot data
clear hFig
% Start parallel pool
gcp;
%--------------------------------------------------------------------------
hFig.ExampleDensity =                   figure('Name', 'ExampleDensity',        'units', 'normalized', 'Position', [0.25 0.25 0.25 0.25], 'Color', 'w');
hFig.AmountFeeding =                    figure('Name', 'AmountFeeding',         'units', 'normalized', 'Position', [0.25 0.25 0.25 0.25], 'Color', 'w');
hFig.Speed =                            figure('Name', 'Speed',                 'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.PreferenceIndex =                  figure('Name', 'PreferenceIndex',       'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.AffectionIndex =                   figure('Name', 'AffectionIndex',        'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.PatchPreferenceIndex =             figure('Name', 'PatchPreferenceIndex',  'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.FeedingBouts_EQ =                  figure('Name', 'FeedingBouts_EQ',       'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.FeedingBouts_UE =                  figure('Name', 'FeedingBouts_UE',       'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.Time2Feeding =                     figure('Name', 'Time2Feeding',          'units', 'normalized', 'Position', [0.25 0.25 0.50 0.25], 'Color', 'w');
hFig.N_Feeding_TC =                     figure('Name', 'N_Feeding_TC',          'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.PropTimeFeeding =                  figure('Name', 'PropTimeFeeding',       'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.NumOfVisits =                      figure('Name', 'NumOfVisits',           'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.VisitFreq =                        figure('Name', 'VisitFreq',             'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.BoutDuration =                     figure('Name', 'BoutDuration',          'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.Grp_vs_BoutDuration =              figure('Name', 'Grp_vs_BoutDuration',   'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.Grp_vs_BoutInterval =              figure('Name', 'Grp_vs_BoutInterval',   'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.Grp_vs_LeaveInterval =             figure('Name', 'Grp_vs_LeaveInterval',  'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.ConsensusHeatmaps =                figure('Name', 'ConsensusHeatmaps',     'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.SurvivalCurves =                   figure('Name', 'SurvivalCurves',        'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
%--------------------------------------------------------------------------





% Plot if set by user
if exist('hFig','var')
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'ExampleDensity')
        figure(hFig.ExampleDensity); hold on
        col = gray(6);
        col = col(1:end-1,:);
        xvec = linspace(-45, 45, 5000);
        Ani01 = normpdf(xvec, -5, 7)/normpdf(-5, -5, 7);
        Ani02 = normpdf(xvec, 2, 7)/normpdf(2, 2, 7);
        Ani03 = normpdf(xvec, 7.5, 7)/normpdf(7.5, 7.5, 7);
        Ani04 = normpdf(xvec, -20, 7)/normpdf(-20, -20, 7);
        Ani05 = normpdf(xvec, 15, 7)/normpdf(15, 15, 7);
        % --- Ani01
        plot(xvec, Ani01, 'color', col(2,:), 'linewidth', 1)
        plot(-5, ...
            normpdf(-5, 2, 7)/normpdf(2, 2, 7)+...
            normpdf(-5, 7.5, 7)/normpdf(7.5, 7.5, 7)+...
            normpdf(-5, -20, 7)/normpdf(-20, -20, 7)+...
            normpdf(-5, 15, 7)/normpdf(15, 15, 7),...
            'o', 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', 'none')
        % --- Ani02
        plot(xvec, Ani02, 'color', col(3,:), 'linewidth', 1)
        plot(2, ...
            normpdf(2, -5, 7)/normpdf(-5, -5, 7)+...
            normpdf(2, 7.5, 7)/normpdf(7.5, 7.5, 7)+...
            normpdf(2, -20, 7)/normpdf(-20, -20, 7)+...
            normpdf(2, 15, 7)/normpdf(15, 15, 7),...
            'o', 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', 'none')
        % --- Ani03
        plot(xvec,Ani03, 'color', col(4,:), 'linewidth', 1)
        plot(7.5, ...
            normpdf(7.5, -5, 7)/normpdf(-5, -5, 7)+...
            normpdf(7.5, 2, 7)/normpdf(2, 2, 7)+...
            normpdf(7.5, -20, 7)/normpdf(-20, -20, 7)+...
            normpdf(7.5, 15, 7)/normpdf(15, 15, 7),...
            'o', 'MarkerFaceColor', col(4,:), 'MarkerEdgeColor', 'none')
        % --- Ani04
        plot(xvec,Ani04, 'color', col(1,:), 'linewidth', 1)
        plot(-20, ...
            normpdf(-20, -5, 7)/normpdf(-5, -5, 7)+...
            normpdf(-20, 2, 7)/normpdf(2, 2, 7)+...
            normpdf(-20, 7.5, 7)/normpdf(7.5, 7.5, 7)+...
            normpdf(-20, 15, 7)/normpdf(15, 15, 7),...
            'o', 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', 'none')
        % --- Ani05
        plot(xvec,Ani05, 'color', col(5,:), 'linewidth', 1)
        plot(15, ...
            normpdf(15, -5, 7)/normpdf(-5, -5, 7)+...
            normpdf(15, 2, 7)/normpdf(2, 2, 7)+...
            normpdf(15, 7.5, 7)/normpdf(7.5, 7.5, 7)+...
            normpdf(15, -20, 7)/normpdf(-20, -20, 7),...
            'o', 'MarkerFaceColor', col(5,:), 'MarkerEdgeColor', 'none')
        % --- Cosmetics
        plot(xvec, Ani01+Ani02+Ani03+Ani04+Ani05, 'k:', 'linewidth', 1)
        xlim([-30 30])
        ylim([0, 3])
        xticks(unique(sort([-30:5:30,[-5 0 2 7.5 -20 25]])))
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'AmountFeeding')
        figure(hFig.AmountFeeding); hold on
        all = [];
        for iCond = 1:length(SET.ConditionNames.Patch)
            for iGrp = 1:length(SET.ConditionNames.Group)
                
                % Pool all for fit
                all = [all;...
                    PooledDATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AmountFeeding_A;...
                    PooledDATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AmountFeeding_B];
                
                % HQ
                plot(...
                    PooledDATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AmountFeeding_A(:,1),...
                    PooledDATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AmountFeeding_A(:,2),...
                    'o',...
                    'MarkerFaceColor', SET.Color.([SET.ConditionNames.Patch{iCond}, '_2']).(SET.ConditionNames.Group{iGrp}),...
                    'MarkerEdgeColor', 'none',...
                    'MarkerSize', 5)
                
                % HQ2 or LQ
                plot(...
                    PooledDATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AmountFeeding_B(:,1),...
                    PooledDATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AmountFeeding_B(:,2),...
                    'o',...
                    'MarkerFaceColor', SET.Color.([SET.ConditionNames.Patch{iCond}, '_1']).(SET.ConditionNames.Group{iGrp}),...
                    'MarkerEdgeColor', 'none',...
                    'MarkerSize', 5)
                
            end%iGrp
        end%iCond
        % Add linear fit
        mdl = fitlm(all(:,1), all(:,2), 'linear');
        x_fit = linspace(min(all(:,1)), max(all(:,1)), 1000)';
        [y_fit,CI] = predict(mdl, x_fit, 'Alpha', 0.05, 'Simultaneous', true);
        plot(x_fit, y_fit, 'k')
        plot(x_fit, CI, 'k:')
        % Cosmetics
        xlim([-25000 max(all(:,1))+25000])
        ylim([0 4])
        ylabel('\Delta patch weight (g)')
        xlabel('cumul. patch density')
        title(['r2 = ', num2str(round(mdl.Rsquared.Ordinary,2))])
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'Speed')
        figure(hFig.Speed)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Speed;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            xlim([xvec(1) xvec(end)])
            ylim([0 4])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Speed;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            xlim([xvec(1) xvec(end)])
            ylim([0 4])
        end
        clear iSub yvec CIs xvec SubPos
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'PreferenceIndex')
        figure(hFig.PreferenceIndex)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PreferenceIndex;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PreferenceIndex;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
        end
    end
    clear CIs iSub SubPos xvec yvec
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'AffectionIndex')
        figure(hFig.AffectionIndex)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).AffectionIndex;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).AffectionIndex;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
        end
    end
    clear CIs iSub SubPos xvec yvec
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'PatchPreferenceIndex')
        figure(hFig.PatchPreferenceIndex)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PatchPreferenceIndex;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PatchPreferenceIndex;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
        end
    end
    clear CIs iSub SubPos xvec yvec
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'FeedingBouts_EQ')
        figure(hFig.FeedingBouts_EQ); hold on
        % --- EQ ---
        FeedingBouts = [];
        for iN = 1:length(SET.ConditionNames.Group)
            FeedingBouts = [...
                FeedingBouts;...
                PooledDATA.EQ.(SET.ConditionNames.Group{iN}).FeedingBouts*iN];
        end%iN
        imagesc(FeedingBouts)
        title('EQ')
        yticks([])
        xticks([1 500 1000])
        set(gca, 'XTickLabels', [0 30 60])
        % Generate colormap
        col = [];
        for iN = length(SET.ConditionNames.Group):-1:1
            col = [col;...
                SET.Color.EQ_1.(SET.ConditionNames.Group{iN})];
        end
        col = [col; 1 1 1];
        for iN = 1:length(SET.ConditionNames.Group)
            col = [col;...
                SET.Color.EQ_2.(SET.ConditionNames.Group{iN})];
        end
        colormap(col)
        clear FeedingBouts col iN
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'FeedingBouts_UE')
        figure(hFig.FeedingBouts_UE); hold on
        % --- UE ---
        FeedingBouts = [];
        for iN = 1:length(SET.ConditionNames.Group)
            FeedingBouts = [...
                FeedingBouts;...
                PooledDATA.UE.(SET.ConditionNames.Group{iN}).FeedingBouts*iN];
        end%iN
        imagesc(FeedingBouts)
        title('UE')
        yticks([])
        xticks([1 500 1000])
        set(gca, 'XTickLabels', [0 30 60])
        % Generate colormap
        col = [];
        for iN = length(SET.ConditionNames.Group):-1:1
            col = [col;...
                SET.Color.UE_1.(SET.ConditionNames.Group{iN})];
        end
        col = [col; 1 1 1];
        for iN = 1:length(SET.ConditionNames.Group)
            col = [col;...
                SET.Color.UE_2.(SET.ConditionNames.Group{iN})];
        end
        colormap(col)
        clear FeedingBouts col iN
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'Time2Feeding')
        figure(hFig.Time2Feeding); hold on
        % --- EQ ---
        SubPos = 1:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Time2Feeding(:), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 1000; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub), 0.25, properties)
            clear properties ViolinData
            
        end%iSub
        % --- UE ---
        SubPos = 2:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Time2Feeding(:, 2), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 1000; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)-0.3, 0.25, properties)
            clear properties ViolinData
            
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Time2Feeding(:, 1), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 1000; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_2.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_2.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)+0.3, 0.25, properties)
            clear properties ViolinData
        end%iSub
        % Cosmetics
        xticks([1.5 5.5 9.5 13.5 17.5])
        xlim([0 19])
        set(gca, 'XTickLabels', SET.ConditionNames.Group)
        clear iSub
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'N_Feeding_TC')
        figure(hFig.N_Feeding_TC)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).N_Feeding_TC;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([0 str2double(SET.ConditionNames.Group{iSub}(2:end))])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).N_Feeding_TC;
            yvec(:,sum(isnan(yvec))==size(yvec,1)) = 0;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([0 str2double(SET.ConditionNames.Group{iSub}(2:end))])
        end
        clear iSub yvec CIs xvec SubPos
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'PropTimeFeeding')
        figure(hFig.PropTimeFeeding); hold on
        % --- EQ ---
        SubPos = 1:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:,2), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)-0.3, 0.25, properties)
            clear properties ViolinData
            
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:,1), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)+0.3, 0.25, properties)
            clear properties ViolinData
            
        end%iSub
        % --- UE ---
        SubPos = 2:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:, 2), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)-0.3, 0.25, properties)
            clear properties ViolinData
            
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:, 1), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_2.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_2.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)+0.3, 0.25, properties)
            clear properties ViolinData
        end%iSub
        xticks([1.5 5.5 9.5 13.5 17.5])
        xlim([0 19])
        set(gca, 'XTickLabels', SET.ConditionNames.Group)
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'NumOfVisits')
        figure(hFig.NumOfVisits)
        hold on
        % --- EQ ---
        SubPos = 1:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, nanmean(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).NumOfVisits, 2), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)-0.1, 0.25, properties)
            clear properties ViolinData
        end%iSub
        % --- UE ---
        SubPos = 2:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, nanmean(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).NumOfVisits, 2), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)+0.1, 0.25, properties)
            clear properties ViolinData
            
        end%iSub
        xticks([1.5 5.5 9.5 13.5 17.5])
        ylim([0 11])
        xlim([0 19])
        set(gca, 'XTickLabels', SET.ConditionNames.Group)
    end
    clear iSub
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'VisitFreq')
        figure(hFig.VisitFreq)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Get data
            yvec_A = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).InstaVisitFreq_A;
            yvec_A(:,sum(isnan(yvec_A))==size(yvec_A,1)) = 0;
            yvec_B = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).InstaVisitFreq_B;
            yvec_B(:,sum(isnan(yvec_B))==size(yvec_B,1)) = 0;
            yvec_Both = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).InstaVisitFreq_Both;
            yvec_Both(:,sum(isnan(yvec_Both))==size(yvec_Both,1)) = 0;
            % Create Subplot
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            % Get avgs
            avg_A = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            avg_B = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            avg_Both = mean(bootstrp(SET.BootSamples, @nanmean, yvec_Both, 'Options', statset('UseParallel', true)));
            % Get CIs
            CIs_A = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            CIs_B = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            CIs_Both = bootci(SET.BootSamples, {@nanmean, yvec_Both}, 'Options', statset('UseParallel', true));
            % --- Plot avg
            plot(xvec, avg_A, 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, avg_B, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, avg_Both, 'Color', [0 0 0], 'LineWidth', 3)
            % --- Plot CIs
            plot(xvec, CIs_A(1, :), 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_A(2, :), 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_Both(1, :), 'Color', [0 0 0], 'LineWidth', 1)
            plot(xvec, CIs_Both(2, :), 'Color', [0 0 0], 'LineWidth', 1)
            % Cosmetics
            xlim([xvec(1) xvec(end)])
            %         ylim([-1 1])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Get data
            yvec_A = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).InstaVisitFreq_A;
            yvec_A(:,sum(isnan(yvec_A))==size(yvec_A,1)) = 0;
            yvec_B = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).InstaVisitFreq_B;
            yvec_B(:,sum(isnan(yvec_B))==size(yvec_B,1)) = 0;
            yvec_Both = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).InstaVisitFreq_Both;
            yvec_Both(:,sum(isnan(yvec_Both))==size(yvec_Both,1)) = 0;
            % Create Subplot
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            % Get avgs
            avg_A = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            avg_B = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            avg_Both = mean(bootstrp(SET.BootSamples, @nanmean, yvec_Both, 'Options', statset('UseParallel', true)));
            % Get CIs
            CIs_A = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            CIs_B = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            CIs_Both = bootci(SET.BootSamples, {@nanmean, yvec_Both}, 'Options', statset('UseParallel', true));
            % --- Plot avg
            plot(xvec, avg_A, 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, avg_B, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, avg_Both, 'Color', [0 0 0], 'LineWidth', 3)
            % --- Plot CIs
            plot(xvec, CIs_A(1, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_A(2, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_Both(1, :), 'Color', [0 0 0], 'LineWidth', 1)
            plot(xvec, CIs_Both(2, :), 'Color', [0 0 0], 'LineWidth', 1)
            % Cosmetics
            xlim([xvec(1) xvec(end)])
            %         ylim([-1 1])
        end
    end
    clear CIs* iSub SubPos xvec yvec*
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'BoutDuration')
        figure(hFig.BoutDuration); hold on
        minVal = [];
        maxVal = [];
        % --- EQ ---
        SubPos = 1:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, log([PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1); PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1)]), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub), 0.25, properties)
            clear properties ViolinData
            
            minVal = [minVal; min([PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1); PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1)])];
            maxVal = [maxVal; max([PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1); PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1)])];
            
        end%iSub
        % --- UE ---
        SubPos = 2:4:(length(SET.ConditionNames.Group)*4);
        for iSub = 1:length(SET.ConditionNames.Group)
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1)), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_1.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)-0.3, 0.25, properties)
            clear properties ViolinData
            
            % Bootstrap data for an better estimate of the distribution
            ViolinData = bootstrp(SET.BootSamples, @nanmean, log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1)), 'Options', statset('UseParallel', true));
            % Properties for violin plot
            properties.NumPoints = 500; % Points at which to evaluate the probability density estimate
            properties.MinVal = min(ViolinData); % Smallest possible value (e.g. errors = 0)
            properties.MaxVal = max(ViolinData); % Biggest possible value
            properties.AvgType = 'mean'; %(Set which measure should be plotted: 'median', 'mean' or 'both')
            properties.EdgeCol = SET.Color.UE_2.(SET.ConditionNames.Group{iSub}); %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
            properties.MeanCol = SET.Color.UE_2.(SET.ConditionNames.Group{iSub}); %(Set colour of nanmean. If it should be depicted)
            properties.MeanWidth = 2; %(Set width of the median's line and/or symbol)
            properties.SeparateOutliers = 0; %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
            % Plot violin
            SubFcn.violinplot_advanced(ViolinData, SubPos(iSub)+0.3, 0.25, properties)
            clear properties ViolinData
            
            minVal = [minVal; min([PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1); PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1)])];
            maxVal = [maxVal; max([PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1); PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1)])];
            
        end%iSub
        xticks([1.5 5.5 9.5 13.5 17.5])
        xlim([0 19])
        set(gca, 'XTickLabels', SET.ConditionNames.Group)
    end
    clear iSub SubPos
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'Grp_vs_BoutDuration')
        figure(hFig.Grp_vs_BoutDuration)
        maxDens_up = [];
        maxDens_side = [];
        
        % *** EQ GrpSize Density ***
        subplot(6, 6, [1 2]); hold on
        for iSub = 2:length(SET.ConditionNames.Group)
            x_A = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 3);
            x_B = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 3);
            [X,Y] = stairs(1:0.5:30, histc([x_A; x_B], 1:0.5:30));
            plot(X,Y/max(Y),'Color',SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}),'LineWidth',2)
        end
        xlim([0.75 15.25])
        ylim([0 1])
        xticks([])
        yticks([])
        box on
        %--------------------------------------------------------------------------
        % *** UE GrpSize Density ***
        for iSub = 2:length(SET.ConditionNames.Group)
            x_A = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 3);
            x_B = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 3);
            [X_A,Y_A] = stairs(1:0.5:30, histc(x_A,1:0.5:30));
            [X_B,Y_B] = stairs(1:0.5:30, histc(x_B,1:0.5:30));
            % --- B
            subplot(6, 6, [19 20]); hold on
            plot(X_B, Y_B/max(Y_B), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 2);
            % --- A
            subplot(6, 6, [22 23]); hold on
            plot(X_A, Y_A/max(Y_A), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 2);
        end
        subplot(6, 6, [19 20]); hold on
        xlim([0.75 15.25])
        ylim([0 1])
        xticks([])
        yticks([])
        box on
        subplot(6, 6, [22 23]); hold on
        xlim([0.75 15.25])
        ylim([0 1])
        xticks([])
        yticks([])
        box on
        %--------------------------------------------------------------------------
        % *** EQ StayDur Density ***
        subplot(6, 6, [9 15]); hold on
        for iSub = 2:length(SET.ConditionNames.Group)
            x_A = log(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1));
            x_B = log(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1));
            [X,Y] = stairs(linspace(min(log(minVal)), max(log(maxVal)), 25), histc([x_A; x_B], linspace(min(log(minVal)), max(log(maxVal)), 25)));
            plot(Y/max(Y),X,'Color',SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}),'LineWidth',2)
        end
        xlim([0 1])
        ylim([min(log(minVal)) max(log(maxVal))])
        xticks([])
        yticks([])
        box on
        %--------------------------------------------------------------------------
        % *** UE StayDur Density ***
        for iSub = 2:length(SET.ConditionNames.Group)
            x_A = log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1));
            x_B = log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1));
            [X_A,Y_A] = stairs(linspace(min(log(minVal)), max(log(maxVal)), 25), histc(x_A, linspace(min(log(minVal)), max(log(maxVal)), 25)));
            [X_B,Y_B] = stairs(linspace(min(log(minVal)), max(log(maxVal)), 25), histc(x_B, linspace(min(log(minVal)), max(log(maxVal)), 25)));
            % --- B
            subplot(6, 6, [27 33]); hold on
            plot(Y_B/max(Y_B), X_B, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 2);
            % --- A
            subplot(6, 6, [30 36]); hold on
            plot(Y_A/max(Y_A), X_A, 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 2);
        end
        subplot(6, 6, [27 33]); hold on
        ylim([min(log(minVal)) max(log(maxVal))])
        xticks([])
        yticks([])
        box on
        subplot(6, 6, [30 36]); hold on
        xlim([0 1])
        ylim([min(log(minVal)) max(log(maxVal))])
        xticks([])
        yticks([])
        box on
        %--------------------------------------------------------------------------
        % *** EQ GrpSize vs StayDur ***
        subplot(6, 6, [7 8 13 14]); hold on
        x = [];
        y = [];
        for iSub = 2:length(SET.ConditionNames.Group)
            x = [x; PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 3)];
            y = [y; log(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1))];
            x = [x; PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 3)];
            y = [y; log(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1))];
        end
        % Bin 2D data
        xVec = 1:0.5:30;
        yVec = linspace(min(log(minVal)), max(log(maxVal)), length(xVec));
        hist3([x, y], 'ctrs', {xVec yVec}, 'CdataMode', 'auto', 'FaceColor', 'interp', 'EdgeColor', 'interp');
        colormap(SET.ColorHeat)
        % Add fit how bout changes with group size
        bins = hist3([x, y], 'ctrs', {xVec yVec})';
        yMedian = nan(length(xVec), 1);
        for iBinx = 1:length(xVec)
            temp = [];
            for iBiny = 1:length(yVec)
                temp = [temp; ones(bins(iBiny,iBinx),1)*yVec(iBiny)];
            end
            yMedian(iBinx) = median(temp);
        end
        maxBinVal = max(max(bins))*2;
        plot3(xVec, yMedian, ones(length(xVec), 1)*maxBinVal, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')
        ft = fittype('poly1');
        w = sum(bins);
        [fit_par, gof, ~] = fit(xVec(~isnan(yMedian))', yMedian(~isnan(yMedian)), ft, 'Weight', w(~isnan(yMedian)));
        f = @(x, p1, p2) p1*x + p2;
        plot3([0 30], f([0 30], fit_par.p1, fit_par.p2), [maxBinVal maxBinVal], 'w')
        title(['r^{2} = ', num2str(round(gof.rsquare, 2)), ' | RMSE = ', num2str(round(gof.rmse, 2))])
        % Cosmetics
        view(2)
        axis square
        % axis equal
        xlim([0.75 15.25])
        ylim([min(log(minVal)) max(log(maxVal))])
        box on
        %--------------------------------------------------------------------------
        % *** UE GrpSize vs StayDur (B) ***
        subplot(6, 6, [26 26 31 32]); hold on
        x = [];
        y = [];
        for iSub = 2:length(SET.ConditionNames.Group)
            x = [x; PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 3)];
            y = [y; log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1))];
        end
        % Bin 2D data
        xVec = 1:0.5:30;
        yVec = linspace(min(log(minVal)), max(log(maxVal)), length(xVec));
        hist3([x, y], 'ctrs', {xVec yVec}, 'CdataMode', 'auto', 'FaceColor', 'interp', 'EdgeColor', 'interp');
        colormap(SET.ColorHeat)
        % Add fit how bout changes with group size
        bins = hist3([x, y], 'ctrs', {xVec yVec})';
        yMedian = nan(length(xVec), 1);
        for iBinx = 1:length(xVec)
            temp = [];
            for iBiny = 1:length(yVec)
                temp = [temp; ones(bins(iBiny,iBinx),1)*yVec(iBiny)];
            end
            yMedian(iBinx) = median(temp);
        end
        maxBinVal = max(max(bins))*2;
        plot3(xVec, yMedian, ones(length(xVec), 1)*maxBinVal, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')
        ft = fittype('poly1');
        w = sum(bins);
        [fit_par, gof, ~] = fit(xVec(~isnan(yMedian))', yMedian(~isnan(yMedian)), ft, 'Weight', w(~isnan(yMedian)));
        f = @(x, p1, p2) p1*x + p2;
        plot3([0 30], f([0 30], fit_par.p1, fit_par.p2), [maxBinVal maxBinVal], 'w')
        title(['r^{2} = ', num2str(round(gof.rsquare, 2)), ' | RMSE = ', num2str(round(gof.rmse, 2))])
        % Cosmetics
        view(2)
        axis square
        % axis equal
        xlim([0.75 15.25])
        ylim([min(log(minVal)) max(log(maxVal))])
        box on
        %--------------------------------------------------------------------------
        % *** UE GrpSize vs StayDur (A) ***
        subplot(6, 6, [28 29 34 35]); hold on
        x = [];
        y = [];
        for iSub = 2:length(SET.ConditionNames.Group)
            x = [x; PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 3)];
            y = [y; log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1))];
        end
        % Bin 2D data
        xVec = 1:0.5:30;
        yVec = linspace(min(log(minVal)), max(log(maxVal)), length(xVec));
        hist3([x, y], 'ctrs', {xVec yVec}, 'CdataMode', 'auto', 'FaceColor', 'interp', 'EdgeColor', 'interp');
        colormap(SET.ColorHeat)
        % Add fit how bout changes with group size
        bins = hist3([x, y], 'ctrs', {xVec yVec})';
        yMedian = nan(length(xVec), 1);
        for iBinx = 1:length(xVec)
            temp = [];
            for iBiny = 1:length(yVec)
                temp = [temp; ones(bins(iBiny,iBinx),1)*yVec(iBiny)];
            end
            yMedian(iBinx) = median(temp);
        end
        maxBinVal = max(max(bins))*2;
        plot3(xVec, yMedian, ones(length(xVec), 1)*maxBinVal, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'none')
        ft = fittype('poly1');
        w = sum(bins);
        [fit_par, gof, ~] = fit(xVec(~isnan(yMedian))', yMedian(~isnan(yMedian)), ft, 'Weight', w(~isnan(yMedian)));
        f = @(x, p1, p2) p1*x + p2;
        plot3([0 30], f([0 30], fit_par.p1, fit_par.p2), [maxBinVal maxBinVal], 'w')
        title(['r^{2} = ', num2str(round(gof.rsquare, 2)), ' | RMSE = ', num2str(round(gof.rmse, 2))])
        % Cosmetics
        view(2)
        axis square
        % axis equal
        xlim([0.75 15.25])
        ylim([min(log(minVal)) max(log(maxVal))])
        box on
        
    end
    clear bins f f_A f_B fit_par ft gof iBin iSub maxBinVal maxDens_side
    clear maxDens_up maxVal minVal output SubPos w x x_A x_B xi_A xi_B xVec
    clear y yVal yVec
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'Grp_vs_BoutInterval')
        figure(hFig.Grp_vs_BoutInterval)
        
        %EQ
        subplot(1, 3, 1); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.EQ.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_A;PooledDATA.EQ.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_B];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, 30));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, 30);
        for iBin = 1:29
            idx = find(bin_A == iBin);
            if isempty(idx)
                poolData_interval(1, iBin) = xVec(iBin);
            elseif length(idx) == 1
                poolData_interval(1, iBin) = xVec(iBin);
                poolData_interval(2:end, iBin) = currData(idx, 1);
            else
                % Bin
                poolData_interval(1, iBin) = xVec(iBin);
                % Avg
                poolData_interval(2, iBin) = mean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                % CI
                poolData_interval(3:4, iBin) = bootci(SET.BootSamples, {@nanmean, currData(idx, 1)}, 'Options', statset('UseParallel', true));
            end
        end%iBin
        plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 2)
        plot(poolData_interval(1, :), poolData_interval(3:4, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 1)
        clear poolData_interval
        
        ylim([0 180])
        xlim([0 18])
        
        % UE-B
        subplot(1, 3, 2); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.UE.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_B];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, 30));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, 30);
        for iBin = 1:29
            idx = find(bin_A == iBin);
            if isempty(idx)
                poolData_interval(1, iBin) = xVec(iBin);
            elseif length(idx) == 1
                poolData_interval(1, iBin) = xVec(iBin);
                poolData_interval(2:end, iBin) = currData(idx, 1);
            else
                % Bin
                poolData_interval(1, iBin) = xVec(iBin);
                % Avg
                poolData_interval(2, iBin) = mean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                % CI
                poolData_interval(3:4, iBin) = bootci(SET.BootSamples, {@nanmean, currData(idx, 1)}, 'Options', statset('UseParallel', true));
            end
        end%iGrp
        plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 2)
        plot(poolData_interval(1, :), poolData_interval(3:4, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 1)
        
        ylim([0 180])
        xlim([0 18])
        
        % UE-A
        subplot(1, 3, 3); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.UE.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_A];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, 30));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, 30);
        for iBin = 1:29
            idx = find(bin_A == iBin);
            if isempty(idx)
                poolData_interval(1, iBin) = xVec(iBin);
            elseif length(idx) == 1
                poolData_interval(1, iBin) = xVec(iBin);
                poolData_interval(2:end, iBin) = currData(idx, 1);
            else
                % Bin
                poolData_interval(1, iBin) = xVec(iBin);
                % Avg
                poolData_interval(2, iBin) = mean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                % CI
                poolData_interval(3:4, iBin) = bootci(SET.BootSamples, {@nanmean, currData(idx, 1)}, 'Options', statset('UseParallel', true));
            end
        end%iGrp
        plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 2)
        plot(poolData_interval(1, :), poolData_interval(3:4, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 1)
        
        ylim([0 180])
        xlim([0 18])
    end
    clear bin_A currData edges_A iBin iCond idx poolData_interval xVec
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'Grp_vs_LeaveInterval')
        figure(hFig.Grp_vs_LeaveInterval)
        
        %EQ
        subplot(1, 3, 1); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.EQ.(SET.ConditionNames.Group{iCond}).LeavingInterval_GrpSize_A; PooledDATA.EQ.(SET.ConditionNames.Group{iCond}).LeavingInterval_GrpSize_B];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, 30));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, 30);
        for iBin = 1:29
            idx = find(bin_A == iBin);
            if isempty(idx)
                poolData_interval(1, iBin) = xVec(iBin);
            elseif length(idx) == 1
                poolData_interval(1, iBin) = xVec(iBin);
                poolData_interval(2:end, iBin) = currData(idx, 1);
            else
                % Bin
                poolData_interval(1, iBin) = xVec(iBin);
                % Avg
                poolData_interval(2, iBin) = mean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                % CI
                poolData_interval(3:4, iBin) = bootci(SET.BootSamples, {@nanmean, currData(idx, 1)}, 'Options', statset('UseParallel', true));
            end
        end%iBin
        plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 2)
        plot(poolData_interval(1, :), poolData_interval(3:4, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 1)
        clear poolData_interval
        
        ylim([0 450])
        xlim([0 18])
        
        % UE-B
        subplot(1, 3, 2); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.UE.(SET.ConditionNames.Group{iCond}).LeavingInterval_GrpSize_B];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, 30));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, 30);
        for iBin = 1:29
            idx = find(bin_A == iBin);
            if isempty(idx)
                poolData_interval(1, iBin) = xVec(iBin);
            elseif length(idx) == 1
                poolData_interval(1, iBin) = xVec(iBin);
                poolData_interval(2:end, iBin) = currData(idx, 1);
            else
                % Bin
                poolData_interval(1, iBin) = xVec(iBin);
                % Avg
                poolData_interval(2, iBin) = mean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                % CI
                poolData_interval(3:4, iBin) = bootci(SET.BootSamples, {@nanmean, currData(idx, 1)}, 'Options', statset('UseParallel', true));
            end
        end%iGrp
        plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 2)
        plot(poolData_interval(1, :), poolData_interval(3:4, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 1)
        
        ylim([0 450])
        xlim([0 18])
        
        % UE-A
        subplot(1, 3, 3); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.UE.(SET.ConditionNames.Group{iCond}).LeavingInterval_GrpSize_A];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, 30));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, 30);
        for iBin = 1:29
            idx = find(bin_A == iBin);
            if isempty(idx)
                poolData_interval(1, iBin) = xVec(iBin);
            elseif length(idx) == 1
                poolData_interval(1, iBin) = xVec(iBin);
                poolData_interval(2:end, iBin) = currData(idx, 1);
            else
                % Bin
                poolData_interval(1, iBin) = xVec(iBin);
                % Avg
                poolData_interval(2, iBin) = mean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                % CI
                poolData_interval(3:4, iBin) = bootci(SET.BootSamples, {@nanmean, currData(idx, 1)}, 'Options', statset('UseParallel', true));
            end
        end%iGrp
        plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 2)
        plot(poolData_interval(1, :), poolData_interval(3:4, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}), 'LineWidth', 1)
        
        ylim([0 450])
        xlim([0 18])
    end
    clear bin_A currData edges_A iBin iCond idx poolData_interval xVec
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'ConsensusHeatmaps')
        figure(hFig.ConsensusHeatmaps)
        
        % Get all data
        all_EQ = [];
        all_UE = [];
        for iN = 1:length(SET.ConditionNames.Group)
            all_EQ = [all_EQ; PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting];
            all_UE = [all_UE; PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting];
        end
        
        
        % Define Grid
        xGrid = 0:0.25:18;
        yGrid = linspace(0, 1, length(xGrid));
        % Preallocation
        binned_EQ = zeros(length(yGrid), length(xGrid));
        binned_UE = zeros(length(yGrid), length(xGrid));
        % Bin data
        for iX = 2:length(xGrid)+1
            % Get corresponding indices.
            % Be aware of edges
            if iX == 2
                idx_UE = find(all_UE(:, 1)>= xGrid(iX-1) & all_UE(:, 1)<= xGrid(iX));
                idx_EQ = find(all_EQ(:, 1)>= xGrid(iX-1) & all_EQ(:, 1)<= xGrid(iX));
            elseif iX == length(xGrid)+1
                idx_UE = find(all_UE(:, 1)>xGrid(iX-1));
                idx_EQ = find(all_EQ(:, 1)>xGrid(iX-1));
            else
                idx_UE = find(all_UE(:, 1)>xGrid(iX-1) & all_UE(:, 1)<= xGrid(iX));
                idx_EQ = find(all_EQ(:, 1)>xGrid(iX-1) & all_EQ(:, 1)<= xGrid(iX));
            end
            % If data are available, compute the normalized probability density
            % --- UE
            if ~isempty(idx_UE)
                binned_UE(:, iX-1) = ksdensity(all_UE(idx_UE, 2), yGrid, 'BoundaryCorrection', 'reflection');
                binned_UE(:, iX-1) = binned_UE(:, iX-1)-min(binned_UE(:, iX-1));
                binned_UE(:, iX-1) = binned_UE(:, iX-1)/max(binned_UE(:, iX-1));
            end
            % --- EQ
            if ~isempty(idx_EQ)
                binned_EQ(:, iX-1) = ksdensity(all_EQ(idx_EQ, 2), yGrid, 'BoundaryCorrection', 'reflection');
                binned_EQ(:, iX-1) = binned_EQ(:, iX-1)-min(binned_EQ(:, iX-1));
                binned_EQ(:, iX-1) = binned_EQ(:, iX-1)/max(binned_EQ(:, iX-1));
            end
        end
        
        % --- EQ
        subplot(1, 2, 1); hold on
        imagesc(binned_EQ)
        plot([0.5, length(xGrid)+0.5], [find(yGrid == 0.5), find(yGrid == 0.5)], 'w:', 'LineWidth', 2)
        axis equal
        xlim([0.5, length(xGrid)+0.5])
        ylim([0.5, length(yGrid)+0.5])
        xticks([find(xGrid == 0), find(xGrid == 3), find(xGrid == 6), find(xGrid == 9), find(xGrid == 12), find(xGrid == 15), find(xGrid == 18)])
        xticklabels(0:3:18)
        yticks([find(yGrid == 0), find(yGrid == 0.5), find(yGrid == 1)])
        yticklabels(0:0.5:1)
        title('EQ')
        colormap(SET.ColorHeat)
        caxis([0 1])
        colorbar
        
        % --- UE
        subplot(1, 2, 2); hold on
        imagesc(binned_UE)
        plot([0.5, length(xGrid)+0.5], [find(yGrid == 0.5), find(yGrid == 0.5)], 'w:', 'LineWidth', 2)
        axis equal
        xlim([0.5, length(xGrid)+0.5])
        ylim([0.5, length(yGrid)+0.5])
        xticks([find(xGrid == 0), find(xGrid == 3), find(xGrid == 6), find(xGrid == 9), find(xGrid == 12), find(xGrid == 15), find(xGrid == 18)])
        xticklabels(0:3:18)
        yticks([find(yGrid == 0), find(yGrid == 0.5), find(yGrid == 1)])
        yticklabels(0:0.5:1)
        title('UE')
        colormap(SET.ColorHeat)
        caxis([0 1])
        colorbar
    end
    clear all_EQ all_UE bin_A binned_EQ binned_UE currData edges_A iBin iCond
    clear idx idx_EQ idx_UE iN iX poolData_interval xGrid xVec yGrid
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    if isfield(hFig, 'SurvivalCurves')
        figure(hFig.SurvivalCurves); hold on
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        SubPos = 1:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec_A = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).SurvivalCurves_A;
            yvec_B = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).SurvivalCurves_B;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg_A = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            avg_B = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            CIs_A = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            CIs_B = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            plot(xvec, avg_A, 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, avg_B, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs_A(1, :), 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_A(2, :), 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            xlim([xvec(1) 5])
            ylim([0 1])
        end
        % --- UE ---
        SubPos = 2:2:(length(SET.ConditionNames.Group)*2);
        for iSub = 1:length(SET.ConditionNames.Group)
            yvec_A = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).SurvivalCurves_A;
            yvec_B = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).SurvivalCurves_B;
            subplot(length(SET.ConditionNames.Group), 2, SubPos(iSub)); hold on
            avg_A = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            avg_B = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            CIs_A = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            CIs_B = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            plot(xvec, avg_A, 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, avg_B, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs_A(1, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_A(2, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs_B(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            xlim([xvec(1) 5])
            ylim([0 1])
        end
    end
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    %--------------------------------------------------------------------------
    % Save everything
    openFigs = fieldnames(hFig);
    for iFig = 1:length(openFigs)
        figure(hFig.(openFigs{iFig}))
        export_fig(['FIG\raw\',openFigs{iFig}], '-pdf')
    end%iFig
    close all
    
    
    
end%if Plot


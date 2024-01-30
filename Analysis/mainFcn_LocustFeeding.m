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
% This is the main analysis script
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
SET.Path2Data = '...\Data\Tracking\';
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
SET.PatchWeights = '...\Data\PatchWeights.csv';
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

% Statisitcs
SET.nBootStat = 5e6;

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

%% Collect and pool data

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
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_A =                  nan(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_B =                  nan(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_A =      cell(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingBoutDuration_B =      cell(length(uniqueAnimals), 1);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PreferenceIndex =            nan(SET.CutAfter, 1);

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

            % Determine walking and standing intervals
            % --- walk
            helper = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion(:, iAni);
            helper(logical(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)+DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)))=0;
            temp = regionprops(bwlabel(helper));
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion_walk(iAni,1) = ...
                nanmedian([temp.Area]/SET.FrameRate);
            % --- stand
            helper = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion(:, iAni)<1;
            helper(logical(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)+DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)))=0;
            temp = regionprops(bwlabel(helper));
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).intermittentMotion_stand(iAni,1) = ...
                nanmedian([temp.Area]/SET.FrameRate);


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
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).FeedingInterval(iAni,1) = nanmedian(temp2);


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

            % Determine the number of visits to either patch
            % --- Patch A
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_A(iAni, 1) = length(find(diff(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_A(:, iAni)) == 1));
            % --- Patch B
            DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).N_Visit_B(iAni, 1) = length(find(diff(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).AtPatch_B(:, iAni)) == 1));

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

        end %iAni

        % Bin 2D data
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_2D = ...
            hist3([currFile.Data.pos_x, currFile.Data.pos_y], 'ctrs', {linspace(-1, 1, SET.HeatMapGrid) linspace(-1, 1, SET.HeatMapGrid)})/(size(currFile.Data,1)/currFile.N);

        % Local density at food patches, fictive patches and animals
        % --- Patch A
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_A = sum(normpdf(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_A, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2)), 2);
        % --- Patch B
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_B = sum(normpdf(DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).Dist2Patch_B, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2)), 2);
        % --- Patch E1 and E2
        % Create fictive patches E1 and E2 in empty quarters. For this,
        % rotate real patch locations py 90deg
        M = [cosd(90) -sind(90);
            sind(90)  cosd(90)];
        Dist2Patch_E1 = nan(SET.CutAfter, length(uniqueAnimals));
        Dist2Patch_E2 = nan(SET.CutAfter, length(uniqueAnimals));
        for iAni = 1:length(uniqueAnimals)
            pos_x = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_x(:,iAni);
            pos_y = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).pos_y(:,iAni);
            helper = [pos_x, pos_y] - DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchA.Loc*M;
            Dist2Patch_E1(:, iAni) = sqrt(sum(helper'.*helper'))';
            helper = [pos_x, pos_y] - DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchB.Loc*M;
            Dist2Patch_E2(:, iAni) = sqrt(sum(helper'.*helper'))';
        end
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_E1 = sum(normpdf(Dist2Patch_E1, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2)), 2);
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_E2 = sum(normpdf(Dist2Patch_E2, 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2)), 2);
        clear helper M Dist2Patch_E1 Dist2Patch_E2 iAni pos_x pos_y

        % PreferenceIndex
        helper_A = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_A;
        helper_B = DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PatchDensity_B;
        DATA.(currFile.Condition).(currFile.N_str).(['Trial_', currFile.FileName]).PreferenceIndex(:, 1) = (helper_A-helper_B)./(helper_A+helper_B);
        clear helper_*

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
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PatchDensity = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PreferenceIndex = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PropTimeFeeding = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_B = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_E = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t1 = zeros(2,2);
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Switching_t2 = zeros(4,2);
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_A = [];
            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_B = [];
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
                    currData.intermittentMotion_walk];

                % --- intermittentMotion_stand
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_stand = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).intermittentMotion_stand;...
                    currData.intermittentMotion_stand];

                % --- FeedingInterval
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval;...
                    currData.FeedingInterval];

                % --- NonFeedingInterval
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NonFeedingInterval = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).NonFeedingInterval;...
                    currData.NonFeedingInterval];

                % --- PatchDensity
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PatchDensity = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PatchDensity;...
                    nanmean(currData.PatchDensity_A), nanmean(currData.PatchDensity_B), nanmean(currData.PatchDensity_E1), nanmean(currData.PatchDensity_E2)];

                % --- PreferenceIndex
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PreferenceIndex = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PreferenceIndex;...
                    nanmean(currData.PreferenceIndex,2)'];

                % --- PropTimeFeeding
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PropTimeFeeding = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).PropTimeFeeding;...
                    sum(sum(currData.AtPatch_A))/sum(sum([currData.AtPatch_A;currData.AtPatch_B])), sum(sum(currData.AtPatch_B))/sum(sum([currData.AtPatch_A;currData.AtPatch_B]))];

                % --- Density_TC
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_A = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_A;...
                    currData.PatchDensity_A(:)'];
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_B = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_B;...
                    currData.PatchDensity_B(:)'];
                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_E = [...
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).Density_TC_E;...
                    max([currData.PatchDensity_E1(:)';currData.PatchDensity_E2(:)'])];

                % --- FeedingInterval_GrpSize_A (col1: interval, col2: density at entry)
                helper_pop_A = sum(currData.AtPatch_A, 2);
                helper_diff = diff(helper_pop_A);
                if range(helper_diff) ~= 0
                    idx = find(helper_diff>0); idx(1) = [];
                    helper_interval = diff(idx)/SET.FrameRate;
                    helper_pop_A = helper_pop_A(idx(2:end));
                    helper_density = currData.PatchDensity_A(idx(2:end));
                    % Correct density (subract the density imposed by joining
                    % animal)
                    helper_density = helper_density - normpdf(currData.PatchA.DetectionRadius, 0, SET.InteractionRange_SD/(SET.dArena/2)) / normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                    % Get the remaining density to correct for depletion of
                    % pool of available animals
                    helper_remainingAnimal = size(currData.AtPatch_A, 2) - (helper_density + currData.PatchDensity_B(idx(2:end)));
                    % Multiply interval with remaining animals
                    helper_interval = helper_interval.*helper_remainingAnimal;
                    % Put everything together
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_A = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_A;...
                        helper_interval, helper_density];
                end
                clear helper_* idx

                % --- FeedingInterval_GrpSize_B (col1: interval, col2: density at entry)
                helper_pop_B = sum(currData.AtPatch_B, 2);
                helper_diff = diff(helper_pop_B);
                if range(helper_diff) ~= 0
                    idx = find(helper_diff>0); idx(1) = [];
                    helper_interval = diff(idx)/SET.FrameRate; %1, interval
                    helper_pop_B = helper_pop_B(idx(2:end)); %2, grp size at entry
                    helper_density = currData.PatchDensity_B(idx(2:end)); %4, density at entry
                    % Correct density (subract the density imposed by joining
                    % animal)
                    helper_density = helper_density - normpdf(currData.PatchB.DetectionRadius, 0, SET.InteractionRange_SD/(SET.dArena/2)) / normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
                    % Get the remaining density to correct for depletion of
                    % pool of available animals
                    helper_remainingAnimal = size(currData.AtPatch_B, 2) - (helper_density + currData.PatchDensity_A(idx(2:end)));
                    % Multiply interval with remaining animals
                    helper_interval = helper_interval.*helper_remainingAnimal;
                    % Put everything together
                    PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_B = [...
                        PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingInterval_GrpSize_B;...
                        helper_interval, helper_density];
                end
                clear helper_* idx

                % --- FeedingBout_GrpSize_A
                % 1: duration of stay
                % 2: number of animals present
                % 3: density at patch
                helper_pop = sum(currData.AtPatch_A, 2);
                for iAni = 1:str2double(SET.ConditionNames.Group{iGroup}(2:end))
                    helper_stay = currData.FeedingBoutDuration_A{iAni, 1};
                    if ~isempty(helper_stay)
                        for iStay = 1:size(helper_stay, 1)
                            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_A = [...
                                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_A;...
                                helper_stay(iStay, 3), nanmean(helper_pop(helper_stay(iStay, 1):helper_stay(iStay, 2))), nanmean(currData.PatchDensity_A(helper_stay(iStay, 1):helper_stay(iStay, 2)))];
                        end%iStay
                    end%if
                end%iAni
                clear helper_pop helper_stay

                % --- FeedingBout_GrpSize_B
                % 1: duration of stay
                % 2: number of animals present
                % 3: density at patch
                helper_pop = sum(currData.AtPatch_B, 2);
                for iAni = 1:str2double(SET.ConditionNames.Group{iGroup}(2:end))
                    helper_stay = currData.FeedingBoutDuration_B{iAni, 1};
                    if ~isempty(helper_stay)
                        for iStay = 1:size(helper_stay, 1)
                            PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_B = [...
                                PooledDATA.(SET.ConditionNames.Patch{iPatch}).(SET.ConditionNames.Group{iGroup}).FeedingBout_GrpSize_B;...
                                helper_stay(iStay, 3), nanmean(helper_pop(helper_stay(iStay, 1):helper_stay(iStay, 2))), nanmean(currData.PatchDensity_B(helper_stay(iStay, 1):helper_stay(iStay, 2)))];
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
                % 1: density summed across both patches
                % 2: abs(PreferenceIndex)
                helper = [currData.PatchDensity_A+currData.PatchDensity_B, abs((currData.PatchDensity_A-currData.PatchDensity_B)./(currData.PatchDensity_A+currData.PatchDensity_B))];
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
clearvars -except DATA PooledDATA SET STATS hFig


%% Plot data
clear hFig
% Start parallel pool
gcp;
%--------------------------------------------------------------------------
% Figure 01
hFig.ExampleDensity =                   figure('Name', 'ExampleDensity',                'units', 'normalized', 'Position', [0.25 0.25 0.25 0.25], 'Color', 'w');

% Figure 02
hFig.Density_TC =                       figure('Name', 'Density_TC',                    'units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');
hFig.PreferenceIndex =                  figure('Name', 'PreferenceIndex',               'units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');
hFig.ConsensusTC =                      figure('Name', 'ConsensusTC',                   'units', 'normalized', 'Position', [0 0 1 1], 'Color', 'w');

Figure 03
hFig.Grp_vs_BoutInterval =              figure('Name', 'Grp_vs_BoutInterval',           'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.Grp_vs_BoutDuration =              figure('Name', 'Grp_vs_BoutDuration',           'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.ConsensusHeatmaps =                figure('Name', 'ConsensusHeatmaps',             'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
% chosing the populated patch is plotted with the script mainFcn_BayesOpt_plot.m

% Figure 04
% model results are plotted with the script mainFcn_BayesOpt_plot.m

% S1
hFig.FeedingBouts_EQ =                  figure('Name', 'FeedingBouts_EQ',               'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.FeedingBouts_UE =                  figure('Name', 'FeedingBouts_UE',               'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.AmountFeeding =                    figure('Name', 'AmountFeeding',                 'units', 'normalized', 'Position', [0.25 0.25 0.25 0.25], 'Color', 'w');
hFig.PropTimeFeeding =                  figure('Name', 'PropTimeFeeding',               'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.EvidenceAccumulationIntervals =    figure('Name', 'EvidenceAccumulationIntervals', 'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
hFig.Grp_vs_BoutInterval_notPooled =      figure('Name', 'Grp_vs_BoutInterval_notPooled', 'units', 'normalized', 'Position', [0.25 0.25 0.50 0.50], 'Color', 'w');
%--------------------------------------------------------------------------





% Plot if set by user
if exist('hFig','var')

    %----------------------------------------------------------------------
    % FIGURE 01
    %----------------------------------------------------------------------
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
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------










    %----------------------------------------------------------------------
    % FIGURE 02
    %----------------------------------------------------------------------
    if isfield(hFig, 'Density_TC')
        figure(hFig.Density_TC)
        ylim_settings = [1 3 6 9 18];
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        for iSub = 1:length(SET.ConditionNames.Group)
            subplot(2,8,iSub); hold on
            % --- patch A
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Density_TC_A;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            % --- patch B
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Density_TC_B;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            % --- fictive patches
            yvec = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Density_TC_E;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            % --- Cosmetics
            ylim([0 ylim_settings(iSub)])
            xlim([xvec(1) xvec(end)])
            plot([5 5],[0 ylim_settings(iSub)],'k:')
            plot([20 20],[0 ylim_settings(iSub)],'k:')
            % --- Averages
            subplot(2,8,[6 7 8]); hold on
            idx = find(xvec>=5 & xvec<20);
            % --- patch A
            clear properties
            yvec_A = nanmean(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Density_TC_A(:,idx),2) / str2double(SET.ConditionNames.Group{iSub}(2:end));
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            plot([iSub-1/6 iSub-1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub-1/6 iSub-1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_A, iSub-1/6, 1/6, properties)
            % --- patch B
            clear properties
            yvec_B = nanmean(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Density_TC_B(:,idx),2)  / str2double(SET.ConditionNames.Group{iSub}(2:end));
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            plot([iSub+1/6 iSub+1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub+1/6 iSub+1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.EQ_2.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_B, iSub+1/6, 1/6, properties)
            % --- patch E
            yvec_E = nanmean(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).Density_TC_E(:,idx),2)  / str2double(SET.ConditionNames.Group{iSub}(2:end));
            avg_E = mean(bootstrp(SET.BootSamples, @nanmean, yvec_E, 'Options', statset('UseParallel', true)));
            plot([iSub-1/6 iSub+1/6] + [-1/12 1/12],[avg_E avg_E],'k:')
            ylim([0 1])

        end%iSub

        % --- UE ---
        for iSub = 1:length(SET.ConditionNames.Group)
            subplot(2,8,iSub+8); hold on
            % --- patch A
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Density_TC_A;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            % --- patch B
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Density_TC_B;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            % --- fictive patches
            yvec = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Density_TC_E;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', [0.5 0.5 0.5], 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
            % --- Cosmetics
            ylim([0 ylim_settings(iSub)])
            xlim([xvec(1) xvec(end)])
            plot([5 5],[0 ylim_settings(iSub)],'k:')
            plot([20 20],[0 ylim_settings(iSub)],'k:')
            % --- Averages
            subplot(2,8,[6 7 8]+8); hold on
            idx = find(xvec>=5 & xvec<20);
            % --- patch A
            clear properties
            yvec_A = nanmean(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Density_TC_A(:,idx),2) / str2double(SET.ConditionNames.Group{iSub}(2:end));
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            plot([iSub-1/6 iSub-1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub-1/6 iSub-1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.UE_2.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_A, iSub-1/6, 1/6, properties)
            % --- patch B
            clear properties
            yvec_B = nanmean(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Density_TC_B(:,idx),2)  / str2double(SET.ConditionNames.Group{iSub}(2:end));
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            plot([iSub+1/6 iSub+1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub+1/6 iSub+1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.UE_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_B, iSub+1/6, 1/6, properties)
            % --- patch E
            yvec_E = nanmean(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).Density_TC_E(:,idx),2)  / str2double(SET.ConditionNames.Group{iSub}(2:end));
            avg_E = mean(bootstrp(SET.BootSamples, @nanmean, yvec_E, 'Options', statset('UseParallel', true)));
            plot([iSub-1/6 iSub+1/6] + [-1/12 1/12],[avg_E avg_E],'k:')
            ylim([0 1])

        end%iSub
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if isfield(hFig, 'PreferenceIndex')
        figure(hFig.PreferenceIndex)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        for iSub = 1:length(SET.ConditionNames.Group)
            subplot(2,8,iSub); hold on
            yvec_EQ = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PreferenceIndex;
            yvec_EQ(:,sum(isnan(yvec_EQ))==size(yvec_EQ,1)) = 0;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_EQ, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_EQ}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
            % --- Averages
            subplot(2,8,[6 7 8]); hold on
            idx = find(xvec>=5 & xvec<20);
            clear properties
            yvec_EQ = nanmean(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PreferenceIndex(:,idx),2);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_EQ, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_EQ}, 'Options', statset('UseParallel', true));
            plot([iSub-1/6 iSub-1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub-1/6 iSub-1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_EQ, iSub-1/6, 1/6, properties)
            ylim([-1 1])

            % --- UE ---
            subplot(2,8,iSub); hold on
            yvec_UE = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PreferenceIndex;
            yvec_UE(:,sum(isnan(yvec_UE))==size(yvec_UE,1)) = 0;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_UE, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_UE}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0 0], 'k:')
            plot([5 5],[-1 1],'k:')
            plot([20 20],[-1 1],'k:')
            xlim([xvec(1) xvec(end)])
            ylim([-1 1])
            % --- Averages
            subplot(2,8,[6 7 8]); hold on
            idx = find(xvec>=5 & xvec<20);
            clear properties
            yvec_UE = nanmean(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PreferenceIndex(:,idx),2);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_UE, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_UE}, 'Options', statset('UseParallel', true));
            plot([iSub+1/6 iSub+1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub+1/6 iSub+1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.UE_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_UE, iSub+1/6, 1/6, properties)
            ylim([-1 1])

        end% iSub
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if isfield(hFig, 'ConsensusTC')
        figure(hFig.ConsensusTC)
        xvec = linspace(0, SET.CutAfter/60/SET.FrameRate, SET.CutAfter);
        % --- EQ ---
        for iSub = 2:length(SET.ConditionNames.Group)
            subplot(2,8,iSub), hold on
            yvec_EQ = abs(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PreferenceIndex);
            yvec_EQ(:,sum(isnan(yvec_EQ))==size(yvec_EQ,1)) = 0;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_EQ, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_EQ}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0.5 0.5], 'k:')
            plot([xvec(1), xvec(end)], [1/3 1/3], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([0 1])
            % --- Averages
            subplot(2,8,[6 7 8]); hold on
            idx = find(xvec>=5 & xvec<20);
            clear properties
            yvec_EQ = nanmean(abs(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PreferenceIndex(:,idx)),2);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_EQ, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_EQ}, 'Options', statset('UseParallel', true));
            plot([iSub-1/6 iSub-1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub-1/6 iSub-1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_EQ, iSub-1/6, 1/6, properties)
            ylim([0 1])
            xlim([0.5 5.5])

            % --- UE ---
            subplot(2,8,iSub); hold on
            yvec_UE = abs(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PreferenceIndex);
            yvec_UE(:,sum(isnan(yvec_UE))==size(yvec_UE,1)) = 0;
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_UE, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_UE}, 'Options', statset('UseParallel', true));
            plot(xvec, avg, 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 3)
            plot(xvec, CIs(1, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot(xvec, CIs(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iSub}), 'LineWidth', 1)
            plot([xvec(1), xvec(end)], [0.5 0.5], 'k:')
            plot([xvec(1), xvec(end)], [1/3 1/3], 'k:')
            xlim([xvec(1) xvec(end)])
            ylim([0 1])
            plot([5 5],[-1 1],'k:')
            plot([20 20],[-1 1],'k:')

            % --- Averages
            subplot(2,8,[6 7 8]); hold on
            idx = find(xvec>=5 & xvec<20);
            clear properties
            yvec_UE = nanmean(abs(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PreferenceIndex(:,idx)),2);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_UE, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_UE}, 'Options', statset('UseParallel', true));
            plot([iSub+1/6 iSub+1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub+1/6 iSub+1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.UE_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_UE, iSub+1/6, 1/6, properties)
            ylim([0 1])
            xlim([0.5 5.5])

        end%iSub
        subplot(2,8,[6 7 8]); hold on
        plot([xvec(1), xvec(end)], [0.5 0.5], 'k:')
        plot([xvec(1), xvec(end)], [1/3 1/3], 'k:')
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------










    %----------------------------------------------------------------------
    % FIGURE 03
    %----------------------------------------------------------------------
    if isfield(hFig, 'Grp_vs_BoutInterval')
        figure(hFig.Grp_vs_BoutInterval)
        SET_Bins = 20;
        % --- EQ ---
        subplot(1, 3, 1); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData;...
                PooledDATA.EQ.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_A;...
                PooledDATA.EQ.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_B];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        SET_Bins = ceil(1+3.322*log(size(currData,1))/2);
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, SET_Bins));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, SET_Bins);
        for iBin = 1:SET_Bins-1
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
        ylim([0 1400])
        xlim([0 18])


        % UE-B
        subplot(1, 3, 2); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.UE.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_B];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        SET_Bins = ceil(1+3.322*log(size(currData,1))/2);
        [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, SET_Bins));
        xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
        poolData_interval = nan(4, SET_Bins);
        for iBin = 1:SET_Bins-1
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
        ylim([0 1400])
        xlim([0 18])


        % UE-A
        subplot(1, 3, 3); hold on
        currData = [];
        for iCond = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [currData; PooledDATA.UE.(SET.ConditionNames.Group{iCond}).FeedingInterval_GrpSize_A];
        end%iCond
        currData(currData(:, 2)<0, 2) = 0;
        SET_Bins = ceil(1+3.322*log(size(currData,1))/2);
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

        ylim([0 1400])
        xlim([0 18])
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if isfield(hFig, 'Grp_vs_BoutDuration')
        figure(hFig.Grp_vs_BoutDuration)

        % Pool data
        x = [];
        y = [];
        for iSub = 2:length(SET.ConditionNames.Group)
            x = [x; PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 3)];
            y = [y; log(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1))];
            x = [x; PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 3)];
            y = [y; log(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1))];
        end
        all_EQ = [x,y];
        all_LQ = [PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 3), log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_B(:, 1))];
        all_HQ = [PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 3), log(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingBout_GrpSize_A(:, 1))];

        % Define Grid
        xGrid = 0 : 0.5 : ceil(max([all_EQ(:,1); all_LQ(:,1); all_HQ(:,1)]));
        yGrid = linspace(floor(min([all_EQ(:,2); all_LQ(:,1); all_HQ(:,2)])), ceil(max([all_EQ(:,2); all_LQ(:,1); all_HQ(:,2)])), 100);

        % Preallocation
        binned_EQ = zeros(length(yGrid), length(xGrid));
        binned_LQ = zeros(length(yGrid), length(xGrid));
        binned_HQ = zeros(length(yGrid), length(xGrid));
        cnt_EQ = zeros(length(xGrid),1);
        cnt_LQ = zeros(length(xGrid),1);
        cnt_HQ = zeros(length(xGrid),1);

        % Bin data
        for iX = 2:length(xGrid)+1
            % Get corresponding indices. Be aware of edges
            if iX == 2
                idx_EQ = find(all_EQ(:, 1)>= xGrid(iX-1) & all_EQ(:, 1)<= xGrid(iX));
                idx_LQ = find(all_LQ(:, 1)>= xGrid(iX-1) & all_LQ(:, 1)<= xGrid(iX));
                idx_HQ = find(all_HQ(:, 1)>= xGrid(iX-1) & all_HQ(:, 1)<= xGrid(iX));
            elseif iX == length(xGrid)+1
                idx_EQ = find(all_EQ(:, 1)>xGrid(iX-1));
                idx_LQ = find(all_LQ(:, 1)>xGrid(iX-1));
                idx_HQ = find(all_HQ(:, 1)>xGrid(iX-1));
            else
                idx_EQ = find(all_EQ(:, 1)>xGrid(iX-1) & all_EQ(:, 1)<= xGrid(iX));
                idx_LQ = find(all_LQ(:, 1)>xGrid(iX-1) & all_LQ(:, 1)<= xGrid(iX));
                idx_HQ = find(all_HQ(:, 1)>xGrid(iX-1) & all_HQ(:, 1)<= xGrid(iX));
            end

            % Count
            cnt_EQ(iX-1) = length(idx_EQ);
            cnt_LQ(iX-1) = length(idx_LQ);
            cnt_HQ(iX-1) = length(idx_HQ);

            % If data are available, compute the normalized probability density
            % --- EQ
            if ~isempty(idx_EQ)
                [a,b] = ksdensity(all_EQ(idx_EQ, 2),'BoundaryCorrection','reflection');
                a = a-min(a);
                a = a/quantile(a, 0.99);
                a(a>1) = 1;
                vq1 = interp1(b,a,yGrid);
                vq1(isnan(vq1)) = 0;
                binned_EQ(:, iX-1) = vq1;
            end%if EQ
            if ~isempty(idx_LQ)
                [a,b] = ksdensity(all_LQ(idx_LQ, 2),'BoundaryCorrection','reflection');
                a = a-min(a);
                a = a/quantile(a, 0.99);
                a(a>1) = 1;
                vq1 = interp1(b,a,yGrid);
                vq1(isnan(vq1)) = 0;
                binned_LQ(:, iX-1) = vq1;
            end%if LQ
            if ~isempty(idx_HQ)
                [a,b] = ksdensity(all_HQ(idx_HQ, 2),'BoundaryCorrection','reflection');
                a = a-min(a);
                a = a/quantile(a, 0.99);
                a(a>1) = 1;
                vq1 = interp1(b,a,yGrid);
                vq1(isnan(vq1)) = 0;
                binned_HQ(:, iX-1) = vq1;
            end%if LQ
            clear a b vq1
        end

        % ---------- Plot EQ ----------
        nexttile; hold on
        % Plot density
        [X,Y] = meshgrid(xGrid, yGrid);
        contourf(X,Y,binned_EQ, 100, 'LineColor', 'none')
        colormap(SET.ColorHeat)
        % Get maxima
        [~, idx] = max(binned_EQ);
        plot(xGrid(cnt_EQ>0), yGrid(idx(cnt_EQ>0)), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
        % Add fit
        ft = fittype('poly1');
        [fit_par, gof, ~] = fit(xGrid(cnt_EQ>0)', yGrid(idx(cnt_EQ>0))', ft, 'Weight', cnt_EQ(cnt_EQ>0));
        f = @(x, p1, p2) p1*x + p2;
        plot([xGrid(1) xGrid(end)], f([xGrid(1) xGrid(end)], fit_par.p1, fit_par.p2), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}))
        title(['r^{2} = ', num2str(round(gof.rsquare, 2)), ' | RMSE = ', num2str(round(gof.rmse, 2))])
        % Cosmetics
        axis square
        % axis equal
        xlim([xGrid(1) xGrid(end)])
        ylim([yGrid(1) yGrid(end)])
        box on
        colorbar

        % ---------- Plot LQ ----------
        nexttile; hold on
        % Plot density
        [X,Y] = meshgrid(xGrid, yGrid);
        contourf(X,Y,binned_LQ, 100, 'LineColor', 'none')
        colormap(SET.ColorHeat)
        % Get maxima
        [~, idx] = max(binned_LQ);
        plot(xGrid(cnt_LQ>0), yGrid(idx(cnt_LQ>0)), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
        % Add fit
        ft = fittype('poly1');
        [fit_par, gof, ~] = fit(xGrid(cnt_LQ>0)', yGrid(idx(cnt_LQ>0))', ft, 'Weight', cnt_LQ(cnt_LQ>0));
        f = @(x, p1, p2) p1*x + p2;
        plot([xGrid(1) xGrid(end)], f([xGrid(1) xGrid(end)], fit_par.p1, fit_par.p2), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}))
        title(['r^{2} = ', num2str(round(gof.rsquare, 2)), ' | RMSE = ', num2str(round(gof.rmse, 2))])
        % Cosmetics
        axis square
        % axis equal
        xlim([xGrid(1) xGrid(end)])
        ylim([yGrid(1) yGrid(end)])
        box on
        colorbar

        % ---------- Plot HQ ----------
        nexttile; hold on
        % Plot density
        [X,Y] = meshgrid(xGrid, yGrid);
        contourf(X,Y,binned_HQ, 100, 'LineColor', 'none')
        colormap(SET.ColorHeat)
        % Get maxima
        [~, idx] = max(binned_HQ);
        plot(xGrid(cnt_HQ>0), yGrid(idx(cnt_HQ>0)), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k')
        % Add fit
        ft = fittype('poly1');
        [fit_par, gof, ~] = fit(xGrid(cnt_HQ>0)', yGrid(idx(cnt_HQ>0))', ft, 'Weight', cnt_HQ(cnt_HQ>0));
        f = @(x, p1, p2) p1*x + p2;
        plot([xGrid(1) xGrid(end)], f([xGrid(1) xGrid(end)], fit_par.p1, fit_par.p2), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{find(strcmp(SET.ConditionNames.Group,'N10'))}))
        title(['r^{2} = ', num2str(round(gof.rsquare, 2)), ' | RMSE = ', num2str(round(gof.rmse, 2))])
        % Cosmetics
        axis square
        % axis equal
        xlim([xGrid(1) xGrid(end)])
        ylim([yGrid(1) yGrid(end)])
        box on
        colorbar
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %--------------------------------------------------------------------------
    if isfield(hFig, 'ConsensusHeatmaps')
        figure(hFig.ConsensusHeatmaps)
        % Get all data
        all_EQ = [];
        all_UE = [];
        cnt = 1;
        for iN = 2:length(SET.ConditionNames.Group)

            all_EQ = [all_EQ; PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting];
            all_UE = [all_UE; PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting];

            % Define Grid
            max_density = ceil(max([PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(:,1); PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(:,1)]));
            xGrid = linspace(0, max_density, 20);
            yGrid = linspace(0, 1, 100);
            % Preallocation
            binned_EQ = zeros(length(yGrid), length(xGrid));
            binned_UE = zeros(length(yGrid), length(xGrid));
            % Bin data
            for iX = 2:length(xGrid)+1
                % Get corresponding indices.
                % Be aware of edges
                if iX == 2
                    idx_UE = find(PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(:, 1)>= xGrid(iX-1) & PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(:, 1)<= xGrid(iX));
                    idx_EQ = find(PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(:, 1)>= xGrid(iX-1) & PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(:, 1)<= xGrid(iX));
                elseif iX == length(xGrid)+1
                    idx_UE = find(PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(:, 1)>xGrid(iX-1));
                    idx_EQ = find(PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(:, 1)>xGrid(iX-1));
                else
                    idx_UE = find(PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(:, 1)>xGrid(iX-1) & PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(:, 1)<= xGrid(iX));
                    idx_EQ = find(PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(:, 1)>xGrid(iX-1) & PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(:, 1)<= xGrid(iX));
                end
                % If data are available, compute the normalized probability density
                % --- UE
                if ~isempty(idx_UE)
                    [a,b] = ksdensity(PooledDATA.UE.(SET.ConditionNames.Group{iN}).Splitting(idx_UE, 2), 'BoundaryCorrection','reflection');
                    a = a-min(a);
                    a = a/quantile(a, 0.99);
                    a(a>1) = 1;
                    vq1 = interp1(b,a,yGrid);
                    vq1(isnan(vq1)) = 0;
                    binned_UE(:, iX-1) = vq1;
                    clear a b vq1
                end
                % --- EQ
                if ~isempty(idx_EQ)
                    [a,b] = ksdensity(PooledDATA.EQ.(SET.ConditionNames.Group{iN}).Splitting(idx_EQ, 2), 'BoundaryCorrection','reflection');
                    a = a-min(a);
                    a = a/quantile(a, 0.99);
                    a(a>1) = 1;
                    vq1 = interp1(b,a,yGrid);
                    vq1(isnan(vq1)) = 0;
                    binned_EQ(:, iX-1) = vq1;
                    clear a b vq1
                end
            end

            % --- EQ
            ax(cnt) = nexttile; hold on
            [X,Y] = meshgrid(xGrid, yGrid);
            contourf(X,Y,binned_EQ, 100, 'LineColor', 'none')
            xticks(0:3:18)
            yticks(0:0.5:1)
            title(['EQ | ', SET.ConditionNames.Group{iN}, ' | ', num2str(max_density)])
            caxis([0 1])
            colorbar
            colormap(ax(cnt), SET.ColorHeat)
            cnt = cnt+1;
            axis off


            % --- UE
            ax(cnt) = nexttile; hold on
            [X,Y] = meshgrid(xGrid, yGrid);
            contourf(X,Y,binned_UE, 100, 'LineColor', 'none')
            xticks(0:3:18)
            yticks(0:0.5:1)
            title(['UE | ', SET.ConditionNames.Group{iN}, ' | ', num2str(max_density)])
            caxis([0 1])
            colormap(ax(cnt), SET.ColorHeat)
            cnt = cnt+1;
            colorbar
            axis off
        end%iN

        % Define Grid
        max_density = ceil(max([all_UE(:,1); all_EQ(:,1)]));
        xGrid = linspace(0, max_density, 20);
        yGrid = linspace(0, 1, 100);
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
                [a,b] = ksdensity(all_UE(idx_UE, 2), 'BoundaryCorrection', 'reflection');
                a = a-min(a);
                a = a/quantile(a, 0.99);
                a(a>1) = 1;
                vq1 = interp1(b,a,yGrid);
                vq1(isnan(vq1)) = 0;
                binned_UE(:, iX-1) = vq1;
            end
            % --- EQ
            if ~isempty(idx_EQ)
                [a,b] = ksdensity(all_EQ(idx_EQ, 2), 'BoundaryCorrection', 'reflection');
                a = a-min(a);
                a = a/quantile(a, 0.99);
                a(a>1) = 1;
                vq1 = interp1(b,a,yGrid);
                vq1(isnan(vq1)) = 0;
                binned_EQ(:, iX-1) = vq1;
            end%if
        end%iX

        % --- EQ
        ax(cnt) = nexttile; hold on
        [X,Y] = meshgrid(xGrid, yGrid);
        contourf(X,Y,binned_EQ, 100, 'LineColor', 'none')
        xticks(0:3:18)
        yticks(0:0.5:1)
        title(['EQ all', ' | ', num2str(max_density)])
        caxis([0 1])
        colormap(ax(cnt), SET.ColorHeat)
        cnt = cnt+1;
        colorbar
        axis off

        % --- UE
        ax(cnt) = nexttile; hold on
        [X,Y] = meshgrid(xGrid, yGrid);
        contourf(X,Y,binned_UE, 100, 'LineColor', 'none')
        xticks(0:3:18)
        yticks(0:0.5:1)
        title(['UE all', ' | ', num2str(max_density)])
        caxis([0 1])
        colormap(ax(cnt), SET.ColorHeat)
        cnt = cnt+1;
        colorbar
        axis off

        % --- UE
        nexttile; hold on
        plot(-mean(binned_UE,2), Y(:,1), 'color', SET.Color.UE_1.N10)
        plot(-mean(binned_EQ,2), Y(:,1), 'color', SET.Color.EQ_1.N10)
        title('avg consensus')
        xlim([-0.9, 0])
        ylim([0 1])
        box on
        xlabel('avg prob. density')
        ylabel('consensus')
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------










    %----------------------------------------------------------------------
    % FIGURE 04
    %----------------------------------------------------------------------
    if isfield(hFig, 'EvidenceAccumulationIntervals')
        figure(hFig.EvidenceAccumulationIntervals)
        for iSub = 1:length(SET.ConditionNames.Group)
            subplot(2,length(SET.ConditionNames.Group), iSub); hold on
            % --- EQ personal
            properties.MarkerType =         'o';
            properties.MarkerFaceColor =    SET.Color.EQ_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize =         4;
            SubFcn.beeswarmplot_advanced(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).FeedingInterval, 0.75, 0.25, properties)
            plot([0.75 0.75],[nanmedian(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingInterval) nanmedian(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingInterval)],'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
            % --- UE personal
            properties.MarkerType =         'o';
            properties.MarkerFaceColor =    SET.Color.UE_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize =         4;
            SubFcn.beeswarmplot_advanced(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingInterval, 1.25, 0.25, properties)
            plot([1.25 1.25],[nanmedian(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingInterval) nanmedian(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).FeedingInterval)],'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
            % --- cosmetics personal
            set(gca, 'YScale', 'log')
            ylim([1e-1 1e4])
            xlim([0.5 1.5])
            ylabel('personal accumulation time (s)')
            set(gca, 'XTick', [0.75 1.25], 'XTickLabel', {'EQ', 'UE'})

            subplot(2,length(SET.ConditionNames.Group), iSub+length(SET.ConditionNames.Group)); hold on
            % --- EQ personal
            properties.MarkerType =         'o';
            properties.MarkerFaceColor =    SET.Color.EQ_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize =         4;
            SubFcn.beeswarmplot_advanced(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).intermittentMotion_stand, 0.75, 0.25, properties)
            plot([0.75 0.75],[nanmedian(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).intermittentMotion_stand) nanmedian(PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).intermittentMotion_stand)],'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
            % --- UE personal
            properties.MarkerType =         'o';
            properties.MarkerFaceColor =    SET.Color.UE_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize =         4;
            SubFcn.beeswarmplot_advanced(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).intermittentMotion_stand, 1.25, 0.25, properties)
            plot([1.25 1.25],[nanmedian(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).intermittentMotion_stand) nanmedian(PooledDATA.UE.(SET.ConditionNames.Group{iSub}).intermittentMotion_stand)],'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none')
            % --- cosmetics personal
            set(gca, 'YScale', 'log')
            ylim([1e-1 1e4])
            xlim([0.5 1.5])
            ylabel('social accumulation time (s)')
            set(gca, 'XTick', [0.75 1.25], 'XTickLabel', {'EQ', 'UE'})


        end%iSub
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------










    %----------------------------------------------------------------------
    % FIGURE SI
    %----------------------------------------------------------------------
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
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
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
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
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
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if isfield(hFig, 'PropTimeFeeding')
        figure(hFig.PropTimeFeeding); hold on
        for iSub = 1:length(SET.ConditionNames.Group)

            % --- EQ ---
            % --- Patch A
            yvec_A = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:,2);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            plot([iSub-2/6 iSub-2/6],CIs,'k', 'LineWidth', 2)
            plot([iSub-2/6 iSub-2/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.EQ_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_A, iSub-2/6, 1/6, properties)

            % --- Patch B
            yvec_B = PooledDATA.EQ.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:,1);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            plot([iSub-1/6 iSub-1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub-1/6 iSub-1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.EQ_2.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_B, iSub-1/6, 1/6, properties)


            % --- UE ---
            % --- Patch A
            yvec_A = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:,2);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_A, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_A}, 'Options', statset('UseParallel', true));
            plot([iSub+1/6 iSub+1/6],CIs,'k', 'LineWidth', 2)
            plot([iSub+1/6 iSub+1/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.UE_1.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_A, iSub+1/6, 1/6, properties)
            % --- Patch B
            yvec_B = PooledDATA.UE.(SET.ConditionNames.Group{iSub}).PropTimeFeeding(:,1);
            avg = mean(bootstrp(SET.BootSamples, @nanmean, yvec_B, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, yvec_B}, 'Options', statset('UseParallel', true));
            plot([iSub+2/6 iSub+2/6],CIs,'k', 'LineWidth', 2)
            plot([iSub+2/6 iSub+2/6] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.UE_2.(SET.ConditionNames.Group{iSub});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(yvec_B, iSub+2/6, 1/6, properties)

        end%iSub

        xticks(1:length(SET.ConditionNames.Group))
        xlim([0 6])
        set(gca, 'XTickLabels', SET.ConditionNames.Group)
    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if isfield(hFig, 'Grp_vs_BoutInterval_notPooled')
        figure(hFig.Grp_vs_BoutInterval_notPooled)

        % --- EQ ---
        subplot(1, 3, 1); hold on
        for iGrp = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = [...
                PooledDATA.EQ.(SET.ConditionNames.Group{iGrp}).FeedingInterval_GrpSize_A;...
                PooledDATA.EQ.(SET.ConditionNames.Group{iGrp}).FeedingInterval_GrpSize_B];
            % Get optimal number of bins
            SET_Bins = ceil(1+3.322*log(size(currData,1))/2);
            currData(currData(:, 2)<0, 2) = 0;
            [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, SET_Bins));
            xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
            poolData_interval = nan(4, SET_Bins);
            for iBin = 1:SET_Bins-1
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
                    poolData_interval(2, iBin) = nanmean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                end
            end%iBin
            plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.EQ_1.(SET.ConditionNames.Group{iGrp}), 'LineWidth', 2)
            clear poolData_interval
        end%iCond
        ylim([0 1600])


        % --- UE B ---
        subplot(1, 3, 2); hold on
        for iGrp = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = PooledDATA.UE.(SET.ConditionNames.Group{iGrp}).FeedingInterval_GrpSize_B;
            % Get optimal number of bins
            SET_Bins = ceil(1+3.322*log(size(currData,1))/2);
            currData(currData(:, 2)<0, 2) = 0;
            [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, SET_Bins));
            xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
            poolData_interval = nan(4, SET_Bins);
            for iBin = 1:SET_Bins-1
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
                    poolData_interval(2, iBin) = nanmean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                end
            end%iBin
            plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.UE_1.(SET.ConditionNames.Group{iGrp}), 'LineWidth', 2)
            clear poolData_interval
        end%iCond
        ylim([0 1600])


        % --- UE A ---
        subplot(1, 3, 3); hold on
        for iGrp = 2:length(SET.ConditionNames.Group)
            % Bin density data
            currData = PooledDATA.UE.(SET.ConditionNames.Group{iGrp}).FeedingInterval_GrpSize_A;
            % Get optimal number of bins
            SET_Bins = ceil(1+3.322*log(size(currData,1))/2);
            currData(currData(:, 2)<0, 2) = 0;
            [~, edges_A, bin_A] = histcounts(currData(:, 2), linspace(0, max(currData(:, 2))+1e-10, SET_Bins));
            xVec = edges_A(1:end-1)+mean(diff(edges_A))/2;
            poolData_interval = nan(4, SET_Bins);
            for iBin = 1:SET_Bins-1
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
                    poolData_interval(2, iBin) = nanmean(bootstrp(SET.BootSamples, @nanmean, currData(idx, 1), 'Options', statset('UseParallel', true)));
                end
            end%iBin
            plot(poolData_interval(1, :), poolData_interval(2, :), 'Color', SET.Color.UE_2.(SET.ConditionNames.Group{iGrp}), 'LineWidth', 2)
            clear poolData_interval
        end%iCond
        ylim([0 1600])

    end
    clearvars -except DATA PooledDATA SET STATS hFig
    %----------------------------------------------------------------------









    %----------------------------------------------------------------------
    % Save everything
    openFigs = fieldnames(hFig);
    for iFig = 1:length(openFigs)
        figure(hFig.(openFigs{iFig}))
        export_fig(['FIG\raw\',openFigs{iFig}], '-pdf')
    end%iFig
    close all

end%if Plot


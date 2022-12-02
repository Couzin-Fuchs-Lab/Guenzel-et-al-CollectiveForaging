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
% This is the main plotting script for a model that predicts locust feeding
% decisions.
%
% Version: 30-Nov-2022 (MATLAB R2022a)

% Tidy up
clear all
close all
clc
% Add paths
addpath(genpath(pwd))

% Load data
load('PooledData.mat')
RESULTS.output = load('22-Nov-2022_BayesOpt_ModelOutput.mat');
RESULTS.paramter = load('22-Nov-2022_BayesOpt_ModelParameters.mat');
SET.nBootStat = 5e6;

%% Pool results by forming the average for each animal

% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        % Get list of IDs
        IDs = RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id;
        uniqueIDs = unique(IDs, 'rows');
        % Preallocation
        % --- no leaky accmumulation
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg = nan(size(uniqueIDs,1),3);
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).boot = nan(size(uniqueIDs,1),3);
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal = nan(size(uniqueIDs,1),3);
        % --- no leaky accmumulation
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).noleaky.avg = nan(size(uniqueIDs,1),3);
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).noleaky.boot = nan(size(uniqueIDs,1),3);
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).noleaky.optimal = nan(size(uniqueIDs,1),3);
        for iAni = 1:size(uniqueIDs,1)
            % Get index position
            idx = find(sum(IDs == uniqueIDs(iAni,:),2) == size(IDs,2));
            % Pool everything
            % --- no leaky accmumulation
            RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg(iAni,:) = [...
                abs(nanmean(RESULTS.output.ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind(idx))),...
                abs(nanmean(RESULTS.output.ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc(idx))),...
                abs(nanmean(RESULTS.output.ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(idx)))];
            RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).boot(iAni,:) = [...
                abs(nanmean(RESULTS.output.ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind(idx))),...
                abs(nanmean(RESULTS.output.ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc(idx))),...
                abs(nanmean(RESULTS.output.ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(idx)))];
            RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal(iAni,:) = [...
                abs(nanmean(RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind(idx))),...
                abs(nanmean(RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc(idx))),...
                abs(nanmean(RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(idx)))];
        end%iAni

        % Get summary statistics as well
        % --- leaky accumulation
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg_summary = [...
            nanmean(RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg);...
            bootci(SET.BootSamples, {@nanmean, RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg})];
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).boot_summary = [...
            nanmean(RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).boot);...
            bootci(SET.BootSamples, {@nanmean, RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).boot})];
        RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal_summary = [...
            nanmean(RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal);...
            bootci(SET.BootSamples, {@nanmean, RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal})];
    end%iGrp
end%iCond

%% Plot results

ParameterSearch = {'avg'};%, 'boot', 'optimal'
% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        % Iterate over all parameter optimization methods
        for iMethod = 1:length(ParameterSearch)

            % Create a figure for each result
            figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on

            % Get data
            dat_ind =  RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod})(:,1);
            dat_soc =  RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod})(:,2);
            dat_both = RESULTSpooled.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod})(:,3);

            % --- IND ---
            avg = mean(bootstrp(SET.BootSamples, @nanmean, dat_ind, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, dat_ind}, 'Options', statset('UseParallel', true));
            plot([1-1/3 1-1/3], CIs, 'k', 'LineWidth', 2)
            plot([1-1/3 1-1/3] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = [0 0 1];
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(dat_ind, 1-1/3, 1/6, properties)

            % --- SOC ---
            avg = mean(bootstrp(SET.BootSamples, @nanmean, dat_soc, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, dat_soc}, 'Options', statset('UseParallel', true));
            plot([1 1], CIs, 'k', 'LineWidth', 2)
            plot([1 1] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = [1 0 0];
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(dat_soc, 1, 1/6, properties)

            % --- Both ---
            avg = mean(bootstrp(SET.BootSamples, @nanmean, dat_both, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmean, dat_both}, 'Options', statset('UseParallel', true));
            plot([1+1/3 1+1/3], CIs, 'k', 'LineWidth', 2)
            plot([1+1/3 1+1/3] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = SET.Color.([SET.ConditionNames.Patch{iCond},'_1']).(SET.ConditionNames.Group{iGrp});
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(dat_both, 1+1/3, 1/6, properties)

            % --- Cosmetics ---
            title([SET.ConditionNames.Patch{iCond}, ' - ', SET.ConditionNames.Group{iGrp}])
            ylim([0 1])
            ylabel('P correct')
            xticks([])
            xlim([0.5 1.5])
            export_fig(['FIG\raw_model\',SET.ConditionNames.Patch{iCond}, '_', SET.ConditionNames.Group{iGrp}, '_swarm_',ParameterSearch{iMethod}], '-pdf')

            % Statistics
            TestStat_onesample = @(x1, L_x1, PredetVal) abs(mean(x1) - PredetVal);
            [... ind vs both
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_both.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_both.s,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_both.TestStatDistribution,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_both.c] = SubFcn.BootstrapHypothesisTesting('two-sample', dat_ind, dat_both, SET.nBootStat, 1234);
            [... soc vs both
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_both.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_both.s,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_both.TestStatDistribution,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_both.c] = SubFcn.BootstrapHypothesisTesting('two-sample', dat_soc, dat_both, SET.nBootStat, 1234);
            [... ind vs chance
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_chance.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_chance.s,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_chance.TestStatDistribution,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).ind_vs_chance.c] = SubFcn.BootstrapHypothesisTesting('one-sample', dat_ind, 0.5, SET.nBootStat, 1234, TestStat_onesample);
            [... soc vs chance
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_chance.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_chance.s,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_chance.TestStatDistribution,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).soc_vs_chance.c] = SubFcn.BootstrapHypothesisTesting('one-sample', dat_soc, 0.5, SET.nBootStat, 1234, TestStat_onesample);
            [... both vs chance
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).both_vs_chance.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).both_vs_chance.s,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).both_vs_chance.TestStatDistribution,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(ParameterSearch{iMethod}).both_vs_chance.c] = SubFcn.BootstrapHypothesisTesting('one-sample', dat_both, 0.5, SET.nBootStat, 1234, TestStat_onesample);

        end%iMethod
    end%iGrp
end%iCond
close all
ModelResults = RESULTS;
ModelResultsPooled = RESULTSpooled;


%% Pool data across patch conditions and group sizes

% Pool all prediction
dat_ind = [];
dat_soc = [];
dat_both = [];
% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        dat_ind = [dat_ind;   RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind(:)];
        dat_soc = [dat_soc;   RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc(:)];
        dat_both = [dat_both; RESULTS.output.ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(:)];
    end%iGrp
end%iCond

% Get events of either reinforcement or balancing
idx_reinforce = find((dat_ind>0.5 & dat_soc>0.5) | (dat_ind<0.5 & dat_soc<0.5));
idx_reinforce_max = find(dat_ind>0.5 & dat_soc>0.5);
idx_reinforce_min = find(dat_ind<0.5 & dat_soc<0.5);
idx_balance =  find((dat_ind>0.5 & dat_soc<0.5) | (dat_ind<0.5 & dat_soc>0.5));

% Get overview of results
% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % --- Better than avg
    better_avg_congruent.(SET.ConditionNames.Patch{iCond}) = [...
        sum(dat_both(idx_reinforce) > mean([dat_ind(idx_reinforce), dat_soc(idx_reinforce)],2)),...
        sum(dat_both(idx_reinforce) < mean([dat_ind(idx_reinforce), dat_soc(idx_reinforce)],2)),...
        round(mean(dat_both(idx_reinforce) > mean([dat_ind(idx_reinforce), dat_soc(idx_reinforce)],2))*100,2)];
    better_avg_incongruent.(SET.ConditionNames.Patch{iCond}) = [...
        sum(dat_both(idx_balance) > mean([dat_ind(idx_balance), dat_soc(idx_balance)],2)),...
        sum(dat_both(idx_balance) < mean([dat_ind(idx_balance), dat_soc(idx_balance)],2)),...
        round(mean(dat_both(idx_balance) > mean([dat_ind(idx_balance), dat_soc(idx_balance)],2))*100,2)];
    % --- Reinforce
    reinforce.(SET.ConditionNames.Patch{iCond}) = zeros(2);
    reinforce.(SET.ConditionNames.Patch{iCond})(1,1) = sum(dat_both(idx_reinforce_max) > max([dat_ind(idx_reinforce_max),dat_soc(idx_reinforce_max)],[],2));
    reinforce.(SET.ConditionNames.Patch{iCond})(1,2) = sum(dat_both(idx_reinforce_max) < max([dat_ind(idx_reinforce_max),dat_soc(idx_reinforce_max)],[],2));
    reinforce.(SET.ConditionNames.Patch{iCond})(2,1) = sum(dat_both(idx_reinforce_min) < min([dat_ind(idx_reinforce_min),dat_soc(idx_reinforce_min)],[],2));
    reinforce.(SET.ConditionNames.Patch{iCond})(2,2) = sum(dat_both(idx_reinforce_min) > min([dat_ind(idx_reinforce_min),dat_soc(idx_reinforce_min)],[],2));
    % --- Balance
    balance.(SET.ConditionNames.Patch{iCond}) = [sum(dat_both(idx_balance)>0.5), sum(dat_both(idx_balance)<0.5), round(mean(dat_both(idx_balance)>0.5)*100,2)];
end%iCond

disp([SET.ConditionNames.Patch{iCond}, ' | bettern than avg (all): ', num2str(100*mean(dat_both(:) > mean([dat_ind(:), dat_soc(:)], 2)))])
disp([SET.ConditionNames.Patch{iCond}, ' | bettern than avg (congruent): ', num2str(100*mean(dat_both(idx_reinforce) > mean([dat_ind(idx_reinforce), dat_soc(idx_reinforce)], 2)))])
disp([SET.ConditionNames.Patch{iCond}, ' | bettern than avg (incongruent): ', num2str(100*mean(dat_both(idx_balance)   > mean([dat_ind(idx_balance),   dat_soc(idx_balance)],   2)))])
disp([SET.ConditionNames.Patch{iCond}, ' | balancing rescues performance: ', num2str(100*mean(dat_both(idx_balance) > 0.5))])
disp([SET.ConditionNames.Patch{iCond}, ' | reinforcement > max: ', num2str(100*mean(dat_both(idx_reinforce_max) > max([dat_ind(idx_reinforce_max), dat_soc(idx_reinforce_max)], [], 2)))])
disp([SET.ConditionNames.Patch{iCond}, ' | reinforcement < min: ', num2str(100*mean(dat_both(idx_reinforce_min) < min([dat_ind(idx_reinforce_min), dat_soc(idx_reinforce_min)], [], 2)))])
disp([SET.ConditionNames.Patch{iCond}, ' | reinforcement > 0.5: ', num2str(100*mean(dat_both(idx_reinforce) > 0.5))])


%% Plot how the two information classes interact

figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on
% Color-code how the integration behaves relative to the individual
% inforamtion classes
subplot(2,1,1); hold on
for i = 1:length(dat_both)
    plot(dat_ind(i), dat_soc(i), '.', 'color', [dat_both(i),0,1-dat_both(i)], 'Markersize', 10)
end
axis equal
xlim([0 1])
ylim([0 1])
xticks([])
yticks([])

% Visualize reinforcement and balancing
subplot(2,1,2); hold on
for i = 1:length(dat_both)
    if dat_ind(i)<0.5 & dat_soc(i)>0.5 & dat_both(i)>dat_ind(i) & dat_both(i)<dat_soc(i)
        if dat_both(i)>mean([dat_ind(i), dat_soc(i)])
            plot(dat_ind(i), dat_soc(i), '.', 'color', [31,119,180]/255, 'Markersize', 10)
        else
            plot(dat_ind(i), dat_soc(i), '.', 'color', [231,41,138]/255, 'Markersize', 10)
        end
    elseif dat_ind(i)>0.5 & dat_soc(i)<0.5 & dat_both(i)<dat_ind(i) & dat_both(i)>dat_soc(i)
        if dat_both(i)>mean([dat_ind(i), dat_soc(i)])
            plot(dat_ind(i), dat_soc(i), '.', 'color', [31,119,180]/255, 'Markersize', 10)
        else
            plot(dat_ind(i), dat_soc(i), '.', 'color', [231,41,138]/255, 'Markersize', 10)
        end
    elseif dat_ind(i)<0.5 & dat_soc(i)<0.5 & dat_both(i)<min([dat_ind(i), dat_soc(i)])
        plot(dat_ind(i), dat_soc(i), '.', 'color', [255,127,14]/255, 'Markersize', 10)
    elseif dat_ind(i)>0.5 & dat_soc(i)>0.5 & dat_both(i)>max([dat_ind(i), dat_soc(i)])
        plot(dat_ind(i), dat_soc(i), '.', 'color', [255,127,14]/255, 'Markersize', 10)
    end
end
axis equal
xlim([0 1])
ylim([0 1])
xticks([])
yticks([])
export_fig('FIG\raw_model\modelPerformance', '-pdf')



%% Plot distribution of patch choices animals in groups

% Pool model data
ModelData_pool.EQ = [];
ModelData_pool.UE = [];
% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 2:length(SET.ConditionNames.Group)
        ModelData_pool.(SET.ConditionNames.Patch{iCond}) = [...
            ModelData_pool.(SET.ConditionNames.Patch{iCond});...
            RESULTS.output.ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal
            ];
    end%iGrp
end%iCond

figure('units', 'normalized', 'Position', [0 0 1 1])


% -------------------------------------------------------------------------
% Get data - EQ
% -------------------------------------------------------------------------
T_EQ = [ModelData_pool.EQ.choice(:),  ModelData_pool.EQ.soc_curr_A(:), ModelData_pool.EQ.soc_curr_B(:)];
% Kick out some choices
T_EQ(find(T_EQ(:,2) == T_EQ(:,3)),:) = [];
% Get events where the animal picked either the less or more populated one.
densChoice_sum = [];
densChoice_diff = [];
densChoice_pop = [];
for iC = 1:size(T_EQ,1)
    if T_EQ(iC,1) == 1 &&  T_EQ(iC,2)>T_EQ(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_EQ(iC,2:3)), 0];
        densChoice_diff = [densChoice_diff;...
            abs(T_EQ(iC,2)-T_EQ(iC,3)), 0];
        densChoice_pop = [densChoice_pop;...
            T_EQ(iC,2), 0];

    elseif T_EQ(iC,1) == 1 &&  T_EQ(iC,2)<T_EQ(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_EQ(iC,2:3)), 1];
        densChoice_diff = [densChoice_diff;...
            abs(T_EQ(iC,2)-T_EQ(iC,3)), 1];
        densChoice_pop = [densChoice_pop;...
            T_EQ(iC,3), 1];

    elseif T_EQ(iC,1) == 0 &&  T_EQ(iC,2)>T_EQ(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_EQ(iC,2:3)), 1];
        densChoice_diff = [densChoice_diff;...
            abs(T_EQ(iC,2)-T_EQ(iC,3)), 1];
        densChoice_pop = [densChoice_pop;...
            T_EQ(iC,2), 1];

    elseif T_EQ(iC,1) == 0 &&  T_EQ(iC,2)<T_EQ(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_EQ(iC,2:3)), 0];
        densChoice_diff = [densChoice_diff;...
            abs(T_EQ(iC,2)-T_EQ(iC,3)), 0];
        densChoice_pop = [densChoice_pop;...
            T_EQ(iC,3), 0];
    end%if
end%iC

% Display
disp(['EQ | prop. more populated patch: ', num2str(100*(1-mean(densChoice_sum(:,2))))])


% -------------------------------------------------------------------------
% Plot results for binning based on the sum of densities - EQ
% -------------------------------------------------------------------------
subplot(2,6,[1 2]); hold on
[counts_pop,centers] = hist(densChoice_sum(densChoice_sum(:,2)==0), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_pop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_pop(iC)], 'FaceColor', [1 0 0], 'EdgeColor', 'none')
end
[counts_nopop,centers] = hist(densChoice_sum(densChoice_sum(:,2)==1), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_nopop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_nopop(iC)], 'FaceColor', [0 0 1], 'EdgeColor', 'none')
end
title('EQ')
xlabel('density at both patches')
ylabel('cnts')
ylim([0 250])
xlim([-0.5 19.5])

% Statistics
[... populated vs less populated patch
    STATS.ModelOutput.EQ.densChoice_sum.p,...
    STATS.ModelOutput.EQ.densChoice_sum.s,...
    STATS.ModelOutput.EQ.densChoice_sum.TestStatDistribution,...
    STATS.ModelOutput.EQ.densChoice_sum.c] = SubFcn.BootstrapHypothesisTesting('one-sample', counts_pop-counts_nopop, 0, SET.nBootStat, 1234);
    STATS.ModelOutput.EQ.densChoice_sum.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);



% -------------------------------------------------------------------------
% Plot results for binning based on the difference in densities - EQ
% -------------------------------------------------------------------------
subplot(2,6,[3 4]); hold on
[counts_pop,centers] = hist(densChoice_diff(densChoice_diff(:,2)==0), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_pop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_pop(iC)], 'FaceColor', [1 0 0], 'EdgeColor', 'none')
end
[counts_nopop,centers] = hist(densChoice_diff(densChoice_diff(:,2)==1), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_nopop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_nopop(iC)], 'FaceColor', [0 0 1], 'EdgeColor', 'none')
end
title('EQ')
xlabel('density difference')
ylabel('cnts')
ylim([0 400])
xlim([-0.5 19.5])

% Statistics
[... populated vs less populated patch
    STATS.ModelOutput.EQ.densChoice_diff.p,...
    STATS.ModelOutput.EQ.densChoice_diff.s,...
    STATS.ModelOutput.EQ.densChoice_diff.TestStatDistribution,...
    STATS.ModelOutput.EQ.densChoice_diff.c] = SubFcn.BootstrapHypothesisTesting('one-sample', counts_pop-counts_nopop, 0, SET.nBootStat, 1234);
    STATS.ModelOutput.EQ.densChoice_diff.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);


% -------------------------------------------------------------------------
% Plot results for binning based on the density at the populated patch - EQ
% -------------------------------------------------------------------------
subplot(2,6,[5 6]); hold on
[counts_pop,centers] = hist(densChoice_pop(densChoice_pop(:,2)==0), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_pop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_pop(iC)], 'FaceColor', [1 0 0], 'EdgeColor', 'none')
end
[counts_nopop,centers] = hist(densChoice_pop(densChoice_pop(:,2)==1), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_nopop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_nopop(iC)], 'FaceColor', [0 0 1], 'EdgeColor', 'none')
end
title('EQ')
xlabel('density at populated patch')
ylabel('cnts')
ylim([0 250])
xlim([-0.5 19.5])

% Statistics
[... populated vs less populated patch
    STATS.ModelOutput.EQ.densChoice_pop.p,...
    STATS.ModelOutput.EQ.densChoice_pop.s,...
    STATS.ModelOutput.EQ.densChoice_pop.TestStatDistribution,...
    STATS.ModelOutput.EQ.densChoice_pop.c] = SubFcn.BootstrapHypothesisTesting('one-sample', counts_pop-counts_nopop, 0, SET.nBootStat, 1234);
    STATS.ModelOutput.EQ.densChoice_pop.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);


% -------------------------------------------------------------------------
% Get data - UE
% -------------------------------------------------------------------------
T_UE = [ModelData_pool.UE.choice(:),  ModelData_pool.UE.soc_curr_A(:), ModelData_pool.UE.soc_curr_B(:)];
% Kick out some choices
T_UE(find(T_UE(:,2) == T_UE(:,3)),:) = [];
% Get events where the animal picked either the less or more populated one.
densChoice_sum = [];
densChoice_diff = [];
densChoice_pop = [];
for iC = 1:size(T_UE,1)
    if T_UE(iC,1) == 1 &&  T_UE(iC,2)>T_UE(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_UE(iC,2:3)), 0];
        densChoice_diff = [densChoice_diff;...
            abs(T_UE(iC,2)-T_UE(iC,3)), 0];
        densChoice_pop = [densChoice_pop;...
            T_UE(iC,2), 0];

    elseif T_UE(iC,1) == 1 &&  T_UE(iC,2)<T_UE(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_UE(iC,2:3)), 1];
        densChoice_diff = [densChoice_diff;...
            abs(T_UE(iC,2)-T_UE(iC,3)), 1];
        densChoice_pop = [densChoice_pop;...
            T_UE(iC,3), 1];

    elseif T_UE(iC,1) == 0 &&  T_UE(iC,2)>T_UE(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_UE(iC,2:3)), 1];
        densChoice_diff = [densChoice_diff;...
            abs(T_UE(iC,2)-T_UE(iC,3)), 1];
        densChoice_pop = [densChoice_pop;...
            T_UE(iC,2), 1];

    elseif T_UE(iC,1) == 0 &&  T_UE(iC,2)<T_UE(iC,3)
        densChoice_sum = [densChoice_sum;...
            sum(T_UE(iC,2:3)), 0];
        densChoice_diff = [densChoice_diff;...
            abs(T_UE(iC,2)-T_UE(iC,3)), 0];
        densChoice_pop = [densChoice_pop;...
            T_UE(iC,3), 0];
    end%if
end%iC

% Display
disp(['UE | prop. more populated patch: ', num2str(100*(1-mean(densChoice_sum(:,2))))])


% -------------------------------------------------------------------------
% Plot results for binning based on the sum of densities - UE
% -------------------------------------------------------------------------
subplot(2,6,[7 8]); hold on
[counts_pop,centers] = hist(densChoice_sum(densChoice_sum(:,2)==0), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_pop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_pop(iC)], 'FaceColor', [1 0 0], 'EdgeColor', 'none')
end
[counts_nopop,centers] = hist(densChoice_sum(densChoice_sum(:,2)==1), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_nopop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_nopop(iC)], 'FaceColor', [0 0 1], 'EdgeColor', 'none')
end
title('UE')
xlabel('density at both patches')
ylabel('cnts')
ylim([0 250])
xlim([-0.5 19.5])

% Statistics
[... populated vs less populated patch
    STATS.ModelOutput.UE.densChoice_sum.p,...
    STATS.ModelOutput.UE.densChoice_sum.s,...
    STATS.ModelOutput.UE.densChoice_sum.TestStatDistribution,...
    STATS.ModelOutput.UE.densChoice_sum.c] = SubFcn.BootstrapHypothesisTesting('one-sample', counts_pop-counts_nopop, 0, SET.nBootStat, 1234);
    STATS.ModelOutput.UE.densChoice_sum.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);


% -------------------------------------------------------------------------
% Plot results for binning based on the difference in densities - UE
% -------------------------------------------------------------------------
subplot(2,6,[9 10]); hold on
[counts_pop,centers] = hist(densChoice_diff(densChoice_diff(:,2)==0), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_pop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_pop(iC)], 'FaceColor', [1 0 0], 'EdgeColor', 'none')
end
[counts_nopop,centers] = hist(densChoice_diff(densChoice_diff(:,2)==1), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_nopop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_nopop(iC)], 'FaceColor', [0 0 1], 'EdgeColor', 'none')
end
title('UE')
xlabel('density difference')
ylabel('cnts')
ylim([0 400])
xlim([-0.5 19.5])

% Statistics
[... populated vs less populated patch
    STATS.ModelOutput.UE.densChoice_diff.p,...
    STATS.ModelOutput.UE.densChoice_diff.s,...
    STATS.ModelOutput.UE.densChoice_diff.TestStatDistribution,...
    STATS.ModelOutput.UE.densChoice_diff.c] = SubFcn.BootstrapHypothesisTesting('one-sample', counts_pop-counts_nopop, 0, SET.nBootStat, 1234);
    STATS.ModelOutput.UE.densChoice_diff.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);


% -------------------------------------------------------------------------
% Plot results for binning based on the density at the populated patch - EQ
% -------------------------------------------------------------------------
subplot(2,6,[11 12]); hold on
[counts_pop,centers] = hist(densChoice_pop(densChoice_pop(:,2)==0), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_pop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_pop(iC)], 'FaceColor', [1 0 0], 'EdgeColor', 'none')
end
[counts_nopop,centers] = hist(densChoice_pop(densChoice_pop(:,2)==1), 0:1:19);
dC = mean([diff(centers)])/2;
for iC = 1:length(counts_nopop)
    rectangle('position', [centers(iC)-dC, 0, 1, counts_nopop(iC)], 'FaceColor', [0 0 1], 'EdgeColor', 'none')
end
title('UE')
xlabel('density at populated patch')
ylabel('cnts')
ylim([0 250])
xlim([-0.5 19.5])

% Statistics
[... populated vs less populated patch
    STATS.ModelOutput.UE.densChoice_pop.p,...
    STATS.ModelOutput.UE.densChoice_pop.s,...
    STATS.ModelOutput.UE.densChoice_pop.TestStatDistribution,...
    STATS.ModelOutput.UE.densChoice_pop.c] = SubFcn.BootstrapHypothesisTesting('one-sample', counts_pop-counts_nopop, 0, SET.nBootStat, 1234);
    STATS.ModelOutput.UE.densChoice_pop.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);




export_fig('FIG\raw_model\patchSelection', '-pdf')

%%
close all
save('PooledData.mat', 'DATA', 'PooledDATA', 'ModelResults', 'ModelResultsPooled', 'STATS', 'SET', '-v7.3')
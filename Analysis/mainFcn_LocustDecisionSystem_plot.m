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
% This is the main plotting script for a model that predicts locust feeding
% decisions.
%
% Version: 15-Jan-2022 (MATLAB R2022a)

% Tidy up
clear all
close all
clc
% Add paths
addpath(genpath(pwd))
mkdir('FIG\raw_model')

% Load data
load('PooledData.mat')
load('230110_LocustDecisionSystem_FitResults.mat');
load('230110_LocustDecisionSystem_FitParamters.mat');
SET.nBootStat = 5e6;

%% Plot results

% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        % Iterate over all parameter optimization methods

            % Create a figure for each result
            figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on

            % Get data            
            dat_ind =  FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind;
            dat_soc =  FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc;
            dat_both = FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both;

            % --- IND ---
            avg = mean(bootstrp(SET.BootSamples, @nanmedian, dat_ind, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmedian, dat_ind}, 'Options', statset('UseParallel', true));
            plot([1-1/3 1-1/3], CIs, 'k', 'LineWidth', 2)
            plot([1-1/3 1-1/3] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = [0 0 1];
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(dat_ind, 1-1/3, 1/6, properties)

            % --- SOC ---
            avg = mean(bootstrp(SET.BootSamples, @nanmedian, dat_soc, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmedian, dat_soc}, 'Options', statset('UseParallel', true));
            plot([1 1], CIs, 'k', 'LineWidth', 2)
            plot([1 1] + [-1/12 1/12] ,[avg avg],'k', 'LineWidth', 2)
            properties.MarkerFaceColor = [1 0 0];
            properties.MarkerSize = 5;
            SubFcn.beeswarmplot_advanced(dat_soc, 1, 1/6, properties)

            % --- Both ---
            avg = mean(bootstrp(SET.BootSamples, @nanmedian, dat_both, 'Options', statset('UseParallel', true)));
            CIs = bootci(SET.BootSamples, {@nanmedian, dat_both}, 'Options', statset('UseParallel', true));
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
            export_fig(['FIG\raw_model\',SET.ConditionNames.Patch{iCond}, '_', SET.ConditionNames.Group{iGrp}, '_swarm'], '-pdf')

            % Statistics
            TestStat_onesample = @(x1, L_x1, PredetVal) abs(mean(x1) - PredetVal);
            [... ind vs both
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.s,...
                ~,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', dat_ind, dat_both, SET.nBootStat, 1234);
            [... soc vs both
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.s,...
                ~,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', dat_soc, dat_both, SET.nBootStat, 1234);
            % Correct for multiple comparison
            STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.p = STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.p*2;
            STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.s = -log2(STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_both.p);
            STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.p = STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.p*2;
            STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.s = -log2(STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_both.p);
            [... ind vs chance
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_chance.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_chance.s,...
                ~,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind_vs_chance.c] = SubFcn.BootstrapHypothesisTesting('one-sample', dat_ind, 0.5, SET.nBootStat, 1234, TestStat_onesample);
            [... soc vs chance
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_chance.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_chance.s,...
                ~,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc_vs_chance.c] = SubFcn.BootstrapHypothesisTesting('one-sample', dat_soc, 0.5, SET.nBootStat, 1234, TestStat_onesample);
            [... both vs chance
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both_vs_chance.p,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both_vs_chance.s,...
                ~,...
                STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both_vs_chance.c] = SubFcn.BootstrapHypothesisTesting('one-sample', dat_both, 0.5, SET.nBootStat, 1234, TestStat_onesample);
            close all
    end%iGrp
end%iCond


%% Pool data across patch conditions and group sizes

% Pool all prediction
dat_ind = [];
dat_soc = [];
dat_both = [];

% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        dat_ind = [dat_ind;   FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind(:)];
        dat_soc = [dat_soc;   FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc(:)];
        dat_both = [dat_both; FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(:)];
    end%iGrp

    % Get events of either reinforcement or balancing
    idx_reinforce = find((dat_ind>0.5 & dat_soc>0.5) | (dat_ind<0.5 & dat_soc<0.5));
    idx_reinforce_max = find(dat_ind>0.5 & dat_soc>0.5);
    idx_reinforce_min = find(dat_ind<0.5 & dat_soc<0.5);
    idx_balance =  find((dat_ind>0.5 & dat_soc<0.5) | (dat_ind<0.5 & dat_soc>0.5));

    % Get overview of results
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

    % Dispplay some statistics
    disp([SET.ConditionNames.Patch{iCond}, ' | bettern than avg (all): ',         num2str(100*mean(dat_both(:) > mean([dat_ind(:), dat_soc(:)], 2)))])
    disp([SET.ConditionNames.Patch{iCond}, ' | bettern than avg (congruent): ',   num2str(100*mean(dat_both(idx_reinforce) > mean([dat_ind(idx_reinforce), dat_soc(idx_reinforce)], 2)))])
    disp([SET.ConditionNames.Patch{iCond}, ' | bettern than avg (incongruent): ', num2str(100*mean(dat_both(idx_balance)   > mean([dat_ind(idx_balance),   dat_soc(idx_balance)],   2)))])
    disp([SET.ConditionNames.Patch{iCond}, ' | balancing rescues performance: ',  num2str(100*mean(dat_both(idx_balance) > 0.5))])
    disp([SET.ConditionNames.Patch{iCond}, ' | reinforcement > max: ',            num2str(100*mean(dat_both(idx_reinforce_max) > max([dat_ind(idx_reinforce_max), dat_soc(idx_reinforce_max)], [], 2)))])
    disp([SET.ConditionNames.Patch{iCond}, ' | reinforcement < min: ',            num2str(100*mean(dat_both(idx_reinforce_min) < min([dat_ind(idx_reinforce_min), dat_soc(idx_reinforce_min)], [], 2)))])
    disp([SET.ConditionNames.Patch{iCond}, ' | reinforcement > 0.5: ',            num2str(100*mean(dat_both(idx_reinforce) > 0.5))])

end%iCond


%% Plot how the two information classes interact

figure('units', 'centimeters', 'Position', [5 5 30 30]); hold on

% Integration rule
subplot(2,2,1); hold on
% --- Create dummy data
p = linspace(0,1,100); q = 1-p;
avg = zeros(100);
P_decision = zeros(100);
P = @(Px_q, Px_s, Py_q, Py_s) (Px_q.*Px_s) ./ (Px_q.*Px_s + Py_q.*Py_s);
for ix = 1:100
    for iy = 1:100
        avg(iy,ix) = mean([p(iy), p(ix)]);
        P_decision(iy,ix) = P(p(iy), p(ix), q(iy), q(ix));
    end
end
% ---Plot a surface
[X,Y] = meshgrid(linspace(0,1,100));
s = surf(X,Y,P_decision);
s.EdgeColor = 'none';
% --- Plot guide lines
plot3(p,p, P(p(:), p(:), q(:), q(:)), 'k', 'LineWidth',2)
p1 = linspace(0,1,100);
p2 = linspace(1,0,100);
plot3(p1,p2, P(p1(:), p2(:), 1-p1(:), 1-p2(:)), 'k', 'LineWidth',2)
plot3([0 1],[0 1],[0 0],'k')
for i = 0.25:0.25:0.75
    plot3([i i],[i i],[0 P(i, i, 1-i, 1-i)], 'k:', 'LineWidth',1)
    plot3(i,i, P(i, i, 1-i, 1-i), 'ko', 'MarkerFaceColor', 'k')
    plot3(i,i, 0, 'ko', 'MarkerFaceColor', 'k')
end
% --- Cosmetics
axis equal
set(gca, 'xtick', [0 0.25 0.5 0.75 1], 'XTickLabel', [])
set(gca, 'ytick', [0 0.25 0.5 0.75 1], 'YTickLabel', [])
set(gca, 'ztick', [0 0.25 0.5 0.75 1], 'ZTickLabel', [])
box on
grid on
xlim([0 1])
ylim([0 1])
zlim([0 1])
set(gca, 'view', [45/2 45/2])
col = [interp1([0 0.5 1],[0 0.5 1],linspace(0,1,100))', interp1([0 0.5 1],[0 0 0],linspace(0,1,100))', interp1([0 0.5 1],[1 0.5 0],linspace(0,1,100))'];
col_val = linspace(0,1,100);
colormap(col)
colorbar
hbar = colorbar;
hbar = colorbar;
hbar.TickDirection = "out";
hbar.Ticks = 0:0.25:1;
caxis([0 1])


% Color-code how the integration behaves relative to the individual
% inforamtion classes
subplot(2,2,3); hold on
for i = 1:length(dat_both)
    [~,col_idx] = min(abs(col_val-dat_both(i)));
    plot(dat_ind(i), dat_soc(i), '.', 'color', col(col_idx,:), 'Markersize', 10)
end
% Cosmetics
axis equal
xlim([0 1])
ylim([0 1])
xticks([])
yticks([])
box on

% Visualize reinforcement and balancing
subplot(2,2,4); hold on
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
% Cosmetics
axis equal
xlim([0 1])
ylim([0 1])
xticks([])
yticks([])
box on
export_fig('FIG\raw_model\modelPerformance', '-pdf')

%% Plot distribution of patch choices animals in groups

% Pool model data
ModelData_pool.EQ = [];
ModelData_pool.UE = [];
% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 2:length(SET.ConditionNames.Group)
        % Get data for current group size
        currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});
        % Pool
        ModelData_pool.(SET.ConditionNames.Patch{iCond}) = [...
            ModelData_pool.(SET.ConditionNames.Patch{iCond});...
            SubFcn_LocustDecisionSystem.InfoIntegration(currData, [...
            FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.tau_ind,...
            FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.tau_soc]);
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
    ~,...
    STATS.ModelOutput.EQ.densChoice_sum.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', counts_pop, counts_nopop, SET.nBootStat, 1234);
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
    ~,...
    STATS.ModelOutput.EQ.densChoice_diff.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', counts_pop, counts_nopop, SET.nBootStat, 1234);
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
    ~,...
    STATS.ModelOutput.EQ.densChoice_pop.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', counts_pop, counts_nopop, SET.nBootStat, 1234);
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
    ~,...
    STATS.ModelOutput.UE.densChoice_sum.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', counts_pop, counts_nopop, SET.nBootStat, 1234);
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
    ~,...
    STATS.ModelOutput.UE.densChoice_diff.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', counts_pop, counts_nopop, SET.nBootStat, 1234);
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
    ~,...
    STATS.ModelOutput.UE.densChoice_pop.c] = SubFcn.BootstrapHypothesisTesting('two-sample-pairs', counts_pop, counts_nopop, SET.nBootStat, 1234);
STATS.ModelOutput.UE.densChoice_pop.binom = SubFcn.BinomTest(counts_pop, counts_pop+counts_nopop, 0.5);
export_fig('FIG\raw_model\patchSelection', '-pdf')


%% ------------------------------------------------------------------------
% Plot results reinforcement and balacing
% -------------------------------------------------------------------------
figure('units', 'normalized', 'Position', [0.25 0.25 0.5 0.5])
% Iterate over both patch conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Preallocation
    incongruent.ind_all = 0;
    incongruent.ind_success = 0;
    incongruent.soc_all = 0;
    incongruent.soc_success = 0;
    congruent_all = 0;
    congruent_success = 0;
    % Iterate over all group sizes and group everything
    for iGrp = 2:length(SET.ConditionNames.Group)
        % Get cases of incongruent cues
        idx_ind = find(FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind<0.5 & FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc>0.5);
        idx_soc = find(FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc<0.5 & FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind>0.5);
        % Get cases of congruent cues
        idx = find(...
            (FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind<0.5 & FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc<0.5)...
            |...
            (FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc>0.5 & FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind>0.5));
        % Pool
        incongruent.ind_all = incongruent.ind_all + length(idx_ind);
        incongruent.soc_all = incongruent.soc_all + length(idx_soc);
        incongruent.ind_success = incongruent.ind_success + sum(FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(idx_ind)>0.5);
        incongruent.soc_success = incongruent.soc_success + sum(FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(idx_soc)>0.5);
        congruent_all = congruent_all + length(idx);
        congruent_success = congruent_success + sum(FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both(idx)>0.5);
    end%iGrp
    % Plot
    subplot(1,2,iCond); hold on
    rectangle('Position',[1, 0, 0.25, incongruent.ind_success/incongruent.ind_all],'FaceColor','b','EdgeColor','none')
    rectangle('Position',[1.5, 0, 0.25, incongruent.soc_success/incongruent.soc_all],'FaceColor','r','EdgeColor','none')
    rectangle('Position',[2, 0, 0.25, congruent_success/congruent_all],'FaceColor','g','EdgeColor','none')

    text(1.125, incongruent.ind_success/incongruent.ind_all, num2str([incongruent.ind_success,incongruent.ind_all]), 'HorizontalAlignment', 'center')
    text(1.625, incongruent.soc_success/incongruent.soc_all, num2str([incongruent.soc_success,incongruent.soc_all]), 'HorizontalAlignment', 'center')
    text(2.125, congruent_success/congruent_all, num2str([congruent_success,congruent_all]), 'HorizontalAlignment', 'center')

    ylabel('prop. P_{integration} > 0.5')
    title(SET.ConditionNames.Patch{iCond})
    plot([0.5 2.75],[0.5 0.5],'k')
    xlim([0.5 2.75])
    ylim([0 1])
    % Also do stats
    STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).balance_ind.binom_p = SubFcn.BinomTest(incongruent.ind_success, incongruent.ind_all, 0.5);
    STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).balance_soc.binom_p = SubFcn.BinomTest(incongruent.soc_success, incongruent.soc_all, 0.5);
    STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).reinforce.binom_p =   SubFcn.BinomTest(congruent_success, congruent_all, 0.5);
    STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).balance_ind.binom_s = floor(-log2(STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).balance_ind.binom_p));
    STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).balance_soc.binom_s = floor(-log2(STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).balance_soc.binom_p));
    STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).reinforce.binom_s =   floor(-log2(STATS.ModelOutput.(SET.ConditionNames.Patch{iCond}).reinforce.binom_p));
end%iCond
export_fig('FIG\raw_model\ReinforceBalance', '-pdf')

%%
close all
save('PooledData.mat', 'DATA', 'PooledDATA', 'STATS', 'SET', '-v7.3')
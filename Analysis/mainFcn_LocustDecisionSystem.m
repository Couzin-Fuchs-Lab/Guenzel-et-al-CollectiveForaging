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
% This is the main script to investigate the locust decision system.
%
% Version: 15-Jan-2022 (MATLAB R2022a)

% Tidy up
clc; clear all; close all
% Add path to external export function
% (https://github.com/altmany/export_fig)
addpath(genpath('altmany-export_fig'))
% Keep things down
warning('off')
% Start processin pool
pPool = gcp;
pPool.IdleTimeout = 720;

%% Settings

% Set base path =
Setting.BasePath = '...\Analysis\';

% Set whether to estimate parameters again
Settings.OptimizeParameters = true;

% Files that for fit paramters and results
Settings.FitParamters_file = [Setting.BasePath, 'yymmdd_LocustDecisionSystem_FitParamters.mat'];
Settings.FitResults_file = [Setting.BasePath,   'yymmdd_LocustDecisionSystem_FitResults.mat'];

% Optimizable parameters
% --- private
Settings.bayesfit.q_x = optimizableVariable('q_x', [1 5]);
Settings.bayesfit.q_y = optimizableVariable('q_y', [1 5]);
Settings.bayesfit.Kq = optimizableVariable('Kq', [0, 1]);
Settings.bayesfit.tau_ind = optimizableVariable('tau_ind', [1, 900]);
% --- social
Settings.bayesfit.s = optimizableVariable('s', [1 5]);
Settings.bayesfit.Ks = optimizableVariable('Ks', [0, 1]);
Settings.bayesfit.tau_soc = optimizableVariable('tau_soc', [1, 900]);

% Number of repetitions
Settings.nRep = 1500;

% Formular for AICc
AICc = @(n, SS, k) n*log(SS/n) + 2*k + ((2*k*(k+1))/(n-k-1));

% Load data
load([Setting.BasePath, 'PooledData.mat']);

%% Optimize parameters, using bayesopt function
if Settings.OptimizeParameters
    % Iterate over both conditions
    for iCond = 1:length(SET.ConditionNames.Patch)
        % Iterate over all group sizes
        for iGrp = 1:length(SET.ConditionNames.Group)
            % Preallocation
            FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar = ...
                array2table(ones(1,7), 'VariableNames', {'q_x', 'q_y', 'Kq', 'tau_ind', 's', 'Ks', 'tau_soc'});
            % Indicate progress
            disp(' ')
            disp('Current trial:')
            disp(['Cond: ', num2str(iCond), '| Grp: ', num2str(iGrp)])
            tic

            % Get data for current trial
            currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});

            % Get optimizable parameters
            % Check whehter it was the EQ condition
            if strcmp(SET.ConditionNames.Patch{iCond}, 'EQ')
                % Check whether it was just one animal
                if iGrp == 1
                    % Fit parameters to model
                    % --- prepare
                    vars_soc = [];
                    vars_ind = [Settings.bayesfit.q_x, Settings.bayesfit.Kq, Settings.bayesfit.tau_ind];
                    fun_soc = @(vars_soc) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [1, 1]), 1, 1, 1, 1, 1);
                    fun_ind = @(vars_ind) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [vars_ind.tau_ind, 1]), vars_ind.q_x, vars_ind.q_x, 1, vars_ind.Kq, 1);
                else % ------------------------------------------------
                    % Fit parameters to model
                    % --- prepare
                    vars_soc = [Settings.bayesfit.s, Settings.bayesfit.Ks, Settings.bayesfit.tau_soc];
                    vars_ind = [Settings.bayesfit.q_x, Settings.bayesfit.Kq, Settings.bayesfit.tau_ind];
                    fun_soc = @(vars_soc) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [1, vars_soc.tau_soc]), 1, 1, vars_soc.s, 1, vars_soc.Ks);
                    fun_ind = @(vars_ind) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [vars_ind.tau_ind, 1]), vars_ind.q_x, vars_ind.q_x, 1, vars_ind.Kq, 1);
                end%if iGrp==1
            else % ----------------------------------------------------
                % Check whether it was just one animal
                if iGrp == 1
                    % Fit parameters to model
                    % --- prepare
                    vars_soc = [];
                    vars_ind = [Settings.bayesfit.q_x, Settings.bayesfit.q_y, Settings.bayesfit.Kq, Settings.bayesfit.tau_ind];
                    fun_soc = @(vars_soc) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [1, 1]), 1, 1, 1, 1, 1);
                    fun_ind = @(vars_ind) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [vars_ind.tau_ind, 1]), vars_ind.q_x, vars_ind.q_y, 1, vars_ind.Kq, 1);
                else % ------------------------------------------------
                    % Fit parameters to model
                    % --- prepare
                    vars_soc = [Settings.bayesfit.s, Settings.bayesfit.Ks, Settings.bayesfit.tau_soc];
                    vars_ind = [Settings.bayesfit.q_x, Settings.bayesfit.q_y, Settings.bayesfit.Kq, Settings.bayesfit.tau_ind];
                    fun_soc = @(vars_soc) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [1, vars_soc.tau_soc]), 1, 1, vars_soc.s, 1, vars_soc.Ks);
                    fun_ind = @(vars_ind) SubFcn_LocustDecisionSystem.evalfun(SubFcn_LocustDecisionSystem.InfoIntegration(currData, [vars_ind.tau_ind, 1]), vars_ind.q_x, vars_ind.q_y, 1, vars_ind.Kq, 1);
                end%if iGrp==1
            end%if EQ

            % Fit best parameters with social information
            if ~isempty(vars_soc)
                rng(1234)
                results_soc = ...
                    bayesopt(fun_soc, vars_soc,...
                    'Verbose',0,...
                    'AcquisitionFunctionName','expected-improvement-plus',...
                    'IsObjectiveDeterministic', true,...
                    'UseParallel', true,...
                    ...'PlotFcn', [],...
                    'MinWorkerUtilization', pPool.NumWorkers-1,...
                    'MaxObjectiveEvaluations', Settings.nRep);
                % Get best parameter set
                zbest = bestPoint(results_soc);
                for iPar = 1:length(vars_soc)
                    FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.(vars_soc(iPar).Name) = zbest.(vars_soc(iPar).Name);
                end%iPar
            end

            % Fit best parameters with private information
            rng(1234)
            results_ind = ...
                bayesopt(fun_ind, vars_ind,...
                'Verbose',0,...
                'AcquisitionFunctionName','expected-improvement-plus',...
                'IsObjectiveDeterministic', true,...
                'UseParallel', true,...
                ...'PlotFcn', [],...
                'MinWorkerUtilization', pPool.NumWorkers-1,...
                'MaxObjectiveEvaluations', Settings.nRep);

            % Get best parameter set
            zbest = bestPoint(results_ind);
            for iPar = 1:length(vars_ind)
                FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.(vars_ind(iPar).Name) = zbest.(vars_ind(iPar).Name);
            end%iPar
            if iCond==1
                FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.q_y = FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.q_x;
            end
            % Save
            FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).nPar = length(vars_soc) + length(vars_ind);
            save(Settings.FitParamters_file, 'FitParamters', 'Settings'); toc
            close all
        end%iGrp
    end%iCond
    clearvars -except Settings DATA FitParamters SET AICc
else
    load(Settings.FitParamters_file)
end% if optimize

%% Get results

% Iterate over both conditions
for iCond = 1:length(SET.ConditionNames.Patch)
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        % Iterate over all trials
        TrialList = fieldnames(DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}));
        % Get parameters
        in = FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar;
        % Get data for current group size
        currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});
        % Get accumulated evidence
        data_both = SubFcn_LocustDecisionSystem.InfoIntegration(currData, [in.tau_ind, in.tau_soc]);
        % Preallocation
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = [];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = [];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = [];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = [];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.ind = [];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.soc = [];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.both = [];
        % Correct entries if group size equals one
        if iGrp == 1
            data_both.soc_leaky_A = zeros(size(data_both,1),1);
            data_both.soc_leaky_B = zeros(size(data_both,1),1);
        end
        % Also, get data sets to estimate the contribution of the
        % two inforamtion classes (personal and social)
        % --- personal info only
        data_ind = data_both;
        data_ind.soc_leaky_A = zeros(size(data_both,1),1);
        data_ind.soc_leaky_B = zeros(size(data_both,1),1);
        % --- social info only
        data_soc = data_both;
        data_soc.ind_leaky_A = zeros(size(data_both,1),1);
        data_soc.ind_leaky_B = zeros(size(data_both,1),1);
        % Get predictions
        out_both = SubFcn_LocustDecisionSystem.exefun(data_both, in.q_x, in.q_y, in.s, in.Kq, in.Ks);
        out_ind = SubFcn_LocustDecisionSystem.exefun(data_ind,   in.q_x, in.q_y, in.s, in.Kq, in.Ks);
        out_soc = SubFcn_LocustDecisionSystem.exefun(data_soc,   in.q_x, in.q_y, in.s, in.Kq, in.Ks);

        % Pool results
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = [...
            FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id;...
            ones(length(out_both),1)*iCond, ones(length(out_both),1)*iGrp, data_both.trial_ID(:)];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = [...
            FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind;...
            out_ind(:)];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = [...
            FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc;...
            out_soc(:)];
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = [...
            FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both;...
            out_both(:)];

        % Get Akqike information criterion (corrected)
        % --- Get number of parameters +1
        n = length(out_both);
        k_both = FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).nPar+1;
        k_ind = k_both-3;
        k_soc = k_both-(3+double(iCond==2));
        if iGrp == 1
            k_both = 0;
            k_ind = FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).nPar+1;
            k_soc = 0;
        end
        % --- Get sum square errors
        SS_ind = sum((ones(size(out_ind))-out_ind).^2);
        SS_soc = sum((ones(size(out_soc))-out_soc).^2);
        SS_both = sum((ones(size(out_both))-out_both).^2);
        % --- % Get AICc
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.ind = AICc(n,k_ind,SS_ind);
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.soc = AICc(n,k_soc,SS_soc);
        FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.both = AICc(n,k_both,SS_both);

    end%iGrp
end%iCond
save(Settings.FitResults_file, 'FitResults', 'Settings')
clearvars -except Settings DATA FitParamters SET AICc FitResults

%% Report mean parameters and corresponding confidence intervals
parameter_names = {'q_x', 'q_y', 'Kq', 'tau_ind', 's', 'Ks', 'tau_soc'};

% Iterate over both patch conditions
clc
for iCond = 1:length(SET.ConditionNames.Patch)
    disp(SET.ConditionNames.Patch{iCond})
    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)
        avg = nan(1,length(parameter_names));
        CI = nan(2,length(parameter_names));
        str_avg=SET.ConditionNames.Group{iGrp};
        str_CI='  ';
        for iPar = 1:length(parameter_names)            
            avg(1, iPar) = nanmean(bootstrp(5000, @nanmean, FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.(parameter_names{iPar})));
            CI(:, iPar) = bootci(5000, {@nanmean, FitParamters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).bestPar.(parameter_names{iPar})});

            str_avg = [str_avg, ' & ', num2str(round(avg(1, iPar),2))];
            str_CI = [str_CI, ' & [',num2str(round(CI(1, iPar),2)), ', ', num2str(round(CI(2, iPar),2)), ']'];
        end%iPar
        disp(str_avg)
        disp(str_CI)      
        disp('\hdashline')
    end%iGrp
end%iCond
clearvars -except Settings DATA FitParamters SET AICc FitResults

%% Report mean AICc and corresponding confidence intervals
parameter_names = {'ind', 'soc', 'both'};

% Iterate over both patch conditions
clc
for iCond = 1:length(SET.ConditionNames.Patch)
    disp(SET.ConditionNames.Patch{iCond})
    % Iterate over all group sizes

    for iPar = 1:length(parameter_names)
        avg.(parameter_names{iPar}) = nan(1,length(SET.ConditionNames.Group));
        CI.(parameter_names{iPar}) = nan(2,length(SET.ConditionNames.Group));
        str_avg.(parameter_names{iPar}) = parameter_names{iPar};
        str_CI.(parameter_names{iPar}) = '  ';
        for iGrp = 1:length(SET.ConditionNames.Group)
            try
                avg.(parameter_names{iPar})(1, iPar) = nanmean(bootstrp(5000, @nanmean, FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.(parameter_names{iPar})));
                CI.(parameter_names{iPar})(:, iPar) = bootci(5000, {@nanmean, FitResults.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).AICc.(parameter_names{iPar})});
            end
            str_avg.(parameter_names{iPar}) = [str_avg.(parameter_names{iPar}), ' & ', num2str(round(avg.(parameter_names{iPar})(1, iPar),2))];
            str_CI.(parameter_names{iPar}) =  [str_CI.(parameter_names{iPar}), ' & [',num2str(round(CI.(parameter_names{iPar})(1, iPar),2)), ', ', num2str(round(CI.(parameter_names{iPar})(2, iPar),2)), ']'];
        end%iGrp
        disp(str_avg.(parameter_names{iPar}))
        disp(str_CI.(parameter_names{iPar}))
        disp('\hdashline')
    end%iPar
end%iCond
clearvars -except Settings DATA FitParamters SET AICc FitResults




















function [ModelData, ModelParameters, ModelOutput_optimal] = mainFcn_BayesOpt_leaky
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
% This is the main analysis script for predicting locust feeding decisions.
%
% Version: 22-Nov-2022 (MATLAB R2022a)

% -------------------------------------------------------------------------
% Get data
% -------------------------------------------------------------------------
load('PooledData.mat');

% Set whether to estimate parameters again
SET.OptimizeParameters = false;

% Paths
SET.OptimizeParametersFile = '22-Nov-2022_BayesOpt_ModelParameters.mat';
SET.OutputFile = '22-Nov-2022_BayesOpt_ModelOutput.mat';

% Set how many times to shuffle
SET.BootSamples = 5000;

%% ------------------------------------------------------------------------
% Determine optimal parameters for each animal
% -------------------------------------------------------------------------
% Check whether to optimize parameters again
if SET.OptimizeParameters

    % Optimizale parameters
    % --- private
    SET.bayesfit.a_x =   optimizableVariable('a_x',[1 10]);
    SET.bayesfit.a_y =   optimizableVariable('a_y',[1 10]);
    SET.bayesfit.Ka =    optimizableVariable('Ka',[0 1]);
    SET.bayesfit.tau_ind =  optimizableVariable('tau_ind',[1e-3 900]);
    % --- social
    SET.bayesfit.s =        optimizableVariable('s',[1 10]);
    SET.bayesfit.Ks =       optimizableVariable('Ks',[0 1]);
    SET.bayesfit.tau_soc =  optimizableVariable('tau_soc',[1e-3 900]);


    % Iterate over both conditions
    for iCond = 1:length(SET.ConditionNames.Patch)

        % Iterate over all group sizes
        for iGrp = 1:length(SET.ConditionNames.Group)

            % Prespecification
            ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [];
            ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [];

            % Get the average parameter set across all trials and animals
            % Get data for current group size
            currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});

            % Get optimizable parameters
            % --- check whehter it was the EQ condition
            if iCond == 1
                % --- check whether it was just one animal
                if iGrp == 1
                    vars = [SET.bayesfit.a_x, SET.bayesfit.Ka, SET.bayesfit.tau_ind];
                    fun = @(vars) evalfun(InfoIntegration_grp(currData, [vars.tau_ind, 1]), vars.a_x, vars.a_x, 1, vars.Ka, 1);
                else % more than one animal
                    vars = [SET.bayesfit.a_x, SET.bayesfit.s, SET.bayesfit.Ka, SET.bayesfit.Ks, SET.bayesfit.tau_ind, SET.bayesfit.tau_soc];
                    fun = @(vars) evalfun(InfoIntegration_grp(currData, [vars.tau_ind, vars.tau_soc]), vars.a_x, vars.a_x, vars.s, vars.Ka, vars.Ks);
                end%if nAni
            else % UE condition
                % --- check whether it was just one animal
                if iGrp == 1
                    vars = [SET.bayesfit.a_x, SET.bayesfit.a_y, SET.bayesfit.Ka, SET.bayesfit.tau_ind];
                    fun = @(vars) evalfun(InfoIntegration_grp(currData, [vars.tau_ind, 1]), vars.a_x, vars.a_y, 1, vars.Ka, 1);
                else% more than one animal
                    vars = [SET.bayesfit.a_x, SET.bayesfit.a_y, SET.bayesfit.s, SET.bayesfit.Ka, SET.bayesfit.Ks, SET.bayesfit.tau_ind, SET.bayesfit.tau_soc];
                    fun = @(vars) evalfun(InfoIntegration_grp(currData, [vars.tau_ind, vars.tau_soc]), vars.a_x, vars.a_y, vars.s, vars.Ka, vars.Ks);
                end%if nAni
            end%if EQ

            % Determine optimal parameters for the current animal
            % using Bayesian hyperparameter optimization
            results = ...
                bayesopt(fun, vars,...
                'Verbose',0,...
                'AcquisitionFunctionName','expected-improvement-plus',...
                'IsObjectiveDeterministic', true,...
                'UseParallel', true,...
                'ExplorationRatio', 0.75,...
                'PlotFcn', [],...
                'MaxObjectiveEvaluations', 750);

            % --- check whehter it was the EQ condition
            if iCond == 1
                % --- check whether it was just one animal
                if iGrp == 1
                    ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                        ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                        results.XAtMinObjective.a_x, results.XAtMinObjective.a_x, NaN, results.XAtMinObjective.Ka, NaN, results.XAtMinObjective.tau_ind, NaN];
                else% more tha one animal
                    ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                        ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                        results.XAtMinObjective.a_x, results.XAtMinObjective.a_x, results.XAtMinObjective.s, results.XAtMinObjective.Ka, results.XAtMinObjective.Ks, results.XAtMinObjective.tau_ind, results.XAtMinObjective.tau_soc];
                end%if nAni
            else% UE condition
                % --- check whether it was just one animal
                if iGrp == 1
                    ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                        ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                        results.XAtMinObjective.a_x, results.XAtMinObjective.a_y, NaN, results.XAtMinObjective.Ka, NaN, results.XAtMinObjective.tau_ind, NaN];
                else% more tha one animal
                    ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                        ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                        results.XAtMinObjective.a_x, results.XAtMinObjective.a_y, results.XAtMinObjective.s, results.XAtMinObjective.Ka, results.XAtMinObjective.Ks, results.XAtMinObjective.tau_ind, results.XAtMinObjective.tau_soc];
                end%if nAni
            end%if EQ

            % Get a list of all trials for the current combination of condition
            % and group size
            TrialList = fieldnames(DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}));

            % Iterate over all trials
            for iTrial = 1:length(TrialList)

                % Get data for current trial
                currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(TrialList{iTrial});

                % Iterate over all animals of the current trial
                nAni = size(currData.pos_x,2);
                for iAni = 1:nAni
                    
                    % Display progress
                    clc
                    disp('Hyperparameters: ')
                    disp(['Cond: ', sprintf('%03d',iCond), ' | Group: ', sprintf('%03d',iGrp), ' | Trial: ', sprintf('%03d',iTrial), ' | Animal: ', sprintf('%03d',iAni)])

                    % Get optimizable parameters
                    % --- check whehter it was the EQ condition
                    if iCond == 1
                        % --- check whether it was just one animal
                        if nAni == 1
                            vars = [SET.bayesfit.a_x, SET.bayesfit.Ka, SET.bayesfit.tau_ind];
                            fun = @(vars) evalfun(InfoIntegration_ani(currData, iAni, iTrial, nAni, [vars.tau_ind, 1]), vars.a_x, vars.a_x, 1, vars.Ka, 1);
                        else % more than one animal
                            vars = [SET.bayesfit.a_x, SET.bayesfit.s, SET.bayesfit.Ka, SET.bayesfit.Ks, SET.bayesfit.tau_ind, SET.bayesfit.tau_soc];
                            fun = @(vars) evalfun(InfoIntegration_ani(currData, iAni, iTrial, nAni, [vars.tau_ind, vars.tau_soc]), vars.a_x, vars.a_x, vars.s, vars.Ka, vars.Ks);
                        end%if nAni
                    else % UE condition
                        % --- check whether it was just one animal
                        if nAni == 1
                            vars = [SET.bayesfit.a_x, SET.bayesfit.a_y, SET.bayesfit.Ka, SET.bayesfit.tau_ind];
                            fun = @(vars) evalfun(InfoIntegration_ani(currData, iAni, iTrial, nAni, [vars.tau_ind, 1]), vars.a_x, vars.a_y, 1, vars.Ka, 1);
                        else% more than one animal
                            vars = [SET.bayesfit.a_x, SET.bayesfit.a_y, SET.bayesfit.s, SET.bayesfit.Ka, SET.bayesfit.Ks, SET.bayesfit.tau_ind, SET.bayesfit.tau_soc];
                            fun = @(vars) evalfun(InfoIntegration_ani(currData, iAni, iTrial, nAni, [vars.tau_ind, vars.tau_soc]), vars.a_x, vars.a_y, vars.s, vars.Ka, vars.Ks);
                        end%if nAni
                    end%if EQ

                    % Determine optimal parameters for the current animal
                    % using Bayesian hyperparameter optimization
                    results = ...
                        bayesopt(fun, vars,...
                        'Verbose',0,...
                        'AcquisitionFunctionName','expected-improvement-plus',...
                        'IsObjectiveDeterministic', true,...
                        'UseParallel', true,...
                        'ExplorationRatio', 0.75,...
                        'PlotFcn', [],...
                        'MaxObjectiveEvaluations', 250);

                    % Store results
                    % --- no choice, no parameters
                    if isempty(results.XAtMinObjective)
                        ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                            ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}); nan(1,7)];
                    else % parameters were estimated
                        % --- check whehter it was the EQ condition
                        if iCond == 1
                            % --- check whether it was just one animal
                            if nAni == 1
                                ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                                    ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                                    results.XAtMinObjective.a_x, results.XAtMinObjective.a_x, NaN, results.XAtMinObjective.Ka, NaN, results.XAtMinObjective.tau_ind, NaN];
                            else% more tha one animal
                                ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                                    ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                                    results.XAtMinObjective.a_x, results.XAtMinObjective.a_x, results.XAtMinObjective.s, results.XAtMinObjective.Ka, results.XAtMinObjective.Ks, results.XAtMinObjective.tau_ind, results.XAtMinObjective.tau_soc];
                            end%if nAni
                        else% UE condition
                            % --- check whether it was just one animal
                            if nAni == 1
                                ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                                    ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                                    results.XAtMinObjective.a_x, results.XAtMinObjective.a_y, NaN, results.XAtMinObjective.Ka, NaN, results.XAtMinObjective.tau_ind, NaN];
                            else% more tha one animal
                                ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}) = [...
                                    ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});...
                                    results.XAtMinObjective.a_x, results.XAtMinObjective.a_y, results.XAtMinObjective.s, results.XAtMinObjective.Ka, results.XAtMinObjective.Ks, results.XAtMinObjective.tau_ind, results.XAtMinObjective.tau_soc];
                            end%if nAni
                        end%if EQ
                    end%if parameter

                end%iAni
            end%iTrial
        end%iGrp
    end%iCond
    % Save results
    save(SET.OptimizeParametersFile, 'ModelParameters', 'ModelParameters_avg', 'SET')
end%if optimize


%% ------------------------------------------------------------------------
% Get average prediction accutracy for each animal, based on OPTIMAL
% parameter
% -------------------------------------------------------------------------
load(SET.OptimizeParametersFile)
% Iterate over both conditions
for iCond = 1:length(SET.ConditionNames.Patch)

    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)

        % Prespecification
        ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal = [];        
        ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = [];
        ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = [];
        ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = [];
        ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = [];

        % Get a list of all trials for the current combination of condition
        % and group size
        TrialList = fieldnames(DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}));

        % Counter for all animals of the current condition x group size
        cnt = 1;

        % Iterate over all trials
        for iTrial = 1:length(TrialList)

            % Get data for current trial
            currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(TrialList{iTrial});

            % Display progress
            clc
            disp('Prediction based on best parameters')
            disp(['Cond: ', sprintf('%03d',iCond), ' | Group: ', sprintf('%03d',iGrp), ' | Trial: ', sprintf('%03d',iTrial)])

            % Iterate over all animals of the current trial
            nAni = size(currData.pos_x,2);
            for iAni = 1:nAni

                % Get parameter set
                if iGrp == 1
                    in.a_x =     ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,1);
                    in.a_y =     ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,2);
                    in.s =       1;
                    in.Ka =      ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,4);
                    in.Ks =      1;
                    in.tau_ind = ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,6);
                    in.tau_soc = 1;
                else
                    in.a_x =     ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,1);
                    in.a_y =     ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,2);
                    in.s =       ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,3);
                    in.Ka =      ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,4);
                    in.Ks =      ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,5);
                    in.tau_ind = ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,6);
                    in.tau_soc = ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(cnt,7);
                end%if single animal

                % Get data
                data_both = InfoIntegration_ani(currData, iAni, iTrial, nAni, [in.tau_ind, in.tau_soc]);
                % Correct entries if group size equals one
                if nAni == 1
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
                out_both =  exefun(data_both, in.a_x, in.a_y, in.s, in.Ka, in.Ks);
                out_ind =  exefun(data_ind, in.a_x, in.a_y, in.s, in.Ka, in.Ks);
                out_soc =  exefun(data_soc, in.a_x, in.a_y, in.s, in.Ka, in.Ks);

                % Store data
                ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal = [...
                    ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).optimal;...
                    data_both];

                % Store results
                % --- both
                ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = [...
                    ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both;...
                    out_both(:)];
                % --- ind only
                ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = [...
                    ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind;...
                    out_ind(:)];
                % --- soc only
                ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = [...
                    ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc;...
                    out_soc(:)];

                % Keep track of which animal was tested
                ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = [...
                    ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id;...
                    ones(size(out_both,1),1)*iCond, ones(size(out_both,1),1)*iGrp, ones(size(out_both,1),1)*iTrial, ones(size(out_both,1),1)*iAni];

                % Update counter
                cnt = cnt+1;

            end%iAni
        end%iTrial
    end%iGrp
end%iCond


%% ------------------------------------------------------------------------
% Get average prediction accuracy for each animal, based on AVERAGE
% parameter
% -------------------------------------------------------------------------
load(SET.OptimizeParametersFile)
% Iterate over both conditions
for iCond = 1:length(SET.ConditionNames.Patch)

    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)

        % Prespecification
        ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg = [];
        ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = [];
        ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = [];
        ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = [];
        ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = [];

        % Get a list of all trials for the current combination of condition
        % and group size
        TrialList = fieldnames(DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}));

        % Get parameter set
        if iGrp == 1
            in.a_x =     ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(1);
            in.a_y =     ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(2);
            in.s =       1;
            in.Ka =      ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(4);
            in.Ks =      1;
            in.tau_ind = ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(6);
            in.tau_soc = 1;
        else
            in.a_x =     ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(1);
            in.a_y =     ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(2);
            in.s =       ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(3);
            in.Ka =      ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(4);
            in.Ks =      ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(5);
            in.tau_ind = ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(6);
            in.tau_soc = ModelParameters_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(7);
        end%if single animal

        % Iterate over all trials
        for iTrial = 1:length(TrialList)

            % Get data for current trial
            currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(TrialList{iTrial});

            % Display progress
            clc
            disp('Prediction based on avg. parameters')
            disp(['Cond: ', sprintf('%03d',iCond), ' | Group: ', sprintf('%03d',iGrp), ' | Trial: ', sprintf('%03d',iTrial)])

            % Iterate over all animals of the current trial
            nAni = size(currData.pos_x,2);
            for iAni = 1:nAni

                % Get data
                data_both = InfoIntegration_ani(currData, iAni, iTrial, nAni, [in.tau_ind, in.tau_soc]);
                % Correct entries if group size equals one
                if nAni == 1
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
                out_both = exefun(data_both, in.a_x, in.a_y, in.s, in.Ka, in.Ks);
                out_ind =  exefun(data_ind, in.a_x, in.a_y, in.s, in.Ka, in.Ks);
                out_soc =  exefun(data_soc, in.a_x, in.a_y, in.s, in.Ka, in.Ks);

                % Store data
                ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg = [...
                    ModelData.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).avg;...
                    data_both];

                % Store results
                % --- both
                ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = [...
                    ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both;...
                    out_both(:)];
                % --- ind only
                ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = [...
                    ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind;...
                    out_ind(:)];
                % --- soc only
                ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = [...
                    ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc;...
                    out_soc(:)];

                % Keep track of which animal was tested
                ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = [...
                    ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id;...
                    ones(size(out_both,1),1)*iCond, ones(size(out_both,1),1)*iGrp, ones(size(out_both,1),1)*iTrial, ones(size(out_both,1),1)*iAni];

            end%iAni
        end%iTrial
    end%iGrp
end%iCond


%% ------------------------------------------------------------------------
% Get average prediction accuracy for each animal, based on SAMPLING
% parameter values from the distributions obtained from animal-centric
% optimal values
% -------------------------------------------------------------------------
load(SET.OptimizeParametersFile)
% Iterate over both conditions
for iCond = 1:length(SET.ConditionNames.Patch)

    % Iterate over all group sizes
    for iGrp = 1:length(SET.ConditionNames.Group)

        % Get info
        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id = ModelOutput_optimal.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id;

        % Get a list of all trials for the current combination of condition
        % and group size
        TrialList = fieldnames(DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}));

        % Display progress
        clc
        disp('Prediction based on shuffled parameters')
        disp(['Cond: ', sprintf('%03d',iCond), ' | Group: ', sprintf('%03d',iGrp)])

        % Preallocation
        boot_both = nan(size(ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both, 1), SET.BootSamples);
        boot_ind =  nan(size(ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both, 1), SET.BootSamples);
        boot_soc =  nan(size(ModelOutput_avg.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both, 1), SET.BootSamples);

        % Create distributions to draw from
        parnames = {'a_x','a_y','s','Ka','Ks','tau_ind','tau_soc'};
        for iPar = 1:length(parnames)
            dist.(parnames{iPar}) = normrnd(...
                mean(bootstrp(1e5, @nanmean, ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(:,iPar))),...
                std(bootstrp(1e5, @nanmean,  ModelParameters.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp})(:,iPar))), [1, 1e6]);
        end%iPar        

        % Run this many times each time with a different combination of
        % parameters
        for iBoot = 1:SET.BootSamples

            % Get parameter set
            parnames = {'a_x','a_y','s','Ka','Ks','tau_ind','tau_soc'};
            for iPar = 1:length(parnames)
                in.(parnames{iPar}) = randsample(dist.(parnames{iPar}), 1);
            end%iPar

            % Check whether this is the EQ trial
            if iCond == 1
                in.a_y = in.a_x;
            end

            % Check whether it is just one animal and adjust parameters
            % accordingly
            if iGrp == 1               
                in.s =       1;
                in.Ks =      1;
                in.tau_soc = 1;
            end

            % Iterate over all trials
            for iTrial = 1:length(TrialList)

                % Get data for current trial
                currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).(TrialList{iTrial});

                % Iterate over all animals of the current trial
                nAni = size(currData.pos_x,2);
                for iAni = 1:nAni

                    % Get data
                    data_both = InfoIntegration_ani(currData, iAni, iTrial, nAni, [in.tau_ind, in.tau_soc]);
                    % Correct entries if group size equals one
                    if nAni == 1
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
                    out_both =  exefun(data_both, in.a_x, in.a_y, in.s, in.Ka, in.Ks);
                    out_ind =  exefun(data_ind, in.a_x, in.a_y, in.s, in.Ka, in.Ks);
                    out_soc =  exefun(data_soc, in.a_x, in.a_y, in.s, in.Ka, in.Ks);

                    % Store result
                    % --- get index
                    idx = find(...
                        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id(:,1) == iCond  & ...
                        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id(:,2) == iGrp   & ...
                        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id(:,3) == iTrial & ...
                        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).id(:,4) == iAni);
                    % --- store
                    boot_both(idx,iBoot) = out_both;
                    boot_ind(idx,iBoot) = out_ind;
                    boot_soc(idx,iBoot) = out_soc;

                end%iAni
            end%iTrial
        end%iBoot
        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).both = nanmean(boot_both, 2);
        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).ind = nanmean(boot_ind, 2);
        ModelOutput_boot.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp}).soc = nanmean(boot_soc, 2);
    end%iGrp
end%iCond

% Save results
save(SET.OutputFile , 'ModelData', 'ModelOutput_optimal', 'ModelOutput_avg', 'ModelOutput_boot', 'ModelParameters', 'SET')

end%FCN


%% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------
function ModelData = InfoIntegration_ani(currData, iAni, iTrial, nAni, tau)
% Extract data in the correct format.

%--------------------------------------------------------------
% Setting
%--------------------------------------------------------------
SET.InteractionRange_SD = 7; %[cm]
SET.dArena = 90; %[cm]
SET.FrameRate = 25;
SET.TauInd = tau(1);
SET.TauSoc = tau(2);
SET.minSpeed = 0.25;
SET.N_bin = 7;
SET.BFnames = {'ind_num_visit'; 'ind_experience_cumul'; 'soc_density_cumul'; 'ind_experience_leaky'; 'soc_density_leaky'};
SET.ProbNames = [SET.BFnames; 'integration_cumul'; 'integration_leaky'];

%--------------------------------------------------------------
% Extract data
%--------------------------------------------------------------
% Prepare table
Tbl = nan(150000, 16);
% Counters
cnt = 1;
% Create time vector
timeVec = linspace(0, size(currData.pos_x,1)/SET.FrameRate - (1/SET.FrameRate), size(currData.pos_x,1));
dt = 1/SET.FrameRate;
% Create helper variables for the different patches
helper_A = currData.AtPatch_A(:,iAni);
helper_B = currData.AtPatch_B(:,iAni);
% Get joining events
helper_A_diff = diff(helper_A);
helper_B_diff = diff(helper_B);
join_A = find(helper_A_diff>0);
join_B = find(helper_B_diff>0);
% Check whether animal was there from the beginning
if helper_A(1) == 1
    join_A = [1; join_A];
end
if helper_B(1) == 1
    join_B = [1; join_B];
end
% Concat joining event
JoinEvents = [join_A(:) ones(length(join_A),1); join_B(:) zeros(length(join_B),1)];
% Sort joining event
[~,idx] = sort(JoinEvents(:,1));
JoinEvents = JoinEvents(idx,:);
% Get interval between joining events
if size(JoinEvents,1)==1
    JoinEvents = [JoinEvents, NaN];
elseif size(JoinEvents,1)>1
    JoinEvents = [JoinEvents, [NaN; diff(JoinEvents(:,1))]];
end
% Patch density
dens_A = currData.PatchDensity_A;
dens_A = dens_A - normpdf(currData.Dist2Patch_A(:,iAni), 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
dens_B = currData.PatchDensity_B;
dens_B = dens_B - normpdf(currData.Dist2Patch_B(:,iAni), 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
% Intermittent motion
StopBout = currData.rawSpeed(:,iAni)<SET.minSpeed;
StopBout = StopBout.*(~(currData.AtPatch_A(:,iAni)+currData.AtPatch_B(:,iAni)));
% Accumulated evidence (include memory extinction)
% PERSONAL
hist_A_cumul = cumsum(helper_A)/SET.FrameRate;
hist_B_cumul = cumsum(helper_B)/SET.FrameRate;
hist_A_leaky = SubFcn.leakyIntegrator(helper_A, SET.TauInd, timeVec, dt);
hist_B_leaky = SubFcn.leakyIntegrator(helper_B, SET.TauInd, timeVec, dt);
%SOCIAL
dens_A_cumul = cumsum(dens_A(:).*StopBout(:));
dens_B_cumul = cumsum(dens_B(:).*StopBout(:));
dens_A_leaky = SubFcn.leakyIntegrator(dens_A(:).*StopBout(:), SET.TauSoc, timeVec, dt);
dens_B_leaky = SubFcn.leakyIntegrator(dens_B(:).*StopBout(:), SET.TauSoc, timeVec, dt);

% Iterate over consecutive joinings and combine
% everything
for iEvent = 1:size(JoinEvents,1)
    % Exclude animals that were feeding from the
    % beginning
    if JoinEvents(iEvent,1)~=1
        % Keep track of actual choice, animal and trial
        % identity
        Tbl(cnt, 1) = JoinEvents(iEvent,2);
        Tbl(cnt, 2) = iAni;
        Tbl(cnt, 3) = iTrial;
        Tbl(cnt, 4) = nAni;
        % -----
        % Get the number of past visits to each shelter
        if cnt == 1
            Tbl(cnt, 5:6) = [0, 0];
        else
            if JoinEvents(iEvent-1,2) == 1
                Tbl(cnt, 5:6) = Tbl(cnt-1,5:6) + [1,0];
            elseif JoinEvents(iEvent-1,2) == 0
                Tbl(cnt, 5:6) = Tbl(cnt-1,5:6) + [0,1];
            end%if
        end%if
        % -----
        % Get the cumulative time an animal has spent
        Tbl(cnt, 7:8) = [hist_A_cumul(JoinEvents(iEvent,1)), hist_B_cumul(JoinEvents(iEvent,1))]; %ind_cumul
        % -----
        % Get the prevailing animal density
        Tbl(cnt, 9:10) = [dens_A_cumul(JoinEvents(iEvent,1)), dens_B_cumul(JoinEvents(iEvent,1))]; %soc_cumul
        % -----
        % Get the past experience of the animal with
        % the patch
        Tbl(cnt, 11:12) = [hist_A_leaky(JoinEvents(iEvent,1)), hist_B_leaky(JoinEvents(iEvent,1))]; %ind_leaky
        % -----
        % Get the leaky animal density
        Tbl(cnt, 13:14) = [dens_A_leaky(JoinEvents(iEvent,1)), dens_B_leaky(JoinEvents(iEvent,1))]; %soc_leaky
        % -----
        % Get the actual animal density
        Tbl(cnt, 15:16) = [dens_A(JoinEvents(iEvent,1)), dens_B(JoinEvents(iEvent,1))]; %soc_leaky
        % -----
        % Update counter
        cnt = cnt+1;
    end
end%iEvent

% Cut
Tbl = Tbl(1:cnt-1,:);

% Prepare output
ModelData = table(...
    Tbl(:, 1),  Tbl(:, 2),  Tbl(:, 3),  Tbl(:, 4),  Tbl(:, 5),...
    Tbl(:, 6),  Tbl(:, 7),  Tbl(:, 8),  Tbl(:, 9),  Tbl(:, 10),...
    Tbl(:, 11), Tbl(:, 12), Tbl(:, 13), Tbl(:, 14), Tbl(:, 15), Tbl(:, 16),...
    'VariableNames',...
    {'choice',...           1
    'animal_ID',...         2
    'trial_ID',...          3
    'group_size',...        4
    'num_visits_A',...      5
    'num_visits_B',...      6
    'ind_cumul_A',...       7
    'ind_cumul_B',...       8
    'soc_cumul_A',...       9
    'soc_cumul_B',...      10
    'ind_leaky_A',...      11
    'ind_leaky_B',...      12
    'soc_leaky_A',...      13
    'soc_leaky_B',...      14
    'soc_curr_A',...       15
    'soc_curr_B'}...       16
    );

end%FCN:InfoIntegration_ani

function ModelData = InfoIntegration_grp(DATA, tau)
% Extract data in the correct format.

%--------------------------------------------------------------
% Setting
%--------------------------------------------------------------
SET.InteractionRange_SD = 7; %[cm]
SET.dArena = 90; %[cm]
SET.FrameRate = 25;
SET.TauInd = tau(1);
SET.TauSoc = tau(2);
SET.minSpeed = 0.25;
SET.N_bin = 7;
SET.BFnames = {'ind_num_visit'; 'ind_experience_cumul'; 'soc_density_cumul'; 'ind_experience_leaky'; 'soc_density_leaky'};
SET.ProbNames = [SET.BFnames; 'integration_cumul'; 'integration_leaky'];

%--------------------------------------------------------------
% Extract data
%--------------------------------------------------------------
% Prepare table
Tbl = nan(150000, 16);
currTrialList = fieldnames(DATA);
% Counters
cnt = 1;
cnt_ani = 0;
cnt_trial = 0;
% Iterate over trials
for iTrial = 1:length(currTrialList)
    cnt_trial = cnt_trial+1;
    % Shortcut to data
    currData = DATA.(currTrialList{iTrial});
    % Get number of animals
    grpSize = size(currData.pos_x,2);
    % Iterate over all animals
    for iAni = 1:size(currData.pos_x,2)
        cnt_ani = cnt_ani+1;
        % Create time vector
        timeVec = linspace(0, size(currData.pos_x,1)/SET.FrameRate - (1/SET.FrameRate), size(currData.pos_x,1));
        dt = 1/SET.FrameRate;
        % Create helper variables for the different patches
        helper_A = currData.AtPatch_A(:,iAni);
        helper_B = currData.AtPatch_B(:,iAni);
        % Get joining events
        helper_A_diff = diff(helper_A);
        helper_B_diff = diff(helper_B);
        join_A = find(helper_A_diff>0);
        join_B = find(helper_B_diff>0);
        % Check whether animal was there from the beginning
        if helper_A(1) == 1
            join_A = [1; join_A];
        end
        if helper_B(1) == 1
            join_B = [1; join_B];
        end
        % Concat joining event
        JoinEvents = [join_A(:) ones(length(join_A),1); join_B(:) zeros(length(join_B),1)];
        % Sort joining event
        [~,idx] = sort(JoinEvents(:,1));
        JoinEvents = JoinEvents(idx,:);
        % Get interval between joining events
        if size(JoinEvents,1)==1
            JoinEvents = [JoinEvents, NaN];
        elseif size(JoinEvents,1)>1
            JoinEvents = [JoinEvents, [NaN; diff(JoinEvents(:,1))]];
        end
        % Patch density
        dens_A = currData.PatchDensity_A;
        dens_A = dens_A - normpdf(currData.Dist2Patch_A(:,iAni), 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
        dens_B = currData.PatchDensity_B;
        dens_B = dens_B - normpdf(currData.Dist2Patch_B(:,iAni), 0, SET.InteractionRange_SD/(SET.dArena/2))/normpdf(0, 0, SET.InteractionRange_SD/(SET.dArena/2));
        % Intermittent motion
        StopBout = currData.rawSpeed(:,iAni)<SET.minSpeed;
        StopBout = StopBout.*(~(currData.AtPatch_A(:,iAni)+currData.AtPatch_B(:,iAni)));
        % Accumulated evidence (include memory extinction)
        % PERSONAL
        hist_A_cumul = cumsum(helper_A)/SET.FrameRate;
        hist_B_cumul = cumsum(helper_B)/SET.FrameRate;
        hist_A_leaky = SubFcn.leakyIntegrator(helper_A, SET.TauInd, timeVec, dt);
        hist_B_leaky = SubFcn.leakyIntegrator(helper_B, SET.TauInd, timeVec, dt);
        %SOCIAL
        dens_A_cumul = cumsum(dens_A(:).*StopBout(:));
        dens_B_cumul = cumsum(dens_B(:).*StopBout(:));
        dens_A_leaky = SubFcn.leakyIntegrator(dens_A(:).*StopBout(:), SET.TauSoc, timeVec, dt);
        dens_B_leaky = SubFcn.leakyIntegrator(dens_B(:).*StopBout(:), SET.TauSoc, timeVec, dt);

        % Iterate over consecutive joinings and combine
        % everything
        for iEvent = 1:size(JoinEvents,1)
            % Exclude animals that were feeding from the
            % beginning
            if JoinEvents(iEvent,1)~=1
                % Keep track of actual choice, animal and trial
                % identity
                Tbl(cnt, 1) = JoinEvents(iEvent,2);
                Tbl(cnt, 2) = cnt_ani;
                Tbl(cnt, 3) = cnt_trial;
                Tbl(cnt, 4) = grpSize;
                % -----
                % Get the number of past visits to each shelter
                if iEvent == 1
                    Tbl(cnt, 5:6) = [0, 0];
                else
                    if JoinEvents(iEvent-1,2) == 1
                        Tbl(cnt, 5:6) = Tbl(cnt-1,5:6) + [1,0];
                    elseif JoinEvents(iEvent-1,2) == 0
                        Tbl(cnt, 5:6) = Tbl(cnt-1,5:6) + [0,1];
                    end%if
                end%if
                % -----
                % Get the cumulative time an animal has spent
                Tbl(cnt, 7:8) = [hist_A_cumul(JoinEvents(iEvent,1)), hist_B_cumul(JoinEvents(iEvent,1))]; %ind_cumul
                % -----
                % Get the prevailing animal density
                Tbl(cnt, 9:10) = [dens_A_cumul(JoinEvents(iEvent,1)), dens_B_cumul(JoinEvents(iEvent,1))]; %soc_cumul
                % -----
                % Get the past experience of the animal with
                % the patch
                Tbl(cnt, 11:12) = [hist_A_leaky(JoinEvents(iEvent,1)), hist_B_leaky(JoinEvents(iEvent,1))]; %ind_leaky
                % -----
                % Get the leaky animal density
                Tbl(cnt, 13:14) = [dens_A_leaky(JoinEvents(iEvent,1)), dens_B_leaky(JoinEvents(iEvent,1))]; %soc_leaky
                % -----
                % Get the actual animal density
                Tbl(cnt, 15:16) = [dens_A(JoinEvents(iEvent,1)), dens_B(JoinEvents(iEvent,1))]; %soc_leaky
                % -----
                % Update counter
                cnt = cnt+1;
            end
        end%iEvent
    end%iAni
end%iTrial

% Cut
Tbl = Tbl(1:cnt-1,:);

% Prepare output
ModelData = table(...
    Tbl(:, 1),  Tbl(:, 2),  Tbl(:, 3),  Tbl(:, 4),  Tbl(:, 5),...
    Tbl(:, 6),  Tbl(:, 7),  Tbl(:, 8),  Tbl(:, 9),  Tbl(:, 10),...
    Tbl(:, 11), Tbl(:, 12), Tbl(:, 13), Tbl(:, 14), Tbl(:, 15), Tbl(:, 16),...
    'VariableNames',...
    {'choice',...           1
    'animal_ID',...         2
    'trial_ID',...          3
    'group_size',...        4
    'num_visits_A',...      5
    'num_visits_B',...      6
    'ind_cumul_A',...       7
    'ind_cumul_B',...       8
    'soc_cumul_A',...       9
    'soc_cumul_B',...      10
    'ind_leaky_A',...      11
    'ind_leaky_B',...      12
    'soc_leaky_A',...      13
    'soc_leaky_B',...      14
    'soc_curr_A',...       15
    'soc_curr_B'}...       16
    );

end%FCN:InfoIntegration_grp

function out = evalfun(data, in_a_x, in_a_y, in_s, in_Ka, in_Ks)
% Evaluate the model with the current parameter set

% Functions
% --- Quality of individual info
A_x = @(a, e_x, e_y, Ka) a.^-(e_x-Ka.*e_y);
A_y = @(a, e_x, e_y, Ka) a.^-(e_y-Ka.*e_x);
% --- Quality of social info
S_x = @(s, n_x, n_y, Ks) s.^-(n_x-Ks.*n_y);
S_y = @(s, n_x, n_y, Ks) s.^-(n_y-Ks.*n_x);
% --- Final probabilities
P_x_good = @(Ax, Sx) (1+Ax.*Sx);
P_y_good = @(Ay, Sy) (1+Ay.*Sy);


% Get data
e_x = data.ind_leaky_A(:); % Experience for patch x
e_y = data.ind_leaky_B(:); % Experience for patch y
n_x = data.soc_leaky_A(:); % Animal density at x
n_y = data.soc_leaky_B(:); % Animal density at x


% Transform data to reduce range of possible values
transform_e_x = sqrt(e_x);
transform_e_y = sqrt(e_y);
transform_n_x = sqrt(n_x);
transform_n_y = sqrt(n_y);

% Calculate parameters
calc.A_x = A_x(in_a_x, transform_e_x, transform_e_y, in_Ka);
calc.A_y = A_y(in_a_y, transform_e_x, transform_e_y, in_Ka);
calc.S_x = S_x(in_s, transform_n_x, transform_n_y, in_Ks);
calc.S_y = S_y(in_s, transform_n_x, transform_n_y, in_Ks);


% Get prediction based on probability matching
overall_P = ...
    (1 + ...
    P_x_good(calc.A_x, calc.S_x) ./ ...
    P_y_good(calc.A_y, calc.S_y)...
    ).^-1;


% Get prob. that the choice was correct
idx = find(data.choice==0);
overall_P(idx) = 1-overall_P(idx);

% Output
out = mean(ones(size(overall_P))-overall_P);

end%FCN:evalfun


function out =  exefun(data, in_a_x, in_a_y, in_s, in_Ka, in_Ks)
% Execute the model with the best parameter set

% Functions
% --- Quality of individual info
A_x = @(a, e_x, e_y, Ka) a.^-(e_x-Ka.*e_y);
A_y = @(a, e_x, e_y, Ka) a.^-(e_y-Ka.*e_x);
% --- Quality of social info
S_x = @(s, n_x, n_y, Ks) s.^-(n_x-Ks.*n_y);
S_y = @(s, n_x, n_y, Ks) s.^-(n_y-Ks.*n_x);
% --- Final probabilities
P_x_good = @(Ax, Sx) (1+Ax.*Sx);
P_y_good = @(Ay, Sy) (1+Ay.*Sy);


% Get data
e_x = data.ind_leaky_A(:); % Experience for patch x
e_y = data.ind_leaky_B(:); % Experience for patch y
n_x = data.soc_leaky_A(:); % Animal density at x
n_y = data.soc_leaky_B(:); % Animal density at x


% Transform data to reduce range of possible values
transform_e_x = sqrt(e_x);
transform_e_y = sqrt(e_y);
transform_n_x = sqrt(n_x);
transform_n_y = sqrt(n_y);


% Calculate parameters
calc.A_x = A_x(in_a_x, transform_e_x, transform_e_y, in_Ka);
calc.A_y = A_y(in_a_y, transform_e_x, transform_e_y, in_Ka);
calc.S_x = S_x(in_s, transform_n_x, transform_n_y, in_Ks);
calc.S_y = S_y(in_s, transform_n_x, transform_n_y, in_Ks);


% Get prediction based on probability matching
overall_P = ...
    (1 + ...
    P_x_good(calc.A_x, calc.S_x) ./ ...
    P_y_good(calc.A_y, calc.S_y)...
    ).^-1;


% Get prob. that the choice was correct
idx = find(data.choice==0);
overall_P(idx) = 1-overall_P(idx);

% Output
out = overall_P;

end%FCN:exefun
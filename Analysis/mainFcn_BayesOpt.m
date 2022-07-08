function [ModelData, ModelParameters, ModelOutput] = mainFcn_BayesOpt
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
% Version: 16-May-2022 (MATLAB R2022a)

% -------------------------------------------------------------------------
% Get data
% -------------------------------------------------------------------------
load('PooledData.mat');


% -------------------------------------------------------------------------
% First, get optimal time constants for each group size under both
% conditions
% -------------------------------------------------------------------------
% Iterate over group sizes
for iGrp = 1:length(SET.ConditionNames.Group)
    BestTau_ind.(SET.ConditionNames.Group{iGrp}) = nanmean([...
        PooledDATA.EQ.(SET.ConditionNames.Group{iGrp}).FeedingInterval;...
        PooledDATA.UE.(SET.ConditionNames.Group{iGrp}).FeedingInterval]);

    BestTau_soc.(SET.ConditionNames.Group{iGrp}) = nanmean([...
        PooledDATA.EQ.(SET.ConditionNames.Group{iGrp}).intermittentMotion_stand;...
        PooledDATA.UE.(SET.ConditionNames.Group{iGrp}).intermittentMotion_stand]);
end%iGrp



% -------------------------------------------------------------------------
% Use best tau to extract data
% -------------------------------------------------------------------------
for iCond = 1:length(SET.ConditionNames.Patch)
    for iGrp = 1:length(SET.ConditionNames.Group)
        % Get data
        currData = DATA.(SET.ConditionNames.Patch{iCond}).(SET.ConditionNames.Group{iGrp});
        % Get data for model
        tempData = InfoIntegration(currData, [BestTau_ind.(SET.ConditionNames.Group{iGrp}) BestTau_soc.(SET.ConditionNames.Group{iGrp})]);
        if iGrp == 1
            tempData.soc_leaky_A = zeros(size(tempData.soc_leaky_A));
            tempData.soc_leaky_B = zeros(size(tempData.soc_leaky_A));
            ModelData.(SET.ConditionNames.Patch{iCond}) = tempData;
        else
            tempData.animal_ID = tempData.animal_ID + ModelData.(SET.ConditionNames.Patch{iCond}).animal_ID(end);
            ModelData.(SET.ConditionNames.Patch{iCond}) = [ModelData.(SET.ConditionNames.Patch{iCond}); tempData];
        end%if
        clear tempData
    end%iGrp
end%iCond



% -------------------------------------------------------------------------
% Determine optimal parameters
% -------------------------------------------------------------------------
% Optimizale parameters
% --- private
SET.bayesfit.a_x =   optimizableVariable('a_x',[1 25]);
SET.bayesfit.a_y =   optimizableVariable('a_y',[1 25]);
SET.bayesfit.Ka =    optimizableVariable('Ka',[0 1]);
% --- social
SET.bayesfit.s =     optimizableVariable('s',[1 25]);
SET.bayesfit.Ks =    optimizableVariable('Ks',[0 1]);
SET.bayesfit.comp =  optimizableVariable('comp',[0 20]);

for iCond = 1:length(SET.ConditionNames.Patch)

    % Get data
    currData = ModelData.(SET.ConditionNames.Patch{iCond});

    % Get list of animals in order to fit parameters to each of them
    AnimalList = unique(currData.animal_ID);

    % Preallocation
    ModelParameters.(SET.ConditionNames.Patch{iCond}) = nan(length(AnimalList), 6);
    ModelOutput.(SET.ConditionNames.Patch{iCond}) = nan(size(currData,1), 5);

    for iAni = 1:length(AnimalList)

        % Get index position of current animal
        idx = find(currData.animal_ID == AnimalList(iAni));
        % Data for both information classes
        tempData_both = currData(idx,:);
        % Data for private info only
        tempData_ind = currData(idx,:);
        tempData_ind.soc_leaky_A = zeros(size(tempData_ind.soc_leaky_A));
        tempData_ind.soc_leaky_B = zeros(size(tempData_ind.soc_leaky_B));
        % Data for social info only
        tempData_soc = currData(idx,:);
        tempData_soc.ind_leaky_A = zeros(size(tempData_soc.ind_leaky_A));
        tempData_soc.ind_leaky_B = zeros(size(tempData_soc.ind_leaky_B));

        % Get optimizable parameters
        if iCond == 1
            vars = [SET.bayesfit.a_x, SET.bayesfit.s, SET.bayesfit.Ka, SET.bayesfit.Ks, SET.bayesfit.comp];
            fun = @(vars) evalfun(tempData_both, vars.a_x, vars.a_x, vars.s, vars.Ka, vars.Ks, vars.comp);
        else
            vars = [SET.bayesfit.a_x, SET.bayesfit.a_y, SET.bayesfit.s, SET.bayesfit.Ka, SET.bayesfit.Ks, SET.bayesfit.comp];
            fun = @(vars) evalfun(tempData_both, vars.a_x, vars.a_y, vars.s, vars.Ka, vars.Ks, vars.comp);
        end

        % -----------------------------------------------------------------
        % Determine optimal parameters
        % -----------------------------------------------------------------
        results = ...
            bayesopt(fun, vars,...
            'Verbose',0,...
            'AcquisitionFunctionName','expected-improvement-plus',...
            'IsObjectiveDeterministic',true,...
            'UseParallel',true,...
            'ExplorationRatio',0.5,...
            'PlotFcn', [],...
            'MaxObjectiveEvaluations', 100);
        if iCond == 1
            ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,:) = [results.XAtMinObjective.a_x, results.XAtMinObjective.a_x, results.XAtMinObjective.s, results.XAtMinObjective.Ka, results.XAtMinObjective.Ks, results.XAtMinObjective.comp];
        else
            ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,:) = [results.XAtMinObjective.a_x, results.XAtMinObjective.a_y, results.XAtMinObjective.s, results.XAtMinObjective.Ka, results.XAtMinObjective.Ks, results.XAtMinObjective.comp];
        end

        % -----------------------------------------------------------------
        % Get predictions
        % -----------------------------------------------------------------
        para.a_x = ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,1);
        para.a_y = ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,2);
        para.s = ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,3);
        para.Ka = ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,4);
        para.Ks = ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,5);
        para.comp = ModelParameters.(SET.ConditionNames.Patch{iCond})(iAni,6);

        ModelOutput.(SET.ConditionNames.Patch{iCond})(idx,1:2) = [tempData_both.animal_ID, tempData_both.group_size];
        % ----- BOTH
        ModelOutput.(SET.ConditionNames.Patch{iCond})(idx,3) = ...
            exefun(tempData_both, para.a_x, para.a_y, para.s, para.Ka, para.Ks, para.comp);
        % ----- IND
        ModelOutput.(SET.ConditionNames.Patch{iCond})(idx,4) = ...
            exefun(tempData_ind, para.a_x, para.a_y, para.s, para.Ka, 1, para.comp);
        % ----- SOC
        ModelOutput.(SET.ConditionNames.Patch{iCond})(idx,5) = ...
            exefun(tempData_soc, para.a_x, para.a_y, para.s, 1, para.Ks, para.comp);

        clc
        disp(['--- ', SET.ConditionNames.Patch{iCond}, ' ---'])
        disp(['Done with: ', num2str(iAni),'/',num2str(length(AnimalList)), ' animals'])

    end%iAni

    save('BayesOpt_results.mat', 'ModelData', 'ModelParameters', 'ModelOutput')

end%iCond
end


function ModelData = InfoIntegration(DATA, tau)
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
Tbl = nan(150000, 14);
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
        hist_A_leaky = SubFcn.leakyIntegrator(helper_A, SET.TauInd, timeVec, dt)/SET.FrameRate;
        hist_B_leaky = SubFcn.leakyIntegrator(helper_B, SET.TauInd, timeVec, dt)/SET.FrameRate;
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
                % Get the prevailing animal density
                Tbl(cnt, 13:14) = [dens_A_leaky(JoinEvents(iEvent,1)), dens_B_leaky(JoinEvents(iEvent,1))]; %soc_leaky
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
    Tbl(:, 11), Tbl(:, 12), Tbl(:, 13), Tbl(:, 14),...
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
    'soc_leaky_B'}...      14
    );

end%FCN:InfoIntegration

function out = evalfun(data, in_a_x, in_a_y, in_s, in_Ka, in_Ks, in_comp)
% Evaluate the model with the current parameter set

% Functions
% --- Quality of individual info
A_x = @(a, e_x, e_y, Ka) a.^-(e_x-Ka.*e_y);
A_y = @(a, e_x, e_y, Ka) a.^-(e_y-Ka.*e_x);
% --- Quality of social info
S_x = @(s, n_x, n_y, Ks, comp) s.^(comp.*(n_x-Ks.*n_y));
S_y = @(s, n_x, n_y, Ks, comp) s.^(comp.*(n_y-Ks.*n_x));
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


% Check whether to avoid competition
comp = -ones(size(transform_e_x));
for iChoice = 1:length(transform_e_x)
    if (transform_e_x(iChoice)>0) && (transform_e_y(iChoice)>0)
        comp(iChoice) = 2*((transform_n_x(iChoice)>in_comp)-0.5);
    end%if personal info
end%iChoice


% Calculate parameters
calc.A_x = A_x(in_a_x, transform_e_x, transform_e_y, in_Ka);
calc.A_y = A_y(in_a_y, transform_e_x, transform_e_y, in_Ka);
calc.S_x = S_x(in_s, transform_n_x, transform_n_y, in_Ks, comp);
calc.S_y = S_y(in_s, transform_n_x, transform_n_y, in_Ks, comp);


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


function out =  exefun(data, in_a_x, in_a_y, in_s, in_Ka, in_Ks, in_comp)
% Execute the model with the best parameter set

% Functions
% --- Quality of individual info
A_x = @(a, e_x, e_y, Ka) a.^-(e_x-Ka.*e_y);
A_y = @(a, e_x, e_y, Ka) a.^-(e_y-Ka.*e_x);
% --- Quality of social info
S_x = @(s, n_x, n_y, Ks, comp) s.^(comp.*(n_x-Ks.*n_y));
S_y = @(s, n_x, n_y, Ks, comp) s.^(comp.*(n_y-Ks.*n_x));
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


% Check whether to avoid competition
comp = -ones(size(transform_e_x));
for iChoice = 1:length(transform_e_x)
    if (transform_e_x(iChoice)>0) && (transform_e_y(iChoice)>0)
        comp(iChoice) = 2*((transform_n_x(iChoice)>in_comp)-0.5);
    end%if personal info
end%iChoice


% Calculate parameters
calc.A_x = A_x(in_a_x, transform_e_x, transform_e_y, in_Ka);
calc.A_y = A_y(in_a_y, transform_e_x, transform_e_y, in_Ka);
calc.S_x = S_x(in_s, transform_n_x, transform_n_y, in_Ks, comp);
calc.S_y = S_y(in_s, transform_n_x, transform_n_y, in_Ks, comp);


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
classdef SubFcn_LocustDecisionSystem
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
    % This is file contains all helper functions for inference of the
    % locust decision system
    %
    % Version: 15-Jan-2022 (MATLAB R2022a)

    properties
    end

    methods(Static)

        function out = evalfun(data, in_q_x, in_q_y, in_s, in_Kq, in_Ks)
            % Evaluate the model with the current parameter set

            % Probability of a patch being a good choice
            % --- Private information
            P_good_q_x = @(q_x, e_x, e_y, Kq) (1 + q_x.^-(e_x-Kq*e_y)).^-1;
            P_good_q_y = @(q_y, e_x, e_y, Kq) (1 + q_y.^-(e_y-Kq*e_x)).^-1;
            % --- Social information
            P_good_s_x = @(s, n_x, n_y, Ks) (1 + s.^-(n_x-Ks*n_y)).^-1;
            P_good_s_y = @(s, n_x, n_y, Ks) (1 + s.^-(n_y-Ks*n_x)).^-1;
            % Probability matching
            P = @(Px_q, Px_s, Py_q, Py_s) (Px_q.*Px_s) ./ (Px_q.*Px_s + Py_q.*Py_s);


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


            % Calculate individual probabiities of a patch being a good
            % choice
            % --- x is good
            Prob.Px_q = P_good_q_x(in_q_x, transform_e_x, transform_e_y, in_Kq);
            Prob.Px_s = P_good_s_x(in_s, transform_n_x, transform_n_y, in_Ks);
            % --- y is good
            Prob.Py_q = P_good_q_y(in_q_y, transform_e_x, transform_e_y, in_Kq);
            Prob.Py_s = P_good_s_y(in_s, transform_n_x, transform_n_y, in_Ks);


            % Get prediction based on probability matching
            overall_P = P(Prob.Px_q, Prob.Px_s, Prob.Py_q, Prob.Py_s);


            % Get prob. that the choice was correct
            idx = find(data.choice==0);
            overall_P(idx) = 1-overall_P(idx);


            % Output
            % out = mean(ones(size(overall_P))-overall_P); % avg
            out = sum((ones(size(overall_P))-overall_P).^2); % RSS
            out = log1p(out);

        end%FCN:evalfun

        function out =  exefun(data, in_q_x, in_q_y, in_s, in_Kq, in_Ks)

            % Execute the model with the best parameter set

            % Probability of a patch being a good choice
            % --- Private information
            P_good_q_x = @(q_x, e_x, e_y, Kq) (1 + q_x.^-(e_x-Kq*e_y)).^-1;
            P_good_q_y = @(q_y, e_x, e_y, Kq) (1 + q_y.^-(e_y-Kq*e_x)).^-1;
            % --- Social information
            P_good_s_x = @(s, n_x, n_y, Ks) (1 + s.^-(n_x-Ks*n_y)).^-1;
            P_good_s_y = @(s, n_x, n_y, Ks) (1 + s.^-(n_y-Ks*n_x)).^-1;
            % Probability matching
            P = @(Px_q, Px_s, Py_q, Py_s) (Px_q.*Px_s) ./ (Px_q.*Px_s + Py_q.*Py_s);


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


            % Calculate individual probabiities of a patch being a good
            % choice
            % --- x is good
            Prob.Px_q = P_good_q_x(in_q_x, transform_e_x, transform_e_y, in_Kq);
            Prob.Px_s = P_good_s_x(in_s, transform_n_x, transform_n_y, in_Ks);
            % --- y is good
            Prob.Py_q = P_good_q_y(in_q_y, transform_e_x, transform_e_y, in_Kq);
            Prob.Py_s = P_good_s_y(in_s, transform_n_x, transform_n_y, in_Ks);


            % Get prediction based on probability matching
            overall_P = P(Prob.Px_q, Prob.Px_s, Prob.Py_q, Prob.Py_s);


            % Get prob. that the choice was correct
            idx = find(data.choice==0);
            overall_P(idx) = 1-overall_P(idx);

            % Output
            out = overall_P;

        end%FCN:exefun

        function out =  exefun2(data, in_b_x, in_b_y, in_a, in_s, in_Ka, in_Ks, in_rand)

            % Execute the model with the best parameter set

            % Functions
            % --- Quality of individual info
            A_x = @(a, e_x, e_y, Ka) a.^-(e_x-Ka.*e_y);
            A_y = @(a, e_x, e_y, Ka) a.^-(e_y-Ka.*e_x);
            % --- Quality of social info
            S_x = @(s, n_x, n_y, Ks) s.^-(n_x-Ks.*n_y);
            S_y = @(s, n_x, n_y, Ks) s.^-(n_y-Ks.*n_x);
            % --- Final probabilities
            P_x_good = @(b, Ax, Sx) (1+b.*Ax.*Sx);
            P_y_good = @(b, Ay, Sy) (1+b.*Ay.*Sy);


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
            calc.A_x = A_x(in_a, transform_e_x, transform_e_y, in_Ka);
            calc.A_y = A_y(in_a, transform_e_x, transform_e_y, in_Ka);
            calc.S_x = S_x(in_s, transform_n_x, transform_n_y, in_Ks);
            calc.S_y = S_y(in_s, transform_n_x, transform_n_y, in_Ks);


            % Get prediction based on probability matching
            overall_P = ...
                (in_rand/2) + (1-in_rand) .* ...
                (1 + ...
                P_x_good(in_b_x, calc.A_x, calc.S_x) ./ ...
                P_y_good(in_b_y, calc.A_y, calc.S_y)...
                ).^-1;

            % Output
            out = overall_P;

        end%FCN:exefun

        function ModelData = InfoIntegration(DATA, tau)
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
            cnt_ani = 0;
            cnt_trial = 0;
            % Iterate over all trials
            TrialList = fieldnames(DATA);
            for iTrial = 1:length(TrialList)
                cnt_trial = cnt_trial+1;
                % Shortcut to data
                currData = DATA.(TrialList{iTrial});
                % Get number of animals
                grpSize = size(currData.pos_x,2);
                % Iterate over all animals
                for iAni = 1:size(currData.pos_x,2)
                    cnt_curr = 1;
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
                    % Keep current density
                    dens_curr_A = dens_A;
                    dens_curr_B = dens_B;
                    % Normalize with distance to patch
                    dens_A = dens_A./currData.Dist2Patch_A(:,iAni);
                    dens_B = dens_B./currData.Dist2Patch_B(:,iAni);
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
                            if cnt_curr == 1
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
                            Tbl(cnt, 15:16) = [dens_curr_A(JoinEvents(iEvent,1)), dens_curr_B(JoinEvents(iEvent,1))]; %soc_leaky
                            % -----
                            % Update counter
                            cnt = cnt+1;
                            cnt_curr = cnt_curr+1;
                        end%if first event
                    end%iEvent
                end%iAni
            end%iTrial
            % Cut
            Tbl = Tbl(1:cnt-1,:);
            %--------------------------------------------------------------
            % Prepare output
            %--------------------------------------------------------------
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
        end%FCN:InfoIntegration_trial

    end%method
end%class
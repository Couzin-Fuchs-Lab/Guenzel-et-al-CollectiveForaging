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
% Version: 16-May-2022 (MATLAB R2022a)

% Tidy up
clear all
close all
clc
% Add paths
addpath(genpath(pwd))
% Load data
load('PooledData.mat')
load('BayesOpt_results.mat')

%% Plot everything

figure('units', 'normalized', 'Position', [0.1 0.1 0.1 0.5])

for iCond = 1:length(SET.ConditionNames.Patch)
    for iGrp = 1:length(SET.ConditionNames.Group)
        
        figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on
        
        idxGrp = find(ModelOutput.(SET.ConditionNames.Patch{iCond})(:,2) == str2double(SET.ConditionNames.Group{iGrp}(2:end)));
        dat_both = ModelOutput.(SET.ConditionNames.Patch{iCond})(idxGrp,3);
        dat_ind =  ModelOutput.(SET.ConditionNames.Patch{iCond})(idxGrp,4);
        dat_soc =  ModelOutput.(SET.ConditionNames.Patch{iCond})(idxGrp,5);
        
        
        boot_ind = bootstrp(5000,@mean, dat_ind);
        boot_soc = bootstrp(5000,@mean, dat_soc);
        boot_both = bootstrp(5000,@mean, dat_both);
        
        properties.NumPoints =          1000;               % Points at which to evaluate the probability density estimate
        properties.MinVal =             min(dat_ind);       % Smallest possible value (e.g. errors = 0)
        properties.MaxVal =             max(dat_ind);       % Biggest possible value
        properties.Orientation =        'vertical';         %(Set how the box plot should be oriented 'vertical' or 'horizontal')
        properties.AvgType =            'mean';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
        properties.EdgeCol =            'r';                %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
        properties.EdgeWidth =          2;                  %(Set the width of violin's edges. If EdgeCol is 'none', this statement will be ignored)
        properties.MeanCol =            'r';                %(Set colour of mean. If it should be depicted)
        properties.MeanSymbol =         'line';             %(Set how the mean should be depicted as symbol (e.g. 'x') or as line: 'line')
        properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
        properties.SeparateOutliers =   0;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
        SubFcn.violinplot_advanced(boot_ind, 1, 0.75, properties)
        plot(1.5, mean(dat_ind>0.5), 'rd')
        
        properties.NumPoints =          1000;               % Points at which to evaluate the probability density estimate
        properties.MinVal =             min(dat_soc);       % Smallest possible value (e.g. errors = 0)
        properties.MaxVal =             max(dat_soc);       % Biggest possible value
        properties.Orientation =        'vertical';         %(Set how the box plot should be oriented 'vertical' or 'horizontal')
        properties.AvgType =            'mean';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
        properties.EdgeCol =            'b';                %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
        properties.EdgeWidth =          2;                  %(Set the width of violin's edges. If EdgeCol is 'none', this statement will be ignored)
        properties.MeanCol =            'b';                %(Set colour of mean. If it should be depicted)
        properties.MeanSymbol =         'line';             %(Set how the mean should be depicted as symbol (e.g. 'x') or as line: 'line')
        properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
        properties.SeparateOutliers =   0;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
        SubFcn.violinplot_advanced(boot_soc, 2, 0.75, properties)
        plot(2.5, mean(dat_soc>0.5), 'bd')
        
        properties.NumPoints =          1000;               % Points at which to evaluate the probability density estimate
        properties.MinVal =             min(dat_both);       % Smallest possible value (e.g. errors = 0)
        properties.MaxVal =             max(dat_both);       % Biggest possible value
        properties.Orientation =        'vertical';         %(Set how the box plot should be oriented 'vertical' or 'horizontal')
        properties.AvgType =            'mean';             %(Set which measure should be plotted: 'median', 'mean' or 'both')
        properties.EdgeCol =            'k';                %(Color classification of the violin's edge. Set to 'none' if no edge should be depicted)
        properties.EdgeWidth =          2;                  %(Set the width of violin's edges. If EdgeCol is 'none', this statement will be ignored)
        properties.MeanCol =            'k';                %(Set colour of mean. If it should be depicted)
        properties.MeanSymbol =         'line';             %(Set how the mean should be depicted as symbol (e.g. 'x') or as line: 'line')
        properties.MeanWidth =          2;                  %(Set width of the median's line and/or symbol)
        properties.SeparateOutliers =   0;                  %(Logical statement whether to exclude outliers from the violin and depict them as individual data points)
        SubFcn.violinplot_advanced(boot_both, 3, 0.75, properties)
        plot(3.5, mean(dat_both>0.5), 'kd')
        
        title([SET.ConditionNames.Patch{iCond}, ' - ', SET.ConditionNames.Group{iGrp}])
        ylim([0.4 0.9])
        ylabel('P correct')
        xticks([])
        xlim([0 4])
        export_fig(['FIG\',SET.ConditionNames.Patch{iCond}, '_', SET.ConditionNames.Group{iGrp}, '_violins'], '-pdf')
        
    end%iGrp
    
    dat_both = ModelOutput.(SET.ConditionNames.Patch{iCond})(:,3);
    dat_ind =  ModelOutput.(SET.ConditionNames.Patch{iCond})(:,4);
    dat_soc =  ModelOutput.(SET.ConditionNames.Patch{iCond})(:,5);
    
    
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
    
    
    col_reinforce = [255 113 12]/255;
    col_balance = [0 142 243]/255;
    
    figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on
    plot(mean([dat_ind(idx_reinforce), dat_soc(idx_reinforce)],2), dat_both(idx_reinforce),'.', 'color', col_reinforce)
    plot(mean([dat_ind(idx_balance), dat_soc(idx_balance)],2), dat_both(idx_balance),'.', 'color', col_balance)
    plot([0 1], [0 1], 'k')
    plot([0.5 0.5], [0 1], 'k')
    plot([0 1], [0.5 0.5], 'k')
    axis equal
    xlim([0 1])
    ylim([0 1])
    xlabel('mean(ind,soc)')
    title(SET.ConditionNames.Patch{iCond})
    ylabel('both')
    export_fig(['FIG\',SET.ConditionNames.Patch{iCond}, '_mean'], '-pdf')
    
    
    figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on
    plot(max([dat_ind(idx_reinforce_max), dat_soc(idx_reinforce_max)],[],2), dat_both(idx_reinforce_max),'.', 'color', col_reinforce)
    plot(min([dat_ind(idx_reinforce_min), dat_soc(idx_reinforce_min)],[],2), dat_both(idx_reinforce_min),'.', 'color', col_reinforce)
    plot(max([dat_ind(idx_balance), dat_soc(idx_balance)],[],2), dat_both(idx_balance),'.', 'color', col_balance)
    plot([0 1], [0 1], 'k')
    plot([0.5 0.5], [0 1], 'k')
    plot([0 1], [0.5 0.5], 'k')
    axis equal
    xlim([0 1])
    ylim([0 1])
    xlabel('max(ind,soc)')
    title(SET.ConditionNames.Patch{iCond})
    ylabel('both')
    export_fig(['FIG\',SET.ConditionNames.Patch{iCond}, '_max'], '-pdf')
    
   figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on
    plot(min([dat_ind(idx_reinforce_min), dat_soc(idx_reinforce_min)],[],2), dat_both(idx_reinforce_min),'.', 'color', col_reinforce)
%     plot(min([dat_ind(idx_balance), dat_soc(idx_balance)],[],2), dat_both(idx_balance),'.', 'color', col_balance)
    plot([0 1], [0 1], 'k')
    plot([0.5 0.5], [0 1], 'k')
    plot([0 1], [0.5 0.5], 'k')
    axis equal
    xlim([0 1])
    ylim([0 1])
    xlabel('min(ind,soc)')
    title(SET.ConditionNames.Patch{iCond})
    ylabel('both')
    export_fig(['FIG\',SET.ConditionNames.Patch{iCond}, '_min'], '-pdf')
    
    figure('units', 'centimeters', 'Position', [5 5 15 30]); hold on
    cmap_val = gray(1000);
    cmap_bin = linspace(0,1,size(cmap_val,1));
    bins = linspace(0,1,size(cmap_val,1));
    map = zeros(size(cmap_val,1));
    for iBin = 1:size(cmap_val,1)
        map(iBin,:) = mean([ones(1,size(cmap_val,1))*bins(iBin); linspace(0,1,size(cmap_val,1))]);
    end
    h = surf(map);
    h.EdgeColor = 'none';
    colormap(gray)
    for iPt = 1:length(dat_both)
        [~, idx_ind] = min(abs(cmap_bin-dat_ind(iPt)));
        [~, idx_soc] = min(abs(cmap_bin-dat_soc(iPt)));
        [~, idx_both] = min(abs(cmap_bin-dat_both(iPt)));
        plot3(idx_ind, idx_soc, 2, 'o', 'MarkerFaceColor', cmap_val(idx_both,:), 'MarkerEdgeColor', 'none')
    end
    axis equal
    xlim([1 1000])
    ylim([1 1000])
    set(gca, 'XTick', [1 500 1000], 'XTickLabel', [0 0.5 1])
    set(gca, 'YTick', [1 500 1000], 'YTickLabel', [0 0.5 1])
    xlabel('ind')
    title(SET.ConditionNames.Patch{iCond})
    ylabel('soc')
    export_fig(['FIG\',SET.ConditionNames.Patch{iCond}, '_colormap'], '-png')
    
end%iCond

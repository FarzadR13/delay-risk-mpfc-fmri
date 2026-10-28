function SPM = glm_heatmap(participant_id, prep_root, out_root, estimate)
% GLM_HEATMAP  Nine categorical regressors per condition: 3 probability bins x 3 distance bins (Figures 5A,B and 6B).
%   cell = 3 * (distance bin - 1) + probability bin; bins: p {.05,.10 | .25,.50,.75 | .90,.95},
%   distance {-30,-20 | -10,0,+10 | +20,+30}. Each contrast sums the two runs of a condition.
if nargin < 4, estimate = true; end
contrasts = cell(18, 2);
k = 0;
for cond = {'immediate', 'delayed'}
    for n = 1:9
        k = k + 1;
        contrasts(k, :) = {sprintf('heatmap_%s_cell%d', cond{1}, n), {sprintf('%s_cell%d\*bf', cond{1}, n), 1}};
    end
end
SPM = run_first_level('heatmap', participant_id, prep_root, fullfile(out_root, 'heatmap', participant_id), contrasts, estimate);
end

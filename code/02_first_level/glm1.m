function SPM = glm1(participant_id, prep_root, out_root, estimate)
% GLM1  Signed distance to CE and objective probability as parametric modulators (Figure 4).
%   SPM = glm1('sub-01', prep_root, out_root, true)   writes <out_root>/glm1/sub-01
if nargin < 4, estimate = true; end
contrasts = {
    'glm1_immediate_distance',               {'immediatexdistance\^1', 1}
    'glm1_delayed_distance',                 {'delayedxdistance\^1', 1}
    'glm1_immediate_minus_delayed_distance', {'immediatexdistance\^1', 1; 'delayedxdistance\^1', -1}
    'glm1_immediate_prob',                   {'immediatexprob\^1', 1}
    'glm1_delayed_prob',                     {'delayedxprob\^1', 1}};
SPM = run_first_level('glm1', participant_id, prep_root, fullfile(out_root, 'glm1', participant_id), contrasts, estimate);
end

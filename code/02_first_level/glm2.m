function SPM = glm2(participant_id, prep_root, out_root, estimate)
% GLM2  Signed distance to CE and subjective probability weight w(p) = CE / 20 (Figures 6 and 7).
%   SPM = glm2('sub-01', prep_root, out_root, true)   writes <out_root>/glm2/sub-01
if nargin < 4, estimate = true; end
contrasts = {
    'glm2_immediate_distance',         {'immediatexdistance\^1', 1}
    'glm2_delayed_distance',           {'delayedxdistance\^1', 1}
    'glm2_immediate_wp',               {'immediatexwp\^1', 1}
    'glm2_delayed_wp',                 {'delayedxwp\^1', 1}
    'glm2_delayed_minus_immediate_wp', {'delayedxwp\^1', 1; 'immediatexwp\^1', -1}};
SPM = run_first_level('glm2', participant_id, prep_root, fullfile(out_root, 'glm2', participant_id), contrasts, estimate);
end

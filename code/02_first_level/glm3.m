function SPM = glm3(participant_id, prep_root, out_root, estimate)
% GLM3  One onset regressor per run, no parametric modulators (Figure 5C,D).
%   The condition main effect sums the two runs of the same temporal condition.
if nargin < 4, estimate = true; end
contrasts = {
    'glm3_immediate_main_effect', {'immediate\*bf', 1}
    'glm3_delayed_main_effect',   {'delayed\*bf', 1}};
SPM = run_first_level('glm3', participant_id, prep_root, fullfile(out_root, 'glm3', participant_id), contrasts, estimate);
end

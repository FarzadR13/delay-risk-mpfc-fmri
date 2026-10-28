function one_sample_ttests(contrast_root, out_root)
% ONE_SAMPLE_TTESTS  Group one-sample t-tests (N = 22) for every first-level contrast reported in the paper.
%   one_sample_ttests(contrast_root, out_root)
%   contrast_root: folder with <participant_id>/<participant_id>_<contrast>.nii.gz (participant-level contrast images)
%   out_root:      one SPM folder per contrast is written to <out_root>/<contrast>
%   Contrast 1 tests a positive effect, contrast 2 a negative effect (the paper's "NEG" maps).
names = {'glm1_immediate_distance', 'glm1_delayed_distance', 'glm1_immediate_minus_delayed_distance', ...
         'glm2_immediate_distance', 'glm2_immediate_wp', 'glm2_delayed_wp', 'glm2_delayed_minus_immediate_wp'};
spm('defaults', 'fmri'); spm_jobman('initcfg');
for i = 1:numel(names)
    d = fullfile(out_root, names{i});
    scans = unzip_inputs(contrast_root, names{i}, fullfile(d, 'inputs'));
    estimate_group(d, scans, [], {'positive', 1; 'negative', -1});
end
end

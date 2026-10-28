function group_tables(sl_root, out_dir)
% GROUP_TABLES  Cluster tables (p < 0.001, cluster-level pFWE) for every group map reported in the paper,
% and the small-volume correction of the immediate distance effect within the meta-analytic mPFC mask.
%   group_tables(sl_root, out_dir)
%   sl_root: output of one_sample_ttests.m and low_p_covariate_model.m
%   Requires results/masks/bartra_mpfc.nii.gz (code/00_masks/build_masks.py).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
% {table name, model folder, contrast index}
items = {'glm1_immediate_distance_positive',               'glm1_immediate_distance', 1
         'glm1_delayed_distance_positive',                 'glm1_delayed_distance', 1
         'glm1_delayed_distance_negative',                 'glm1_delayed_distance', 2
         'glm1_immediate_minus_delayed_distance_positive', 'glm1_immediate_minus_delayed_distance', 1
         'glm2_delayed_wp_negative',                       'glm2_delayed_wp', 2
         'glm2_immediate_wp_positive',                     'glm2_immediate_wp', 1
         'glm2_immediate_wp_negative',                     'glm2_immediate_wp', 2
         'glm2_delayed_more_negative_wp',                  'glm2_delayed_minus_immediate_wp', 2
         'glm2_delayed_less_negative_wp',                  'glm2_delayed_minus_immediate_wp', 1
         'glm2_low_p_delta_w_covariate',                   'glm2_low_p_delta_w_covariate', 1};
for i = 1:size(items, 1)
    cluster_table(fullfile(sl_root, items{i, 2}), items{i, 3}, fullfile(out_dir, [items{i, 1} '_clusters.csv']));
end
mask = gunzip(fullfile(root, 'results', 'masks', 'bartra_mpfc.nii.gz'), tempdir);
cluster_table(fullfile(sl_root, 'glm1_immediate_distance'), 1, ...
              fullfile(out_dir, 'glm1_immediate_distance_positive_svc_bartra.csv'), mask{1});
end

function low_p_covariate_model(contrast_root, out_root)
% LOW_P_COVARIATE_MODEL  Whole-brain covariate analysis of low-p delta w (Figure 7B), N = 22, df = 20.
%   Input: GLM2 delayed w(p) contrast images. Covariate: low-p delta w, the mean delayed-minus-immediate
%   empirical weight (CE / 20) at p = 0.05 and 0.10, mean-centred. The paper's map tests the covariate
%   on the NEGATIVE delayed w(p) contrast; with positive images this is contrast 'negative_wp_by_low_p'.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
P = readtable(fullfile(root, 'data', 'participants.tsv'), 'FileType', 'text', 'TextType', 'string');
ids = P.participant_id(P.in_fmri_sample == "yes");
W = readtable(fullfile(root, 'data', 'behavior', 'prescan_ce.csv'), 'TextType', 'string');
dw = zeros(numel(ids), 1);
for i = 1:numel(ids)
    w = @(c, p) W.w_empirical(W.participant_id == ids(i) & W.condition == c & abs(W.probability - p) < 1e-9);
    dw(i) = mean([w("delayed", 0.05) - w("immediate", 0.05), w("delayed", 0.10) - w("immediate", 0.10)]);
end
spm('defaults', 'fmri'); spm_jobman('initcfg');
d = fullfile(out_root, 'glm2_low_p_delta_w_covariate');
scans = unzip_inputs(contrast_root, 'glm2_delayed_wp', fullfile(d, 'inputs'));
estimate_group(d, scans, dw, {'negative_wp_by_low_p', [0 -1]; 'positive_wp_by_low_p', [0 1]});
end

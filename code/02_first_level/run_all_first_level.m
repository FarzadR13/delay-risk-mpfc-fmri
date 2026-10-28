function run_all_first_level(prep_root, out_root)
% RUN_ALL_FIRST_LEVEL  Estimate the four first-level models for the 22 participants of the fMRI sample.
%   run_all_first_level('/path/to/preprocessed', '/path/to/output')
%   Requires SPM12 on the MATLAB path. Preprocessed images are not shared (see README).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
P = readtable(fullfile(root, 'data', 'participants.tsv'), 'FileType', 'text', 'TextType', 'string');
ids = P.participant_id(P.in_fmri_sample == "yes");
spm('defaults', 'fmri'); spm_jobman('initcfg');
for i = 1:numel(ids)
    id = char(ids(i));
    glm1(id, prep_root, out_root); glm2(id, prep_root, out_root);
    glm3(id, prep_root, out_root); glm_heatmap(id, prep_root, out_root);
end
end

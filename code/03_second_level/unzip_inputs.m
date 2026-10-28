function scans = unzip_inputs(contrast_root, name, tmp_dir)
% UNZIP_INPUTS  Decompress the 22 participant images of one contrast (SPM12 reads .nii only).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
P = readtable(fullfile(root, 'data', 'participants.tsv'), 'FileType', 'text', 'TextType', 'string');
ids = P.participant_id(P.in_fmri_sample == "yes");
if ~exist(tmp_dir, 'dir'), mkdir(tmp_dir); end
scans = cell(numel(ids), 1);
for i = 1:numel(ids)
    f = gunzip(fullfile(contrast_root, char(ids(i)), sprintf('%s_%s.nii.gz', ids(i), name)), tmp_dir);
    scans{i} = [f{1} ',1'];
end
end

function SPM = run_first_level(model, participant_id, prep_root, out_dir, contrasts, estimate)
% RUN_FIRST_LEVEL  Specify (and optionally estimate) one first-level SPM12 model.
%
%   model          'glm1' | 'glm2' | 'glm3' | 'heatmap' (see make_events.m)
%   participant_id release ID, e.g. 'sub-01'
%   prep_root      folder with preprocessed data, laid out as
%                  <prep_root>/<participant_id>/run-<r>/swr*.nii  (smoothed, normalized 3D volumes)
%                  <prep_root>/<participant_id>/run-<r>/rp_*.txt  (six realignment parameters)
%                  These images are not shared (see README).
%   out_dir        output folder for this participant and model
%   contrasts      cell array {name, {regexp, weight; ...}}: regressors whose name matches regexp get weight
%   estimate       false: design only (SPM.mat with xX.X); true: estimate and compute contrasts
%
%   Settings, read from the original SPM.mat files: canonical HRF, stick functions (duration 0), units
%   seconds, TR = 2 s, microtime resolution 42 and onset 21, no orthogonalization of parametric
%   modulators, six motion regressors per run, 128-s high-pass filter, AR(1), implicit mask threshold
%   0.8, no global normalisation.

if ~exist(out_dir, 'dir'), mkdir(out_dir); end
sess = make_events(participant_id, model);

job.dir = {out_dir};
job.timing.units = 'secs';
job.timing.RT = 2;
job.timing.fmri_t = 42;
job.timing.fmri_t0 = 21;
for k = 1:numel(sess)
    run_dir = fullfile(prep_root, participant_id, sprintf('run-%d', sess(k).run));
    job.sess(k).scans = cellstr(spm_select('ExtFPList', run_dir, '^swr.*\.nii$', 1));
    for c = 1:numel(sess(k).names)
        job.sess(k).cond(c).name = sess(k).names{c};
        job.sess(k).cond(c).onset = sess(k).onsets{c};
        job.sess(k).cond(c).duration = sess(k).durations{c};
        job.sess(k).cond(c).tmod = 0;
        job.sess(k).cond(c).pmod = sess(k).pmod;
        job.sess(k).cond(c).orth = 0;
    end
    job.sess(k).multi = {''};
    job.sess(k).regress = struct('name', {}, 'val', {});
    job.sess(k).multi_reg = cellstr(spm_select('FPList', run_dir, '^rp_.*\.txt$'));
    job.sess(k).hpf = 128;
end
job.fact = struct('name', {}, 'levels', {});
job.bases.hrf.derivs = [0 0];
job.volt = 1;
job.global = 'None';
job.mthresh = 0.8;
job.mask = {''};
job.cvi = 'AR(1)';

matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec = job;
spm_jobman('run', matlabbatch);
load(fullfile(out_dir, 'SPM.mat'), 'SPM');
if ~estimate, return; end

clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(out_dir, 'SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{2}.spm.stats.con.spmmat = {fullfile(out_dir, 'SPM.mat')};
for c = 1:size(contrasts, 1)
    w = zeros(1, numel(SPM.xX.name));
    rules = contrasts{c, 2};
    for r = 1:size(rules, 1)
        w(~cellfun(@isempty, regexp(SPM.xX.name, rules{r, 1}))) = rules{r, 2};
    end
    matlabbatch{2}.spm.stats.con.consess{c}.tcon.name = contrasts{c, 1};
    matlabbatch{2}.spm.stats.con.consess{c}.tcon.weights = w;
    matlabbatch{2}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
end
matlabbatch{2}.spm.stats.con.delete = 1;
spm_jobman('run', matlabbatch);
load(fullfile(out_dir, 'SPM.mat'), 'SPM');
end

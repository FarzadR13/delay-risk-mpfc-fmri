function estimate_group(out_dir, scans, covariate, contrasts)
% ESTIMATE_GROUP  One-sample t-test (optionally with one mean-centred covariate), estimation and t contrasts.
if exist(fullfile(out_dir, 'SPM.mat'), 'file'), delete(fullfile(out_dir, 'SPM.mat')); end
d.dir = {out_dir};
d.des.t1.scans = scans;
d.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
if ~isempty(covariate)
    d.cov(1) = struct('c', covariate, 'cname', 'low_p_delta_w', 'iCFI', 1, 'iCC', 1);
end
d.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
d.masking.tm.tm_none = 1;
d.masking.im = 1;
d.masking.em = {''};
d.globalc.g_omit = 1;
d.globalm.gmsca.gmsca_no = 1;
d.globalm.glonorm = 1;
matlabbatch{1}.spm.stats.factorial_design = d;
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(out_dir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat = {fullfile(out_dir, 'SPM.mat')};
for c = 1:size(contrasts, 1)
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = contrasts{c, 1};
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = contrasts{c, 2};
    matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
end
matlabbatch{3}.spm.stats.con.delete = 1;
spm_jobman('run', matlabbatch);
end

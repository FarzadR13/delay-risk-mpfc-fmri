function T = cluster_table(spm_dir, contrast, out_csv, svc_mask)
% CLUSTER_TABLE  SPM results table at p < 0.001 uncorrected (k = 0), written as CSV.
%   T = cluster_table(spm_dir, contrast, out_csv)            whole-brain table
%   T = cluster_table(spm_dir, contrast, out_csv, svc_mask)  small-volume correction within a mask image
%   Columns: cluster-level pFWE and extent k (voxels), peak-level pFWE, T, Z, uncorrected p, MNI x y z.
%   The table lists up to 3 peaks per cluster at least 8 mm apart (SPM defaults).
spm('defaults', 'fmri'); spm_get_defaults('cmdline', true);
x = struct('swd', spm_dir, 'title', '', 'Ic', contrast, 'n', 1, 'Im', [], 'pm', [], 'Ex', 0, ...
           'u', 0.001, 'k', 0, 'thresDesc', 'none');
[SPM, xSPM] = spm_getSPM(x);
if nargin == 4
    xSPM = small_volume(SPM, xSPM, svc_mask);
end
tab = spm_list('Table', xSPM);
rows = {};
for i = 1:size(tab.dat, 1)
    r = tab.dat(i, :);
    rows(end + 1, :) = {val(r{3}), val(r{5}), val(r{7}), val(r{9}), val(r{10}), val(r{11}), ...
                        r{12}(1), r{12}(2), r{12}(3)}; %#ok<AGROW>
end
if isempty(rows), rows = cell(0, 9); end
T = cell2table(rows, 'VariableNames', {'cluster_pFWE', 'cluster_k', 'peak_pFWE', 'peak_T', 'peak_Z', ...
                                       'peak_p_unc', 'x', 'y', 'z'});
writetable(T, out_csv);
end

function xSPM = small_volume(SPM, xSPM, mask_file)
% Restrict the search volume to a mask image, as SPM's small-volume correction (spm_VOI) does:
% voxels and resels of the mask define the search volume for peak- and cluster-level FWE p values.
D = spm_vol(mask_file);
inmask = @(xyz) spm_sample_vol(D, xyz(1, :), xyz(2, :), xyz(3, :), 0) > 0;
vox = D.mat \ [xSPM.XYZmm; ones(1, size(xSPM.XYZmm, 2))];
j = find(inmask(vox));                                  % suprathreshold voxels inside the mask
all_mm = SPM.xVol.M(1:3, :) * [SPM.xVol.XYZ; ones(1, SPM.xVol.S)];
k = find(inmask(D.mat \ [all_mm; ones(1, size(all_mm, 2))]));   % analysed voxels inside the mask
fwhm = xSPM.FWHM .* (xSPM.VOX ./ sqrt(sum(D.mat(1:3, 1:3) .^ 2)));
xSPM.S = numel(k);
xSPM.R = spm_resels(fwhm, D, 'I');
xSPM.Z = xSPM.Z(j); xSPM.XYZ = xSPM.XYZ(:, j); xSPM.XYZmm = xSPM.XYZmm(:, j);
if isfield(xSPM, 'Ps') && ~isempty(xSPM.Ps), xSPM.Ps = xSPM.Ps(k); end
[up, xSPM.Pp] = spm_uc_peakFDR(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, xSPM.Vspm, k, xSPM.u);
uu = spm_uc(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, xSPM.S);
[uc, xSPM.Pc, ue] = spm_uc_clusterFDR(0.05, xSPM.df, xSPM.STAT, xSPM.R, xSPM.n, xSPM.Vspm, k, ...
                                      1 / prod(xSPM.FWHM(xSPM.DIM > 1)), xSPM.u);
xSPM.uc = [uu up ue uc];
end

function v = val(x)
if isempty(x), v = NaN; else, v = x; end
end

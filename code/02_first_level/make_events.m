function sess = make_events(participant_id, model)
% MAKE_EVENTS  Session-wise conditions, onsets and parametric modulators for one participant.
%
%   sess = make_events('sub-01', 'glm1')
%
%   model: 'glm1' | 'glm2' | 'glm3' | 'heatmap'
%   Reads data/behavior/scan_choices.csv (and prescan_ce.csv for glm2).
%   Sessions are ordered as in the original models: delayed runs (ascending run number), then
%   immediate runs. Each output element has fields run, condition, names, onsets, durations, pmod.
%
%   Trial selection, as in the analyses that produced the paper:
%   - glm2, glm3, heatmap: valid responses to scheduled trials (included_in_glm = 1). Re-presented
%     trials are not modeled. Trials with a sure amount above 20 EUR are kept.
%   - glm1: every valid response (included_in_glm1 = 1); re-presented trials enter with distance = 0
%     and probability = 0, as in the original onset file used for this model.
%   Parametric modulators are entered raw; SPM does not mean-centre them when orthogonalization is off.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
T = readtable(fullfile(root, 'data', 'behavior', 'scan_choices.csv'), 'TextType', 'string');
T = T(T.participant_id == participant_id, :);
if strcmp(model, 'glm1')
    T = T(T.included_in_glm1 == 1, :);
else
    T = T(T.included_in_glm == 1, :);
end

if strcmp(model, 'glm2')   % w(p) = CE / 20 from the pre-scan certainty equivalents
    W = readtable(fullfile(root, 'data', 'behavior', 'prescan_ce.csv'), 'TextType', 'string');
    W = W(W.participant_id == participant_id, :);
end

runs = [unique(T.run(T.condition == "delayed"))', unique(T.run(T.condition == "immediate"))'];
sess = struct('run', {}, 'condition', {}, 'names', {}, 'onsets', {}, 'durations', {}, 'pmod', {});
for k = 1:numel(runs)
    R = sortrows(T(T.run == runs(k), :), 'onset_s');
    cond = char(R.condition(1));
    s.run = runs(k);
    s.condition = cond;
    switch model
        case {'glm1', 'glm2', 'glm3'}
            s.names = {cond};
            s.onsets = {R.onset_s};
            s.durations = {0};
            s.pmod = struct('name', {}, 'param', {}, 'poly', {});
            if strcmp(model, 'glm1')
                dist = R.distance_level; prob = R.probability;
                dist(R.is_replay == 1) = 0; prob(R.is_replay == 1) = 0;
                s.pmod(1) = struct('name', 'distance', 'param', dist, 'poly', 1);
                s.pmod(2) = struct('name', 'prob', 'param', prob, 'poly', 1);
            elseif strcmp(model, 'glm2')
                wp = zeros(height(R), 1);
                for i = 1:height(R)
                    hit = W.condition == R.condition(i) & abs(W.probability - R.probability(i)) < 1e-9;
                    wp(i) = W.w_empirical(hit);
                end
                s.pmod(1) = struct('name', 'distance', 'param', R.distance_level, 'poly', 1);
                s.pmod(2) = struct('name', 'wp', 'param', wp, 'poly', 1);
            end
        case 'heatmap'
            % cell = 3 * (distance bin - 1) + probability bin
            pbin = 1 + (R.probability > 0.10 + 1e-9) + (R.probability > 0.75 + 1e-9);
            dbin = 1 + (R.distance_level > -20) + (R.distance_level > 10);
            cell = 3 * (dbin - 1) + pbin;
            s.names = arrayfun(@(n) sprintf('%s_cell%d', cond, n), 1:9, 'UniformOutput', false);
            s.onsets = arrayfun(@(n) R.onset_s(cell == n), 1:9, 'UniformOutput', false);
            s.durations = repmat({0}, 1, 9);
            s.pmod = struct('name', {}, 'param', {}, 'poly', {});
        otherwise
            error('unknown model %s', model);
    end
    sess(k) = s; %#ok<AGROW>
end
end

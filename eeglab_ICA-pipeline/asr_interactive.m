function [EEG, R] = asr_interactive(EEG_in, saveSet, stem, logMsg, bc_list, zCut)
% ASR_INTERACTIVE  Sweep ASR thresholds, summarize, and let user choose.
%
% Usage:
%   [EEG, R] = asr_interactive(EEG, saveSet, stem, logMsg)
%   [EEG, R] = asr_interactive(EEG, saveSet, stem, logMsg, bc_list, zCut)
%
% Inputs:
%   EEG_in   : EEGLAB dataset before ASR.
%   saveSet  : function handle, e.g. @(desc,EEG) pop_saveset(EEG, ...).
%   stem     : filename stem (e.g., 'sub-01_task-XYZ').
%   logMsg   : logger handle like @(fmt,varargin) fprintf(fid,fmt,varargin{:})
%              (if empty or missing, logging is skipped safely).
%   bc_list  : vector of BurstCriterion values to try (default [30 25 20 15 12]).
%   zCut     : robust z threshold for "heavily changed" samples (default 3).
%
% Outputs:
%   EEG      : dataset after applying the chosen ASR setting.
%   R        : struct array with results for each bc (EEG, repair_frac, etc).

    if nargin < 6 || isempty(zCut),    zCut = 3; end
    if nargin < 5 || isempty(bc_list), bc_list = [30 25 20 15 12]; end
    if nargin < 4, logMsg = []; end

    assert(exist('pop_clean_rawdata','file')==2, ...
        'clean_rawdata/pop_clean_rawdata not found on path.');

    EEG0 = EEG_in;  % reference (pre-ASR)
    X0   = double(EEG0.data);
    % robust per-channel scale (MAD -> ~sigma)
    scale = 1.4826 * mad(X0, 1, 2);
    scale(scale==0) = eps;

    R = struct('bc',[],'EEG',[],'repair_frac',[],'repair_by_chan',[]);
    R(numel(bc_list)).bc = [];  % pre-allocate

    fprintf('\n=== ASR sweep (BurstCriterion) =====================================\n');
    for k = 1:numel(bc_list)
        bc = bc_list(k);
        EEGk = pop_clean_rawdata(EEG0, ...
            'Highpass','off', 'FlatlineCriterion','off', 'LineNoiseCriterion','off', ...
            'ChannelCriterion','off', 'NoisyChannelsCriterion','off', ...
            'BurstCriterion', bc, 'WindowCriterion','off', 'BurstRejection','off', ...
            'PlotFigures','off');

        % repair fraction (per-sample "any channel changed > zCut")
        X1 = double(EEGk.data);
        dz = abs(X1 - X0) ./ scale;             % [ch x time]
        changed_any = any(dz > zCut, 1);        % [1 x time]
        repair_frac = mean(changed_any);        % scalar fraction

        % channel-wise fraction (handy for QA)
        repair_by_chan = mean(dz > zCut, 2);    % [ch x 1]

        R(k).bc            = bc;
        R(k).EEG           = EEGk;
        R(k).repair_frac   = repair_frac;
        R(k).repair_by_chan= repair_by_chan;

        fprintf('  [%d] bc=%-3g  repaired ~%.2f%% of samples (z>%g)\n', ...
            k, bc, 100*repair_frac, zCut);
    end

    % print a tiny legend / hint
    fprintf('---------------------------------------------------------------------\n');
    fprintf('Rule of thumb (ERP): aim ~1–5%% repaired; err on the conservative side.\n');
    fprintf('You can enter the INDEX (e.g., 2) or the bc value itself (e.g., 20).\n');

    % auto-pick suggestion (1–5%% band, else minimum)
    idx_auto = find([R.repair_frac] >= 0.01 & [R.repair_frac] <= 0.05, 1, 'first');
    if isempty(idx_auto), [~, idx_auto] = min([R.repair_frac]); end
    bc_auto  = R(idx_auto).bc;

    prompt = sprintf('Choose index or bc (Enter for bc=%g): ', bc_auto);
    raw = strtrim(input(prompt,'s'));

    % interpret choice
    EEG = R(idx_auto).EEG;   % default
    chosen_idx = idx_auto; chosen_bc = bc_auto;

    if ~isempty(raw)
        num = str2double(raw);
        if ~isnan(num) && isfinite(num)
            % Could be an index OR a bc value
            if any(abs((1:numel(R)) - num) < 1e-12) && round(num)==num
                ii = round(num);
                chosen_idx = max(1, min(numel(R), ii));
                chosen_bc  = R(chosen_idx).bc;
                EEG        = R(chosen_idx).EEG;
            else
                % treat as bc
                [tf,ii] = ismember(num, [R.bc]);
                if tf
                    chosen_idx = ii; chosen_bc = num;
                    EEG        = R(ii).EEG;
                else
                    fprintf(2,'Unrecognized entry. Using bc=%g (index %d).\n', bc_auto, idx_auto);
                end
            end
        else
            fprintf(2,'Invalid entry. Using bc=%g (index %d).\n', bc_auto, idx_auto);
        end
    end

    fprintf('=> Selected bc=%g (index %d), repaired ~%.2f%% of samples.\n', ...
        chosen_bc, chosen_idx, 100*R(chosen_idx).repair_frac);

    % Save checkpoint
    if ~isempty(saveSet) && isa(saveSet,'function_handle')
        tag = sprintf('asrBC%d', round(chosen_bc));
        saveSet(tag, EEG);
    end

    % Log
    if ~isempty(logMsg)
        logMsg('### ASR (burst repair)');
        logMsg('- Candidates tried: %s', strjoin(string(bc_list),' '));
        logMsg('- Heavily-changed zCut: %.2f', zCut);
        for k = 1:numel(R)
            logMsg('- bc=%-3g -> repaired ~%.2f%% of samples', R(k).bc, 100*R(k).repair_frac);
        end
        logMsg('- Selected bc=%g (index %d) -> repaired ~%.2f%%', ...
            chosen_bc, chosen_idx, 100*R(chosen_idx).repair_frac);
        if exist('tag','var')
            logMsg('- Saved as: %s_desc-%s_eeg.set', stem, tag);
        end
    end

    % A little metadata trail
    if ~isfield(EEG,'etc') || ~isstruct(EEG.etc), EEG.etc = struct; end
    EEG.etc.asr = struct( ...
        'bc_tried', bc_list, ...
        'zCut', zCut, ...
        'results_repair_frac', [R.repair_frac], ...
        'chosen_bc', chosen_bc, ...
        'chosen_idx', chosen_idx);

end

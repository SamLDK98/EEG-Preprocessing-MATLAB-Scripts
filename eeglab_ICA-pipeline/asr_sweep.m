function [EEG, T] = bad_channel_triage(EEG, logFile, saveSet, stem, logMsg)

% BAD_CHANNEL_TRIAGE: Interactive clean_rawdata triage + diagnostics + removal UI.

% Fallback logger

if nargin < 5 || isempty(logMsg)
    fidLM = fopen(logFile,'a');
    logMsg = @(varargin) fprintf(fidLM, '%s\n', sprintf(varargin{:}));
    cleanupObj = onCleanup(@() fclose(fidLM)); %#ok<NASGU>
end

assert(exist('pop_clean_rawdata','file')==2, ...
    'clean_rawdata (clean_rawdata/pop_clean_rawdata) not found on path.');

% Freeze a reference copy for diagnostics & viewing
EEG0    = EEG;
labels0 = {EEG0.chanlocs.labels};
fs      = EEG0.srate;

% Settings
thresholds    = [0.85 0.80 0.75];   % channel-corr sweep
hfSD          = 4;                  % High-freq noise SD (clean_rawdata)
flatlineSec   = 5;                  % Dead flatline threshold (s)

% "Dead" rules
lowVarFrac    = 0.05;               % <5% of median channel STD
clipPctThr    = 5;                  % >=5% samples at exact min or max

% "Noisy-but-not-dead" rules / heuristics
hfBand        = [70 120];           % Hz
neigh_k       = 4;                  % nearest neighbors (if 3D positions exist)
neigh_r_poor  = 0.60;               % median corr threshold (poor if below)
hfZ_noisy     = +4;                 % robust z threshold for HF noise
bridge_flag   = 0.10;               % var(diff)/var(ch) << 0.1 suggests bridge/short

% Ocular/slow-drift heuristics
lfBand        = [0.1 2];            % Hz
lfZ_high      = +3;                 % robust z high LF
blink_r_high  = 0.50;               % strong blink coupling

% Dry-run triage (no changes)
flaggedUnion = {};
trialSumm = struct('thr',{},'removed',{});

for k = 1:numel(thresholds)
    thr = thresholds(k);
    EEGt = EEG0;
    EEGt = pop_clean_rawdata(EEGt, ...
        'FlatlineCriterion',      flatlineSec, ...
        'Highpass',               'off', ...
        'LineNoiseCriterion',     'off', ...
        'ChannelCorrelation',     thr, ...
        'NoisyChannelsCriterion', hfSD, ...
        'BurstCriterion',         'off', ...
        'WindowCriterion',        'off', ...
        'BurstRejection',         'off', ...
        'PlotFigures',            'off');

    kept    = {EEGt.chanlocs.labels};
    removed = setdiff(labels0, kept, 'stable');
    trialSumm(end+1) = struct('thr',thr,'removed',{removed}); %#ok<SAGROW>
    flaggedUnion = union(flaggedUnion, removed, 'stable');
end

fprintf('\n=== clean_rawdata triage (dry) ==============================\n');
for k = 1:numel(trialSumm)
    fprintf(' threshold=%.2f -> flagged %d: %s\n', ...
        trialSumm(k).thr, numel(trialSumm(k).removed), joinStrings(trialSumm(k).removed));
end
fprintf(' Union of flagged channels across trials: %d\n', numel(flaggedUnion));

if isempty(flaggedUnion)
    fprintf('   (none flagged — skipping bad-channel step; EEG unchanged.)\n');
    logMsg('### Bad-channel triage');
    logMsg('- Thresholds: %s; HF SD=%g; flatline=%gs', ...
           strjoin(string(thresholds),' '), hfSD, flatlineSec);
    logMsg('- Flagged union: none — skipped removal step.');
    T = table(); % nothing to report
    return;
end

% Diagnostics on flagged set
X = double(EEG0.data);         
nSamp = size(X,2);

% Dead-channel diagnostic
chStd  = std(X, 0, 2);
refStd = median(chStd(chStd > 0));         
flatThreshSmp = round(flatlineSec * fs);

deadResults = table('Size',[numel(flaggedUnion) 7], ...
    'VariableTypes', {'string','logical','double','double','double','double','double'}, ...
    'VariableNames', {'Label','IsDead','FlatRunsGE5','FlatSecs','StdVal','StdFrac','ClipPct'});

for i = 1:numel(flaggedUnion)
    lab = flaggedUnion{i};
    ci  = find(strcmp(labels0, lab), 1);
    x   = X(ci,:);

    dx = diff(x);
    isConst = [false, dx==0];
    d  = diff([0, isConst, 0]);
    runStarts = find(d==1);
    runEnds   = find(d==-1)-1;
    runLens   = runEnds - runStarts + 1;
    longRuns  = runLens >= flatThreshSmp;
    nFlat     = sum(longRuns);
    flatSecs  = sum(runLens(longRuns))/fs;

    s      = std(x);
    sFrac  = s / max(refStd, eps);
    xmin   = min(x); xmax = max(x);
    clipPct= 100*( mean(x==xmin) + mean(x==xmax) );

    isDead = (nFlat > 0) || (sFrac < lowVarFrac) || (clipPct >= clipPctThr);

    deadResults(i,:) = {string(lab), isDead, nFlat, flatSecs, s, sFrac, clipPct};
end

% Noisy-but-not-dead: HF power z, neighbor corr, bridging
win = min(4*fs, floor(nSamp/8)); win = max(win, 2*fs);
noverlap = floor(win/2);
nfft = max(2^nextpow2(win), win);
[pxx,f] = pwelch(X', win, noverlap, nfft, fs); 
pxx = pxx';                                    

hfMask = f >= hfBand(1) & f <= hfBand(2);
hfPow  = trapz(f(hfMask), pxx(:,hfMask), 2);
hfZ    = (hfPow - median(hfPow)) ./ mad(hfPow, 1);

hasPos = all(isfield(EEG0.chanlocs, {'X','Y','Z'}));
if hasPos
    pos = [[EEG0.chanlocs.X]' [EEG0.chanlocs.Y]' [EEG0.chanlocs.Z]'];
else
    pos = [];
end
winSec = 2; stepSec = 1;
winS   = max(1, round(winSec*fs));
stepS  = max(1, round(stepSec*fs));

neighMedR = nan(EEG0.nbchan,1);
for c = 1:EEG0.nbchan
    x = X(c,:);
    if ~isempty(pos)
        d = vecnorm(pos - pos(c,:), 2, 2);
        [~,ord] = sort(d, 'ascend');
        neighIdx = ord(2:min(neigh_k+1, EEG0.nbchan));
    else
        neighIdx = setdiff(1:EEG0.nbchan, c); 
    end
    y = mean(X(neighIdx,:),1);
    rs = [];
    for s = 1:stepS:(nSamp - winS + 1)
        xi = x(s:s+winS-1); yi = y(s:s+winS-1);
        if std(xi)>0 && std(yi)>0
            r = corr(xi', yi');
            if ~isnan(r), rs(end+1) = r; end 
        end
    end
    if ~isempty(rs)
        neighMedR(c) = median(rs);
    end
end

bridgeRatio = nan(EEG0.nbchan,1);
if ~isempty(pos)
    D = squareform(pdist(pos));
    for c = 1:EEG0.nbchan
        [~,ord] = sort(D(c,:),'ascend');
        nb = ord(find(ord~=c,1));       
        diffSig = X(c,:) - X(nb,:);
        bridgeRatio(c) = var(diffSig) / max(var(X(c,:)), eps);
    end
end

noisyResults = table('Size',[numel(flaggedUnion) 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Label','HFz_70_120','NeighMedR','BridgeRatio'});

for i = 1:numel(flaggedUnion)
    ci = find(strcmp(labels0, flaggedUnion{i}),1);
    br = bridgeRatio(ci);
    noisyResults(i,:) = {string(flaggedUnion{i}), hfZ(ci), neighMedR(ci), br};
end

% Ocular / slow drift
lfMask = f >= lfBand(1) & f <= lfBand(2);
lfPow  = trapz(f(lfMask), pxx(:,lfMask), 2);
lfZ    = (lfPow - median(lfPow)) ./ mad(lfPow, 1);

upperL = upper(string(labels0));
eogNames = ["VEOG","VEOG1","VEOG2","HEOG","HEOG1","HEOG2","EOG","FPZ"];
eogIdx = find(ismember(upperL, eogNames));
if isempty(eogIdx)
    pool = find(ismember(upperL, ["FP1","FP2","AFZ"]));
    if isempty(pool)
        pool = find(startsWith(upperL, ["FP","AF"]));
        pool = pool(1:min(3, numel(pool)));
    end
else
    pool = eogIdx;
end
if isempty(pool), pool = 1:min(2, EEG0.nbchan); end
eogSig = mean(X(pool,:),1);

w = max(1, round(0.5*fs));
eogBL = eogSig - movmean(eogSig, w);

blinkMedR = nan(EEG0.nbchan,1);
for c = 1:EEG0.nbchan
    x = X(c,:) - movmean(X(c,:), w);
    rs = [];
    for s = 1:stepS:(nSamp - winS + 1)
        xi = x(s:s+winS-1); yi = eogBL(s:s+winS-1);
        if std(xi)>0 && std(yi)>0
            r = corr(xi', yi'); if ~isnan(r), rs(end+1)=r; end 
        end
    end
    if ~isempty(rs), blinkMedR(c) = median(rs); end
end

ocularResults = table('Size',[numel(flaggedUnion) 3], ...
    'VariableTypes', {'string','double','double'}, ...
    'VariableNames', {'Label','LFz_0_1_2','BlinkMedR'});

for i = 1:numel(flaggedUnion)
    ci = find(strcmp(labels0, flaggedUnion{i}),1);
    ocularResults(i,:) = {string(flaggedUnion{i}), lfZ(ci), blinkMedR(ci)};
end

% Summarize and provide suggestions
T = outerjoin(deadResults, noisyResults, 'Keys','Label','MergeKeys',true);
T = outerjoin(T, ocularResults, 'Keys','Label','MergeKeys',true);

suggest = strings(height(T),1);
for i = 1:height(T)
    if T.IsDead(i)
        suggest(i) = "REMOVE: dead";
    elseif (T.HFz_70_120(i) > hfZ_noisy && T.NeighMedR(i) < neigh_r_poor)
        suggest(i) = "CONSIDER REMOVE: noisy & low corr";
    elseif (~isnan(T.BridgeRatio(i)) && T.BridgeRatio(i) < bridge_flag)
        suggest(i) = "CONSIDER REMOVE: bridging";
    elseif (T.LFz_0_1_2(i) > lfZ_high || T.BlinkMedR(i) > blink_r_high)
        suggest(i) = "KEEP (likely ocular/slow; handle via ICA)";
    else
        suggest(i) = "BORDERLINE — manual review";
    end
end
T.Suggestion = suggest;

fprintf('\n=== Flagged channel evidence ==================================\n');
printTable(T(:, {'Label','IsDead','StdFrac','ClipPct','HFz_70_120', ...
                 'NeighMedR','BridgeRatio','LFz_0_1_2','BlinkMedR','Suggestion'}));

% Log diagnostics
ts = datestr(now,'yyyy-mm-dd HH:MM:SS');
fid = fopen(logFile,'a');
fprintf(fid, "### Bad-channel diagnostics — %s\n", ts);
fprintf(fid, "- Triage params: corr=%s; HF SD=%.1f; flatline=%gs\n", ...
        strjoin(string(thresholds),'/'), hfSD, flatlineSec);
fprintf(fid, "| Label | Dead | STD frac (%%) | Clip (%%) | HF z | Neigh r | Bridge | LF z | Blink r | Suggestion |\n");
fprintf(fid, "|---|:---:|---:|---:|---:|---:|---:|---:|---:|---|\n");
for i = 1:height(T)
    if isnan(T.BridgeRatio(i)), brStr = 'NA'; else, brStr = sprintf('%.2f', T.BridgeRatio(i)); end
    labStr  = char(string(T.Label(i)));
    deadStr = char(string(T.IsDead(i)));
    suggStr = char(string(T.Suggestion(i)));
    fprintf(fid, '| %s | %s | %.1f | %.2f | %.2f | %.2f | %s | %.2f | %.2f | %s |\n', ...
        labStr, deadStr, 100*T.StdFrac(i), T.ClipPct(i), ...
        T.HFz_70_120(i), T.NeighMedR(i), brStr, ...
        T.LFz_0_1_2(i), T.BlinkMedR(i), suggStr);
end
fprintf(fid, '\n');
fclose(fid);

% Open full raw viewer 
viewer_opened = false;
try
    fprintf('\nOpening raw viewer for ALL channels...\n');
    eegplot(EEG0.data, ...
        'srate', EEG0.srate, ...
        'eloc_file', EEG0.chanlocs, ...
        'winlength', 8, ...                     
        'dispchans', min(32, EEG0.nbchan), ...  
        'spacing', 50, ...                      
        'title', sprintf('Raw data: %s', EEG0.setname));
    viewer_opened = true;
    input('Close the viewer window, then press ENTER to continue...','s');
catch ME
    warning('Could not open eegplot: %s', ME.message);
end



% PART 2:
%  - Show the selection menu (or free-form label entry)
%  - Parse numeric choice or labels
%  - Apply removal, restore coordinates, save checkpoint
%  - Final logging, metadata, and helper subfunctions


% Interactive user choice
fprintf('\n=== Selection menu =============================================\n');
fprintf('Options:\n');
fprintf('  [1] Remove only channels suggested as "REMOVE: dead"\n');
fprintf('  [2] Remove "REMOVE: dead" + "CONSIDER REMOVE" suggestions\n');
fprintf('  [3] Manual selection by label (comma- or space-separated)\n');
fprintf('  [0] Do nothing (skip removal)\n');

rawChoice = strtrim(input('Enter choice [1/2/3/0] OR labels (e.g., "AF7 F7"): ','s'));
labels_remove = {};

% Try to interpret as a numeric menu choice
choiceNum = str2double(rawChoice);
isNumericChoice = ~isempty(rawChoice) && ~isnan(choiceNum) && isfinite(choiceNum);

if isNumericChoice
    choice = round(choiceNum);
    switch choice
        case 1
            labels_remove = cellstr(T.Label(T.Suggestion=="REMOVE: dead"))';
        case 2
            labels_remove = cellstr(T.Label(ismember(T.Suggestion, ...
                ["REMOVE: dead","CONSIDER REMOVE: noisy & low corr","CONSIDER REMOVE: bridging"])))';
        case 3
            % Ask for labels now
            raw = strtrim(input('Enter labels to remove (e.g., "Fp1, T7, POz"): ','s'));
            labels_remove = tokenizeLabels(raw);
        otherwise
            % choice == 0 or anything else -> do nothing
            labels_remove = {};
    end
else
    % Treat the entry as a label list directly
    raw = rawChoice;
    labels_remove = tokenizeLabels(raw);
end

% Sanity-check labels against current dataset
if ~isempty(labels_remove)
    bad = setdiff(labels_remove, labels0);
    if ~isempty(bad)
        fprintf(2,'Warning: unknown labels ignored: %s\n', strjoin(bad, ', '));
    end
    labels_remove = intersect(labels_remove, labels0, 'stable');
end

labels_remove = unique(labels_remove, 'stable');
fprintf('\nSelected for removal (%d): %s\n', numel(labels_remove), joinStrings(labels_remove));

if isempty(labels_remove)
    fprintf('No channels removed. EEG unchanged.\n');
    logMsg('### Bad-channel removal (interactive)');
    logMsg('- User entry="%s"; No channels removed.', rawChoice);
else
    % Preserve original coordinates before removal (pop_select may strip fields)
    orig_chanlocs = EEG.chanlocs;
    orig_labels   = {orig_chanlocs.labels};

    % Apply removal
    EEG = pop_select(EEG, 'nochannel', labels_remove);

    % Restore coordinates for remaining channels
    remaining_labels = {EEG.chanlocs.labels};
    coord_fields = {'X','Y','Z','theta','radius','sph_theta','sph_phi','sph_radius'};
    for i = 1:numel(remaining_labels)
        oi = find(strcmp(orig_labels, remaining_labels{i}), 1);
        if ~isempty(oi)
            for j = 1:numel(coord_fields)
                f = coord_fields{j};
                if isfield(orig_chanlocs, f)
                    EEG.chanlocs(i).(f) = orig_chanlocs(oi).(f);
                end
            end
        end
    end

    % Save & log
    tag = sprintf('bcrem_%dch', numel(labels_remove));
    saveSet(tag, EEG);
    logMsg('### Bad-channel removal (interactive)');
    logMsg('- Removed (%d): %s', numel(labels_remove), joinStrings(labels_remove));
    logMsg('- Saved as: %s_desc-%s_eeg.set', stem, tag);

    % Final log block / summary
    try
        fid = fopen(logFile,'a');
        ts  = datestr(now,'yyyy-mm-dd HH:MM:SS');
        fprintf(fid, "### Bad channel removal (final) — %s\n", ts);

        if ~isempty(flaggedUnion)
            fprintf(fid, "- Flagged candidates (union across trials): %s\n", strjoin(flaggedUnion, ', '));
        else
            fprintf(fid, "- Flagged candidates (union across trials): none\n");
        end

        if ~isempty(labels_remove)
            fprintf(fid, "- Removed now (pre-ICA): %s\n", strjoin(labels_remove, ', '));
        else
            fprintf(fid, "- Removed now (pre-ICA): none\n");
        end

        fprintf(fid, "- Total channels after removal: %d\n", EEG.nbchan);
        fprintf(fid, "- Saved as: %s_desc-%s_eeg.set\n\n", stem, tag);
        fclose(fid);
    catch
    end
end

% Dataset metadata (for downstream steps)
if ~isfield(EEG,'etc') || ~isstruct(EEG.etc), EEG.etc = struct; end
if ~isfield(EEG.etc,'preICA') || ~isstruct(EEG.etc.preICA), EEG.etc.preICA = struct; end

EEG.etc.preICA.flagged_union = flaggedUnion(:)';     
if exist('labels_remove','var'), EEG.etc.preICA.removed_preICA = labels_remove(:)'; else, EEG.etc.preICA.removed_preICA = {}; end
EEG.etc.preICA.decision_time = datestr(now,'yyyy-mm-dd HH:MM:SS');
EEG.etc.preICA.triage_params = struct( ...
    'corr_thresholds', thresholds, ...
    'hfSD', hfSD, 'flatlineSec', flatlineSec, ...
    'hfBand', hfBand, 'lfBand', lfBand, ...
    'neigh_k', neigh_k, 'neigh_r_poor', neigh_r_poor, ...
    'hfZ_noisy', hfZ_noisy, 'bridge_flag', bridge_flag, ...
    'lfZ_high', lfZ_high, 'blink_r_high', blink_r_high);

% Also keep original triage table in the dataset for traceability
EEG.etc.ica_prep.badch = struct( ...
    'triage_thresholds', thresholds, ...
    'hf_sd', hfSD, 'flatline_s', flatlineSec, ...
    'dead_rules', struct('lowVarFrac',lowVarFrac,'clipPctThr',clipPctThr), ...
    'noisy_rules', struct('hfBand',hfBand,'hfZ_noisy',hfZ_noisy, ...
                          'neigh_k',neigh_k,'neigh_r_poor',neigh_r_poor, ...
                          'bridge_flag',bridge_flag), ...
    'ocular_rules', struct('lfBand',lfBand,'lfZ_high',lfZ_high,'blink_r_high',blink_r_high), ...
    'flagged_union',{flaggedUnion}, ...
    'table', T);

end % ====================== end main function ==============================


% ====================== Local helper functions ============================
function s = joinStrings(C)
    if isempty(C), s = '(none)'; return; end
    if isstring(C), C = cellstr(C); end
    s = strjoin(C, ', ');
end

function out = tokenizeLabels(raw)
    if isempty(raw), out = {}; return; end
    raw = strrep(raw, ',', ' ');
    parts = strsplit(strtrim(raw));
    out = cellfun(@strtrim, parts, 'UniformOutput', false);
    out(cellfun(@isempty,out)) = [];
end

function printTable(T)
% CLI-friendly table printer (robust to strings/NaNs)
    varNames = string(T.Properties.VariableNames);
    C = cell(height(T), numel(varNames));
    for j = 1:numel(varNames)
        for i = 1:height(T)
            val = T{i,j};
            if ismissing(val) || (isfloat(val) && isnan(val))
                C{i,j} = "NA";
            elseif islogical(val)
                C{i,j} = string(val);
            elseif isstring(val) || ischar(val)
                C{i,j} = string(val);
            elseif isfloat(val)
                name = varNames(j);
                x = val;
                if contains(name, "StdFrac")
                    C{i,j} = sprintf('%.1f', 100*x);
                elseif contains(name, "ClipPct")
                    C{i,j} = sprintf('%.2f', x);
                elseif contains(name, ["HFz","LFz"])
                    C{i,j} = sprintf('%.2f', x);
                elseif contains(name, ["NeighMedR","BlinkMedR","BridgeRatio"])
                    C{i,j} = sprintf('%.2f', x);
                else
                    C{i,j} = sprintf('%.3g', x);
                end
                C{i,j} = string(C{i,j});
            else
                try, C{i,j} = string(val); catch, C{i,j} = "<obj>"; end
            end
        end
    end
    W = strlength(varNames);
    for j = 1:numel(varNames)
        W(j) = max(W(j), max(strlength(string(C(:,j)))));
        W(j) = max(W(j), 6);
    end
    for j = 1:numel(varNames)
        fmt = sprintf('%%-%ds  ', W(j)); fprintf(fmt, varNames(j));
    end
    fprintf('\n');
    for j = 1:numel(varNames), fprintf('%s  ', repmat('-',1, W(j))); end
    fprintf('\n');
    for i = 1:height(T)
        for j = 1:numel(varNames)
            fmt = sprintf('%%-%ds  ', W(j));
            fprintf(fmt, C{i,j});
        end
        fprintf('\n');
    end
end

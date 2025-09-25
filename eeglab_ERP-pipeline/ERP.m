%% ========================================================================
%%                           ERP branch
%% ========================================================================

% User inputs
sub  = strtrim(input('Enter subject ID (e.g., sub-01): ','s'));
task = strtrim(input('Enter task name (e.g., CogAss1): ','s'));
stem = sprintf('%s_task-%s', sub, task);

% ERP derivatives root
erpDerivRoot = 'C:\Users\samki\Documents\EEG-datasets\krigolson-data\derivatives\eeglab_ERP-pipeline';
derivName    = 'eeglab_ERP-pipeline';   % for log header

% Derive the ICA derivatives root from the ERP root's parent \derivatives\
[derivParent, ~, ~] = fileparts(erpDerivRoot);                 % ...\derivatives
icaDerivRoot = fullfile(derivParent, 'eeglab-ICA_pipeline');   % ICA branch name

% Paths
srcDatasetDir = fullfile(icaDerivRoot, sub, 'eeg');   % read from here
erpDatasetDir = fullfile(erpDerivRoot, sub, 'eeg');   % save to here
figsdir_erp   = fullfile(erpDatasetDir, 'figs');
logdir_erp    = fullfile(erpDatasetDir, 'log');

if ~exist(erpDatasetDir,'dir'), mkdir(erpDatasetDir); end
if ~exist(figsdir_erp,'dir'),   mkdir(figsdir_erp);   end
if ~exist(logdir_erp,'dir'),    mkdir(logdir_erp);    end

% ERP log
logPath_erp = fullfile(logdir_erp, 'preprocessing_log.md');

if ~exist(logPath_erp,'file')
    fid0 = fopen(logPath_erp,'w');
    ts0  = datestr(now,'yyyy-mm-dd HH:MM:SS');
    fprintf(fid0, "# ERP Preprocessing Log — %s (%s)\n\n", sub, ts0);
    fprintf(fid0, "- Derivatives branch: %s\n", derivName);
    fprintf(fid0, "- Task: %s\n\n", task);
    fclose(fid0);
end

% Open in append mode and make a helper
fid = fopen(logPath_erp, 'a');
assert(fid > 0, 'Failed to open ERP log for appending: %s', logPath_erp);
logMsg = @(fmt,varargin) fprintf(fid, [fmt '\n'], varargin{:});
logHdr = @(hdr) fprintf(fid, '\n### %s\n', hdr);

tsStart = datestr(now,'yyyy-mm-dd HH:MM:SS');
logHdr(sprintf('Attach ICA — start (%s)', tsStart));
logMsg('- Subject: %s', sub);
logMsg('- Task: %s', task);
logMsg('- Source (ICA) dir: %s', srcDatasetDir);
logMsg('- Dest (ERP) dir: %s', erpDatasetDir);

% Non-EEG labels for avg-ref exclusion
nonEEG_labels = upper(string( ...
    {'HEOG','VEOG','EOG','HEOG1','HEOG2','VEOG1','VEOG2', ...
     'EMG','EMG1','EMG2','ECG','M1','M2','A1','A2'}));

try

    % --------------------------------------------------------------------
    % 1) Load ERP seed (desc-linrm50/60) from ICA branch
    % --------------------------------------------------------------------
    candidates = { ...
        sprintf('%s_desc-linrm50_eeg.set', stem), ...
        sprintf('%s_desc-linrm60_eeg.set', stem) ...
    };
    seedFile = '';
    for i = 1:numel(candidates)
        if exist(fullfile(srcDatasetDir, candidates{i}), 'file')
            seedFile = candidates{i};
            break;
        end
    end
    assert(~isempty(seedFile), 'No desc-linrmXX seed found in %s for %s.', srcDatasetDir, stem);

    EEG = pop_loadset('filename', seedFile, 'filepath', srcDatasetDir);
    logHdr('Load seed');
    logMsg('- Loaded: %s', fullfile(srcDatasetDir, seedFile));
    logMsg('- srate: %.3f Hz; chans: %d; range: [%.3f %.3f] s', EEG.srate, EEG.nbchan, EEG.xmin, EEG.xmax);

    % --------------------------------------------------------------------
    % 2) Load ICA owner (frozen weights dataset) 
    % --------------------------------------------------------------------
    ownerFile = sprintf('%s_desc-icaMODEL_eeg.set', stem);
    assert(exist(fullfile(srcDatasetDir, ownerFile),'file')==2, ...
        'ICA owner dataset not found: %s', fullfile(srcDatasetDir, ownerFile));
    EEG_owner = pop_loadset('filename', ownerFile, 'filepath', srcDatasetDir);

    % Summarize owner
    ownerRefStr = '<unknown>';
    if isfield(EEG_owner,'ref') && ~isempty(EEG_owner.ref)
        if ischar(EEG_owner.ref) || isstring(EEG_owner.ref)
            ownerRefStr = char(EEG_owner.ref);
        elseif isnumeric(EEG_owner.ref)
            ownerRefStr = 'numeric';
        end
    end
    logHdr('Load ICA owner');
    logMsg('- Loaded: %s', fullfile(srcDatasetDir, ownerFile));
    logMsg('- srate: %.3f Hz; chans: %d; ref: %s', EEG_owner.srate, EEG_owner.nbchan, ownerRefStr);

    % --------------------------------------------------------------------
    % 3) Match channels: drop extras; reorder to match owner 
    % --------------------------------------------------------------------
    labs_erp   = upper(string({EEG.chanlocs.labels}));
    labs_owner = upper(string({EEG_owner.chanlocs.labels}));

    missing = setdiff(labs_owner, labs_erp);
    assert(isempty(missing), 'ERP seed is missing owner channels: %s', strjoin(missing, ', '));

    extra = setdiff(labs_erp, labs_owner);
    if ~isempty(extra)
        logMsg('- Dropping extra channels: %s', strjoin(extra, ', '));
        keepMask = ismember(labs_erp, labs_owner);
        EEG = pop_select(EEG, 'channel', find(keepMask));
        labs_erp = upper(string({EEG.chanlocs.labels}));
    end

    [~, idxERP] = ismember(labs_owner, labs_erp);
    assert(all(idxERP>0), 'Failed to align ERP channels to owner order.');
    EEG = pop_select(EEG, 'channel', idxERP);

    assert(isequal(upper(string({EEG.chanlocs.labels})), labs_owner), ...
        'Channel order mismatch after reordering.');

    logHdr('Channel alignment');
    logMsg('- Matched labels/order to owner (%d channels).', EEG.nbchan);

    % --------------------------------------------------------------------
    % 4) Apply the same reference as owner 
    % --------------------------------------------------------------------
    if isfield(EEG_owner,'ref') && ~isempty(EEG_owner.ref) && ...
            (ischar(EEG_owner.ref) || isstring(EEG_owner.ref)) && strcmpi(EEG_owner.ref,'average')

        labsNow   = upper(string({EEG.chanlocs.labels}));
        excludeIx = find(ismember(labsNow, nonEEG_labels));
        EEG = pop_reref(EEG, [], 'exclude', excludeIx); % average ref with same excludes
        logHdr('Reference');
        exclList = '(none)'; if ~isempty(excludeIx), exclList = strjoin(labsNow(excludeIx), ', '); end
        logMsg('- Applied: average reference (exclude: %s)', exclList);

    elseif isfield(EEG_owner,'ref') && ~isempty(EEG_owner.ref) && isnumeric(EEG_owner.ref)

        EEG = pop_reref(EEG, EEG_owner.ref);           % replicate numeric ref
        logHdr('Reference');
        logMsg('- Applied: numeric reference (owner.ref size: %s)', mat2str(size(EEG_owner.ref)));

    else
        % Fallback: average with excludes (warn + log)
        labsNow   = upper(string({EEG.chanlocs.labels}));
        excludeIx = find(ismember(labsNow, nonEEG_labels));
        EEG = pop_reref(EEG, [], 'exclude', excludeIx);
        logHdr('Reference');
        exclList = '(none)'; if ~isempty(excludeIx), exclList = strjoin(labsNow(excludeIx), ', '); end
        logMsg('- Owner ref unknown — defaulted to average (exclude: %s)', exclList);
    end

    % --------------------------------------------------------------------
    % 5) Compatibility checks: chans/order/count/srate/reference descriptor
    % --------------------------------------------------------------------
    refA = '<unknown>'; refB = '<unknown>';
    if isfield(EEG,'ref') && ~isempty(EEG.ref)
        if ischar(EEG.ref) || isstring(EEG.ref), refA = char(EEG.ref);
        elseif isnumeric(EEG.ref), refA = 'numeric'; else, refA = class(EEG.ref); end
    end
    if isfield(EEG_owner,'ref') && ~isempty(EEG_owner.ref)
        if ischar(EEG_owner.ref) || isstring(EEG_owner.ref), refB = char(EEG_owner.ref);
        elseif isnumeric(EEG_owner.ref), refB = 'numeric'; else, refB = class(EEG_owner.ref); end
    end

    assert(EEG.nbchan == EEG_owner.nbchan, ...
        'Channel count differs (ERP %d vs owner %d).', EEG.nbchan, EEG_owner.nbchan);
    assert(abs(EEG.srate - EEG_owner.srate) < 1e-6, ...
        'Sampling rate differs (ERP %.6f vs owner %.6f).', EEG.srate, EEG_owner.srate);
    assert(strcmpi(refA, refB), 'Reference differs (ERP: %s, owner: %s).', refA, refB);
    assert(isequal(upper(string({EEG.chanlocs.labels})), labs_owner), ...
        'Final label/order mismatch before attaching weights.');

    logHdr('Compatibility check');
    logMsg('- OK: chans=%d, srate=%.3f Hz, ref=%s', EEG.nbchan, EEG.srate, refA);

    % --------------------------------------------------------------------
    % 6) Attach ICA weights
    % --------------------------------------------------------------------
    assert(isfield(EEG_owner,'icaweights') && ~isempty(EEG_owner.icaweights), 'Owner: icaweights missing.');
    assert(isfield(EEG_owner,'icasphere')  && ~isempty(EEG_owner.icasphere),  'Owner: icasphere missing.');

    EEG.icaweights  = EEG_owner.icaweights;
    EEG.icasphere   = EEG_owner.icasphere;
    EEG.icachansind = EEG_owner.icachansind;
    EEG = eeg_checkset(EEG);

    logHdr('Attach ICA');
    logMsg('- Weights: %dx%d (W), sphere: %dx%d', size(EEG.icaweights), size(EEG.icasphere));

    % --------------------------------------------------------------------
    % 7) Save checkpoint to ERP derivatives folder 
    % --------------------------------------------------------------------
    outName = sprintf('%s_desc-erpAttach_eeg.set', stem);
    EEG = pop_saveset(EEG, 'filename', outName, 'filepath', erpDatasetDir);

    logHdr('Save checkpoint');
    logMsg('- Saved: %s', fullfile(erpDatasetDir, outName));

    tsEnd = datestr(now,'yyyy-mm-dd HH:MM:SS');
    logHdr(sprintf('Attach ICA — end (%s)', tsEnd));
    logMsg('---');

catch ME
    % Log the error and rethrow for visibility
    logHdr('ERROR');
    logMsg('- %s: %s', ME.identifier, ME.message);
    fclose(fid);
    rethrow(ME);
end

% Close the log
fclose(fid);


%% choose spectrum range (up to 120 Hz, but <= Nyquist) 
nyq   = EEG.srate/2;
fMax  = min(120, floor(nyq) - 1);            
if fMax < 60, fMax = max(10, floor(nyq)-1); end
freqRange = [1 fMax];

%% ------------------------------------------------------------------------
%% Review only the "To Consider" ICs (sorted by suspicion) 
%% ------------------------------------------------------------------------

if isempty(idxConsider)
    fprintf('No ICs in "To Consider" — nothing to review.\n');
    if exist('logPath_erp','var')
        fid = fopen(logPath_erp, 'a'); if fid>0, logHdr('IC review'); logMsg('- No ICs to consider.'); fclose(fid); end
    end
else
    % Sort "To Consider" by max non-brain (excl Other), descending
    [scoreAll, winnerAll] = max(nonBrainNoOther, [], 2);
    scoreConsider  = scoreAll(idxConsider);
    [~, srt] = sort(scoreConsider, 'descend');
    idxConsiderSorted  = idxConsider(srt);

    % Ensure icawinv exists and is the right size
    if ~isfield(EEG,'icawinv') || isempty(EEG.icawinv) || size(EEG.icawinv,1) ~= EEG.nbchan
        EEG = eeg_checkset(EEG, 'ica');  % compute icawinv from weights/sphere if needed
    end

    % Drop any ICs with invalid maps 
    validIC = find(all(isfinite(EEG.icawinv),1));
    idxConsiderValid = intersect(idxConsiderSorted, validIC);
    dropped = setdiff(idxConsiderSorted, idxConsiderValid);

    fprintf('\n"To Consider" (sorted by max non-brain prob): %d ICs\n', numel(idxConsiderSorted));
    if ~isempty(dropped)
        fprintf('  Skipping %d IC(s) with invalid maps: %s\n', numel(dropped), strtrim(sprintf('%d ', dropped)));
    end

    if exist('logPath_erp','var')
        fid = fopen(logPath_erp, 'a');
        if fid>0
            logHdr('IC review');
            logMsg('- Opening properties for "To Consider" ICs only (n=%d), freqrange=[%d %d] Hz', ...
                   numel(idxConsiderValid), freqRange(1), freqRange(2));
            if ~isempty(dropped), logMsg('- Skipped invalid maps: %s', strtrim(sprintf('%d ', dropped))); end
            fclose(fid);
        end
    end

    % Try multi-IC viewer; fall back to one-by-one if needed
    if exist('pop_viewprops','file')==2
        try
            fprintf('Opening pop_viewprops for %d "To Consider" IC(s) with freqrange [%g %g] Hz...\n', ...
                    numel(idxConsiderValid), freqRange(1), freqRange(2));
            pop_viewprops(EEG, 0, idxConsiderValid, {'freqrange', freqRange});  
        catch ME
            fprintf('pop_viewprops failed (%s). Falling back to pop_prop...\n', ME.message);
            for k = idxConsiderValid(:)'
                try
                    pop_prop(EEG, 0, k, NaN, {'freqrange', freqRange});
                catch
                    fprintf('  pop_prop failed for IC %d (skipping)\n', k);
                end
            end
        end
    else
        fprintf('pop_viewprops not found; opening pop_prop for each "To Consider" IC with freqrange [%g %g] Hz...\n', freqRange(1), freqRange(2));
        for k = idxConsiderValid(:)'
            try
                pop_prop(EEG, 0, k, NaN, {'freqrange', freqRange});
            catch
                fprintf('  pop_prop failed for IC %d (skipping)\n', k);
            end
        end
    end
end

%% ------------------------------------------------------------------------
%% Open a scrollable time-series viewer of all IC activations
%% ------------------------------------------------------------------------
try
    fprintf('Opening IC activations (time series) for ALL ICs...\n');
    pop_eegplot(EEG, 0, 1, 0); 
catch
    try
        act = eeg_getica(EEG);  % computes ICA activations if needed
        eegplot(act, 'srate', EEG.srate, ...
                     'title', sprintf('%s — IC activations (all)', stem), ...
                     'limits', [EEG.xmin EEG.xmax]*1000);
    catch
        fprintf('Could not open IC time series viewer (both methods failed).\n');
    end
end

%% ------------------------------------------------------------------------
%% Topomap grid with ICLabel percentages
%% ------------------------------------------------------------------------
if ~exist('idxConsiderValid','var') || isempty(idxConsiderValid)
    showSet = 1:size(P,1);     % fall back to all
    gridName = sprintf('%s — ICLabel topo grid (ALL ICs)', stem);
else
    showSet = idxConsiderValid; % show Consider ICs
    gridName = sprintf('%s — ICLabel topo grid (To Consider, n=%d)', stem, numel(showSet));
end

% Create compact grid
nShow = numel(showSet);
ncol  = min(10, nShow);
nrow  = ceil(nShow/ncol);
figure('Name', gridName, 'Color','w');
for ii = 1:nShow
    c = showSet(ii);
    subplot(nrow, ncol, ii);
    try
        topoplot(EEG.icawinv(:,c), EEG.chanlocs, 'electrodes','on', 'numcontour',5, ...
                 'emarker',{'.','k',7,1});
    catch
        axis off
        title(sprintf('IC%d (map err)', c), 'FontSize',8);
        continue
    end
    pr = 100 * P(c,:);  
    ttl = sprintf('IC%d  Br%2.0f Ey%2.0f Mu%2.0f He%2.0f Ln%2.0f Ch%2.0f O%2.0f', ...
                  c, pr(iBrain), pr(iEye), pr(iMuscle), pr(iHeart), pr(iLine), pr(iChNoise), pr(iOther));
    title(ttl, 'FontSize',8);
end

% Log these extra viewers
if exist('logPath_erp','var')
    fid = fopen(logPath_erp, 'a');
    if fid>0
        logHdr('Extra reviewers');
        logMsg('- Opened IC activations time-series (all ICs).');
        logMsg('- Opened IC property viewers with freqrange=[%d %d] Hz.', freqRange(1), freqRange(2));
        logMsg('- Opened ICLabel topo grid: %s', gridName);
        fclose(fid);
    end
end

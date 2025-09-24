%% ========================================================================
%%                            session setup
%% ========================================================================

% Self-rooting
st = dbstack('-completenames');
thisFile = st(1).file;
codeDir  = fileparts(thisFile);      % .../project-name/code
root     = fileparts(codeDir);       % .../project-name  (project root)

% pin pipeline code for this project only
addpath(fullfile(root,'eeglab_ICA-pipeline')); 

% Per-run inputs
sub  = strtrim(input('Enter subject ID (e.g., sub-01): ','s'));
task = strtrim(input('Enter task name (e.g., CogAss1): ','s'));

% Project-level settings
pipeline = 'eeglab-ICA_pipeline';

% Paths
datasetDir = fullfile(root, 'derivatives', pipeline, sub, 'eeg');
figsDir    = fullfile(datasetDir, 'figs');
logDir     = fullfile(datasetDir, 'log');
logFile    = fullfile(logDir, 'preprocessing_log.md');

mkdir(datasetDir);
mkdir(figsDir);
mkdir(logDir);

% Filenames/stems
stem      = sprintf('%s_task-%s', sub, task);
startfile = sprintf('%s_desc-rawimport_eeg.set', stem);  % output from your import/audit stage

% Load starting dataset
startpath = fullfile(datasetDir, startfile);
EEG = pop_loadset('filename', startfile, 'filepath', datasetDir);

% Save chanlocs coords
MASTER_CHANLOCS = EEG.chanlocs;

% Saver for checkpoints
saveSet = @(desc,EEGvar) pop_saveset(EEGvar, ...
    'filename', sprintf('%s_desc-%s_eeg.set', stem, desc), ...
    'filepath', datasetDir);
saveFig = @(suffix) saveas(gcf, fullfile(figsDir, sprintf('%s_%s.png', stem, suffix)));

% log
fid = fopen(logFile, 'a');
fprintf(fid, '\n\n');
fprintf(fid, '# Preprocessing for ERP: ICA Preparation for %s_task-%s\n', sub, task);
fprintf(fid, 'Started: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
logMsg = @(varargin) ...
    fprintf(fid, '%s%s\n', ...
        startsWith(string(varargin{1}), '#') * newline, ... % add \n if header
        sprintf(varargin{:}));

%% ========================================================================
%%                           ica workflow
%% ========================================================================

%% 1. Record original montage
EEG.etc.orig_chanlocs = EEG.chanlocs;
logMsg('- Recorded original montage (saved in EEG.etc.orig_chanlocs)');


%% 2. Trim settling period - enter 1st stimulus index!

ev   = EEG.event;
tSec = ([ev.latency]-1)/EEG.srate + EEG.xmin;
k    = min(10, numel(ev));
fprintf('\nFirst %d events:\n', k);
for i=1:k, fprintf('  %3d : %s @ %.3f\n', i, string(ev(i).type), tSec(i)); end

idx = input('Index of first real stimulus (counting from 1): ');
assert(isscalar(idx) && idx==round(idx) && idx>=1 && idx<=numel(ev), 'Bad index.');

pre    = 20;
first  = tSec(idx);
orig   = [EEG.xmin EEG.xmax];
start  = first - pre;
didTrim = start > 0 + 1/EEG.srate;

if didTrim
    EEG = pop_select(EEG, 'time', [max(start, EEG.xmin) EEG.xmax]);
    if exist('saveSet','var'), saveSet('trim20s', EEG); end
end

EEG.etc.ica_prep.trim = struct( ...
    'policy', sprintf('%ds pre-stim', pre), ...
    'idx', idx, 'type', string(ev(idx).type), ...
    'first_stim_s', first, 'applied', didTrim, ...
    'orig_range_s', orig, 'post_range_s', [EEG.xmin EEG.xmax]);

% Log

if exist('logMsg','var')
    logMsg('### Trim (manual)');
    logMsg('- First stim: idx=%d, type=%s, t=%.3f s', idx, string(ev(idx).type), first);
    if didTrim
        logMsg('- Applied: start=%.3f s (kept ~%d s pre-stim)', start, pre);
        logMsg('- Saved: %s_desc-trim20s_eeg.set', sprintf('%s_task-%s', sub, task));
    else
        logMsg('- Skipped: start=%.3f s <= 0 (assumed pre-trimmed)', start);
    end
end

%% 3. Downsampling

origSrate = EEG.srate;
fprintf('Current sampling rate: %.3f Hz\n', origSrate);
ansHz = input('New sampling rate in Hz (Enter to skip): ');

if isempty(ansHz)
    target = origSrate;   % skip
else
    target = round(ansHz);  % enforce new Hz
end

if abs(origSrate - target) <= 1e-6
    if exist('logMsg','var')
        logMsg('### Downsample');
        logMsg('- Original sampling rate already %.3f Hz → step skipped', origSrate);
    end
else
    EEG = pop_resample(EEG, target);     % anti-aliasing included
    tag = sprintf('ds%d', target);
    if exist('saveSet','var'), saveSet(tag, EEG); end

    if exist('logMsg','var')
        logMsg('### Downsample');
        logMsg('- Original sampling rate: %.3f Hz', origSrate);
        logMsg('- New sampling rate: %.3f Hz', EEG.srate);
        logMsg('- Saved as: %s_desc-%s_eeg.set', stem, tag);
    end
end

EEG.etc.ica_prep.downsample = struct( ...
    'applied', abs(origSrate-target) > 1e-6, ...
    'orig_hz', origSrate, 'target_hz', target, 'new_hz', EEG.srate);


%% 4. Line noise removal

fprintf('Line noise removal. Current sampling rate: %.3f Hz\n', EEG.srate);
f0 = input('Mains frequency [50 or 60]: ');
assert(ismember(f0,[50 60]), 'Enter 50 or 60.');
assert(exist('pop_cleanline','file')==2, 'CleanLine plugin not found on the path.');

nyq   = EEG.srate/2;
harms = f0:f0:nyq;

EEG = pop_cleanline(EEG, 'linefreqs', harms, 'chanlist', 1:EEG.nbchan);

tag = sprintf('linrm%d', f0);
if exist('saveSet','var'), saveSet(tag, EEG); end
if exist('logMsg','var')
    logMsg('### Line noise');
    logMsg('- Base: %d Hz; harmonics: %s', f0, strjoin(string(harms), ' '));
    logMsg('- Saved as: %s_desc-%s_eeg.set', stem, tag);
end

EEG.etc.ica_prep.line_noise = struct('method','cleanline','base_hz',f0,'harmonics',harms);


%% 5. High-pass - 1 Hz

orig = [EEG.xmin EEG.xmax];
EEG  = pop_eegfiltnew(EEG, 1, []); % zero-phase FIR HP at 1 Hz

if exist('saveSet','var'), saveSet('hp1', EEG); end
if exist('logMsg','var')
    logMsg('### High-pass (ICA branch)');
    logMsg('- Applied: 1 Hz high-pass (pop_eegfiltnew)');
    logMsg('- Saved as: %s_desc-hp1_eeg.set', stem);
end
EEG.etc.ica_prep.highpass = struct('cutoff_hz',1,'orig_range_s',orig,'post_range_s',[EEG.xmin EEG.xmax]);


%% 6. Gross artifact cleanup - MANUAL CHECK NEEDED

% a. channel pruning 
[EEG, removedChans] = channel_prune_workbench(EEG, saveSet, logMsg, stem);

% b. segment pruning 
[EEG, newBounds] = manual_segment_prune(EEG, saveSet, logMsg, stem);

%% 7. Bad channel removal - clean_rawdata

[EEG, T_bct] = bad_channel_triage(EEG, logFile, saveSet, stem, logMsg);

%% 8. ASR

[EEG, R] = asr_interactive(EEG, saveSet, stem, logMsg);


%% 9. Average reference

% Identify channels to EXCLUDE from the average (don't let these bias the ref)
nonEEG_labels = upper(string( ...
    {'HEOG','VEOG','EOG','HEOG1','HEOG2','VEOG1','VEOG2', ...
     'EMG','EMG1','EMG2','ECG','M1','M2','A1','A2'}));  % add your lab's tags

allLabs  = upper(string({EEG.chanlocs.labels}));
exclude  = find(ismember(allLabs, nonEEG_labels));

% Keep a copy of chanlocs (paranoia; reref should not touch them)
chanlocs_before = EEG.chanlocs;

% Apply average reference over EEG-only channels
EEG = pop_reref(EEG, [], 'exclude', exclude);  % avg-ref, exclude listed chans

% Restore chanlocs just in case (no-op if unchanged)
EEG.chanlocs = chanlocs_before;

% Rank estimate after average-ref
if exist('eeg_rank','file') == 2
    pcaRank = eeg_rank(double(EEG.data));
else
    pcaRank = min(rank(double(EEG.data)'), EEG.nbchan - 1);
end

% Save checkpoint
saveSet('avgRef', EEG);

% Log
if exist('logMsg','var')
    if ~isempty(exclude)
        excl_list = strjoin(string({chanlocs_before(exclude).labels}), ', ');
    else
        excl_list = '(none)';
    end
    logMsg('### Re-reference');
    logMsg('- Reference: average (excluding: %s)', excl_list);
    logMsg('- Estimated rank after reref: %d', pcaRank);
    logMsg('- Saved as: %s_desc-avgRef_eeg.set', stem);
end

fprintf('Avg-ref done. Excluded %d channel(s). Estimated rank=%d\n', numel(exclude), pcaRank);


%% 10. Run ICA 

% Warn if not avg-ref
if ~isfield(EEG,'ref') || ~(ischar(EEG.ref) && strcmpi(EEG.ref,'average'))
    warning('EEG.ref does not indicate average reference. ICA is typically run after avg-ref.');
end

% Ensure pcaRank exists
if ~exist('pcaRank','var') || isempty(pcaRank) || ~isscalar(pcaRank)
    try
        dataRank = rank(double(EEG.data'));
    catch
        dataRank = EEG.nbchan;
    end
    pcaRank = min([EEG.nbchan-1, dataRank]);
    fprintf('pcaRank not provided — using %d.\n', pcaRank);
end

fprintf('Running ICA: runica (extended=1), PCA=%d ...\n', pcaRank);
tic;
EEG = pop_runica(EEG, 'icatype','runica', 'extended',1, 'pca', pcaRank, 'interrupt','on');
t_ica = toc;
fprintf('ICA finished in %.1f min.\n', t_ica/60);

% Save an ICA checkpoint (weights embedded)
saveSet('ica', EEG);

% Log
if exist('logMsg','var')
    logMsg('### ICA');
    logMsg('- Algorithm: runica (extended)');
    logMsg('- PCA dim: %d', pcaRank);
    logMsg('- Duration: %.1f min', t_ica/60);
    logMsg('- Saved as: %s_desc-ica_eeg.set', stem);
end

%% =======================================================================
%%            save ICA weights (frozen owner) for transfer later
%% =======================================================================

% Safety checks
assert(isfield(EEG,'icaweights') && ~isempty(EEG.icaweights), 'ICA weights missing.');
assert(isfield(EEG,'icasphere')  && ~isempty(EEG.icasphere),  'ICA sphere missing.');
if ~isfield(EEG,'icachansind') || isempty(EEG.icachansind)
    EEG.icachansind = 1:EEG.nbchan;
end

% Save a canonical "ICA owner" dataset
saveSet('icaMODEL', EEG);

% Also save core fields to a .mat for later attachment
icaOut = struct();
icaOut.icaweights   = EEG.icaweights;
icaOut.icasphere    = EEG.icasphere;
icaOut.icachansind  = EEG.icachansind;
icaOut.labels       = {EEG.chanlocs.labels};
icaOut.nbchan       = EEG.nbchan;
icaOut.srate        = EEG.srate;

% Safe 'ref' extraction (no getfield default!)
refStr = 'unknown';
if isfield(EEG,'ref') && ~isempty(EEG.ref)
    if ischar(EEG.ref) || isstring(EEG.ref)
        refStr = char(EEG.ref);
    elseif isnumeric(EEG.ref)
        refStr = mat2str(EEG.ref);
    else
        refStr = class(EEG.ref);
    end
end
icaOut.ref          = refStr;

icaOut.pcaDim       = size(EEG.icaweights,1);
icaOut.setname      = EEG.setname;

weightsMat = fullfile(datasetDir, sprintf('%s_desc-icaWeights.mat', stem));
save(weightsMat, '-struct', 'icaOut');

% Log
if exist('logMsg','var')
    logMsg('### ICA model save');
    logMsg('- Purpose: frozen owner of ICA weights (for transfer later)');
    logMsg('- Saved dataset: %s_desc-icaMODEL_eeg.set', stem);
    logMsg('- Saved fields: %s', strjoin(fieldnames(icaOut)', ', '));
    logMsg('- Training conditions: HP ~1 Hz, avg-ref, bad channels removed pre-ICA');
    logMsg('- ICA algorithm: runica (extended)');
    logMsg('- PCA rank used: %d', pcaRank);
end

%% =======================================================================
%%        create second ERP branch (derive paths from current file)
%% =======================================================================

[eegdir_erp, stem_erp, saveStepERP, logPath_erp] = create_erp_branch(EEG, ...
    'root', root, 'stem', stem);


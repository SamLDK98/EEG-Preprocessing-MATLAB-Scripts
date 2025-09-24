%% ========================================================================
%%                               raw import
%% ========================================================================

% krigolson-dataset

% Per-run inputs (robust to "01" vs "sub-01")
sub_in = strtrim(input('Enter subject ID (e.g., sub-01 OR 01): ','s'));
if isempty(sub_in), error('Subject ID cannot be empty.'); end
if startsWith(sub_in, 'sub-')
    sub = sub_in;
else
    num = regexprep(sub_in, '^\D+', ''); % strip any leading non-digits
    if isempty(num), error('Subject ID must be like "sub-01" or a number like "01".'); end
    sub = sprintf('sub-%s', num);
end

task = strtrim(input('Enter task name (e.g., CogAss1): ','s'));
if isempty(task), error('Task name cannot be empty.'); end

% Project-level settings
pipeline = 'eeglab_ICA-pipeline';  % <-- renamed per request
root     = getProjectRootSimple();
addpath(fullfile(root,'code',pipeline));

% Helper to safely make a directory (creates parents too)
ensure_dir = @(p) assert(mkdir(p) || exist(p,'dir')==7, ...
                         'Failed to create directory: %s', p);

% --- Derivatives layout (create everything if missing) ---
derivRoot  = fullfile(root, 'derivatives');  ensure_dir(derivRoot);
pipelineDir = fullfile(derivRoot, pipeline); ensure_dir(pipelineDir);
makeDerivativeDescription(pipelineDir);

datasetDir  = fullfile(pipelineDir, sub, 'eeg');  % <root>/derivatives/pipeline/sub-XX/eeg
ensure_dir(datasetDir);

% Subfolders used by this script
figsDir = fullfile(datasetDir, 'figs');
logDir  = fullfile(datasetDir, 'log');
ensure_dir(figsDir); ensure_dir(logDir);

% Filenames
fileBase      = sprintf('%s_task-%s', sub, task);
rawimportName = sprintf('%s_desc-rawimport_eeg.set', fileBase);
logFile       = fullfile(logDir, 'preprocessing_log.md');

% ------------------------------------------------------------------------
% Locate raw data (allow multiple formats)
rawDir = fullfile(root, sub, 'eeg');

% Preferred name patterns first, else auto-detect in this priority:
% 1) BrainVision .vhdr
% 2) EEGLAB .set
% 3) MATLAB .mat
% 4) EDF .edf  (if BIOSIG available)
cands = struct( ...
    'vhdr', fullfile(rawDir, sprintf('%s_eeg.vhdr', fileBase)), ...
    'set',  fullfile(rawDir, sprintf('%s_eeg.set',  fileBase)), ...
    'mat',  fullfile(rawDir, sprintf('%s_eeg.mat',  fileBase)), ...
    'edf',  fullfile(rawDir, sprintf('%s_eeg.edf',  fileBase)) );

% Auto-detect function
function p = pickFirstExisting(rawDir, specificPath, pattern)
    if exist(specificPath, 'file')
        p = specificPath; return;
    end
    d = dir(fullfile(rawDir, pattern));
    if isempty(d), p = ''; else, p = fullfile(rawDir, d(1).name); end
end

vhdrPath = pickFirstExisting(rawDir, cands.vhdr, '*.vhdr');
setPath  = pickFirstExisting(rawDir, cands.set,  '*.set');
matPath  = pickFirstExisting(rawDir, cands.mat,  '*.mat');
edfPath  = pickFirstExisting(rawDir, cands.edf,  '*.edf');

% Import + coords (format-agnostic)
locFile = fullfile(root, 'sourcedata', 'Standard-10-20-Cap81.ced');

if ~isempty(vhdrPath)
    fprintf('Importing BrainVision: %s\n', vhdrPath);
    [vhdrFolder, vhdrName, vhdrExt] = fileparts(vhdrPath);
    EEG = pop_loadbv(vhdrFolder, [vhdrName vhdrExt]);   % <- folder + filename
    EEG = eeg_checkset(EEG);

elseif ~isempty(setPath)
    fprintf('Importing EEGLAB .set: %s\n', setPath);
    [setFolder, setName, setExt] = fileparts(setPath);
    EEG = pop_loadset('filename', [setName setExt], 'filepath', setFolder); % <- split path
    EEG = eeg_checkset(EEG);

elseif ~isempty(matPath)
    fprintf('Importing from .mat: %s\n', matPath);
    EEG = import_from_mat(matPath);
    EEG = eeg_checkset(EEG);

elseif ~isempty(edfPath)
    fprintf('Importing EDF: %s\n', edfPath);
    try
        EEG = pop_biosig(edfPath);  % full path is fine here
        EEG = eeg_checkset(EEG);
    catch ME
        error(['EDF import requires BIOSIG in path. Original error: ' ME.message]);
    end
else
    error('No supported EEG files found in %s (looked for .vhdr/.set/.mat/.edf).', rawDir);
end


% Attach/lookup channel locations if needed
needsLocs = (~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs) || ...
             all(arrayfun(@(c) isempty(c.X) || isempty(c.Y) || isempty(c.Z), EEG.chanlocs)));
if needsLocs && exist(locFile,'file')==2
    EEG = pop_chanedit(EEG, 'lookup', locFile);
    EEG = eeg_checkset(EEG);
end

% Save WITH coords (or best effort)
pop_saveset(EEG, 'filename', rawimportName, 'filepath', datasetDir);


%% ========================================================================
%%                              audit
%% ========================================================================

% Open log and write header
fid = fopen(logFile, 'a');  % append; creates file if missing
fprintf(fid, '# Audit Report for %s_task-%s\n', sub, task);
fprintf(fid, 'Created: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% >>>>> 1. Channel info

fprintf(fid, '## Channel Information\n\n');

% Counts, flats, NaN/Inf
fprintf(fid, '- Channels: %d\n', EEG.nbchan);
isFlat  = std(double(EEG.data), 0, 2) == 0;
flatIdx = find(isFlat);
if ~isempty(flatIdx)
    fprintf(fid, '- Flat channels: %s\n', strjoin({EEG.chanlocs(flatIdx).labels}, ', '));
else
    fprintf(fid, '- Flat channels: none\n');
end
anyNaN = any(isnan(EEG.data(:)));
anyInf = any(isinf(EEG.data(:)));
fprintf(fid, '- NaNs present: %d\n', anyNaN);
fprintf(fid, '- Infs present: %d\n', anyInf);

% Rank vs channels
r = rank(double(EEG.data'));
fprintf(fid, '- Rank: %d vs nbchan %d\n', r, EEG.nbchan);

% Duplicate labels
labels = {EEG.chanlocs.labels};
hasDup = numel(labels) ~= numel(unique(labels));
fprintf(fid, '- Duplicate labels: %s\n', string(hasDup));

% Montage present (XYZ)
hasXYZ = arrayfun(@(c) ~isempty(c.X) && ~isempty(c.Y) && ~isempty(c.Z), EEG.chanlocs);
fprintf(fid, '- Channels with 3D coords: %d/%d\n', sum(hasXYZ), EEG.nbchan);

% 2D layout figure (save to figs/, link from log/ with ../figs/)
layoutName = sprintf('%s_layout2D.png', fileBase);
fig1       = fullfile(figsDir, layoutName);
figure('Visible','off');
topoplot([], EEG.chanlocs, 'style','blank', 'electrodes','labelpoint');
title('Channel layout (2D)'); axis off;
saveas(gcf, fig1); close(gcf);
fprintf(fid, '- Channel layout: ![layout](../figs/%s)\n', layoutName);

%% >>>>> 2. Recording parameters

fprintf(fid, '\n## Recording Parameters\n\n');
fprintf(fid, '- Trials: %d (1 = continuous)\n', EEG.trials);
fprintf(fid, '- Sampling rate: %g Hz\n', EEG.srate);
fprintf(fid, '- Samples: %d\n', EEG.pnts);

durSec = EEG.pnts / EEG.srate;
hh = floor(durSec/3600);
mm = floor(mod(durSec,3600)/60);
ss = round(mod(durSec,60));
fprintf(fid, '- Duration: %.1f s (≈ %02d:%02d:%02d)\n', durSec, hh, mm, ss);

%% >>>>> 3. Event info

fprintf(fid, '\n## Event Information\n\n');
nEv = numel(EEG.event);
fprintf(fid, '- Event count: %d\n', nEv);

if nEv > 0
    % Types & counts
    types = {EEG.event.type};
    [utypes, ~, idx] = unique(types);
    counts = accumarray(idx(:), 1);
    fprintf(fid, '- Types & counts:\n');
    for i = 1:numel(utypes)
        fprintf(fid, '  - %s: %d\n', string(utypes{i}), counts(i));
    end

    % Histogram over time
    evTimes   = [EEG.event.latency] / EEG.srate;
    evHistName = sprintf('%s_eventHist.png', fileBase);
    fig2       = fullfile(figsDir, evHistName);
    figure('Visible','off');
    histogram(evTimes, 50);
    xlabel('Time (s)'); ylabel('Count'); title('Event distribution');
    saveas(gcf, fig2); close(gcf);
    fprintf(fid, '- Event histogram: ![events](../figs/%s)\n', evHistName);

    % Sorted & in-bounds
    lat = [EEG.event.latency];
    sortedOK = issorted(lat) && lat(1) >= 1 && lat(end) <= EEG.pnts;
    fprintf(fid, '- Latencies sorted & in-bounds: %s\n', string(sortedOK));

    % Boundary markers
    nBoundary = sum(strcmpi({EEG.event.type}, 'boundary'));
    fprintf(fid, '- Boundary markers: %d\n', nBoundary);

    % Visibility field (BrainVision sometimes)
    if isfield(EEG.event, 'visible')
        visVals = [EEG.event.visible];
        if all(visVals==0)
            visStr = 'all hidden';
        elseif all(visVals==1)
            visStr = 'all visible';
        else
            visStr = 'mixed';
        end
        fprintf(fid, '- Event visibility: %s\n', visStr);
    end

    % First few event types (peek)
    peekN = min(5, nEv);
    fprintf(fid, '- First %d events: %s\n', peekN, strjoin(string(types(1:peekN)), ', '));
end

%% >>>>> 4. Sanity Checks

fprintf(fid, '\n## Sanity Checks\n\n');

% Average PSD (quick glance)
[spec, freqs] = spectopo(EEG.data, 0, EEG.srate, 'plot','off');
avgSpec = mean(spec,1);
psdName = sprintf('%s_avgPSD.png', fileBase);
fig3    = fullfile(figsDir, psdName);
figure('Visible','off');
plot(freqs, avgSpec, 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Log Power (dB)'); title('Average PSD'); grid on;
saveas(gcf, fig3); close(gcf);
fprintf(fid, '- Average PSD: ![psd](../figs/%s)\n', psdName);

% Amplitude sanity (temp hp1 on a copy)
EEGtmp = pop_eegfiltnew(EEG, 1, []);
pcts = prctile(abs(double(EEGtmp.data(:))), [95 99.9]);
fprintf(fid, '- Amplitude (hp1): 95th=%.1f µV, 99.9th=%.1f µV\n', pcts(1), pcts(2));

% Structure consistency
EEG = eeg_checkset(EEG);
fprintf(fid, '- eeg_checkset run (see console for warnings)\n');

%% >>>>> End
fprintf(fid, '\n---\nAudit complete. Saved at %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fclose(fid);


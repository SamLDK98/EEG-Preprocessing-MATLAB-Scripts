function [eegdir_erp, stem_erp, saveStepERP, logPath_erp] = create_erp_branch(EEG, varargin)
% CREATE_ERP_BRANCH  Prepare a parallel ERP derivatives branch.
%
% Usage:
%   [eegdir_erp, stem_erp, saveStepERP, logPath_erp] = create_erp_branch(EEG, 'root', root, 'stem', stem)
%
% Inputs (name/value, all optional but recommended):
%   'root' : BIDS project root (string). If omitted, inferred from EEG.filepath or pwd.
%   'stem' : filename stem like 'sub-01_task-CogAss1'. If omitted, inferred from EEG or current dir.
%   'derivName' : derivatives folder name (default 'eeg-ERP').
%
% Outputs:
%   eegdir_erp : path to derivatives/<derivName>/<sub>/eeg
%   stem_erp   : '<sub>_task-<task>' stem used in ERP branch
%   saveStepERP: @(desc,EEG) -> saves into ERP branch with description tag
%   logPath_erp: path to ERP preprocessing_log.md

% -------- Parse inputs
p = inputParser;
addParameter(p, 'root', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'stem', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'derivName', 'eeg-ERP', @(s)ischar(s)||isstring(s));
parse(p, varargin{:});
root      = string(p.Results.root);
stem_in   = string(p.Results.stem);
derivName = string(p.Results.derivName);

% -------- Infer a candidate path to locate BIDS root
% Try EEG.filepath first, then pwd
if ~isempty(EEG) && isfield(EEG,'filepath') && ~isempty(EEG.filepath)
    here = string(EEG.filepath);
else
    here = string(pwd);
end

% If root wasn't provided, infer it by trimming up to /derivatives/
if root == ""
    parts = strsplit(here, filesep);
    idxDeriv = find(strcmpi(parts,'derivatives'), 1, 'last');
    if ~isempty(idxDeriv)
        % BIDS root is everything before 'derivatives'
        root = fullfile(parts{1:idxDeriv-1});
    else
        % Fall back to going up two levels (project root guess)
        root = fileparts(fileparts(here));
    end
end
assert(isfolder(root), 'Could not determine a valid BIDS root.');

% -------- Determine stem, sub, task
% If stem is supplied, use it. Else try EEG.setname or a *.set in folder.
stem = stem_in;
if stem == ""
    if isfield(EEG,'setname') && ~isempty(EEG.setname)
        % Try to pull sub/task from setname; if that fails we'll fall back
        stem = erase(string(EEG.setname), "_desc-ica_eeg.set");
    else
        d = dir(fullfile(here, '*.set'));
        assert(~isempty(d), 'No .set file found to infer subject/task.');
        stem = erase(string(d(1).name), "_desc-ica_eeg.set");
        stem = erase(stem, ".set");
    end
end

% Extract tokens robustly
tokSub  = regexp(stem, '(sub-[A-Za-z0-9]+)', 'tokens', 'once');
tokTask = regexp(stem, '(?<=task-)[A-Za-z0-9]+', 'match', 'once');
assert(~isempty(tokSub),  'Could not infer subject label from stem.');
assert(~isempty(tokTask), 'Could not infer task label from stem.');
sub  = string(tokSub{1});     % e.g., 'sub-01'
task = string(tokTask);       % e.g., 'CogAss1'

% -------- Derivatives/ERP paths
derivRoot  = fullfile(root, 'derivatives', derivName);
subDir     = fullfile(derivRoot, sub);
eegdir_erp = fullfile(subDir, 'eeg');
logdir_erp = fullfile(eegdir_erp, 'log');
figsdir_erp= fullfile(eegdir_erp, 'figs');

cellfun(@(p) ~exist(p,'dir') && mkdir(p), {derivRoot, subDir, eegdir_erp, logdir_erp, figsdir_erp});

% Create dataset_description.json if missing
ddjson = fullfile(derivRoot, 'dataset_description.json');
if ~exist(ddjson, 'file')
    S = struct;
    S.Name        = char(derivName);
    S.BIDSVersion = '1.8.0';
    S.DatasetType = 'derivative';
    S.GeneratedBy = struct('Name','ERP preprocessing (EEGLAB+ERPLAB)', ...
                           'Version','1.0.0', ...
                           'Description','ERP branch: gentle filters, ICA transfer, ERPLAB epoching/averaging', ...
                           'CodeURL','');
    S.SourceDatasets = struct('Name','', 'URL','', 'DOI','');
    txt = jsonencode(S, 'PrettyPrint', true);
    fid = fopen(ddjson,'w'); assert(fid>0, 'Cannot write dataset_description.json');
    fwrite(fid, txt, 'char'); fclose(fid);
    fprintf('Created dataset_description.json at %s\n', ddjson);
else
    fprintf('Found dataset_description.json (kept): %s\n', ddjson);
end

% ERP log header (if new)
logPath_erp = fullfile(logdir_erp, 'preprocessing_log.md');
if ~exist(logPath_erp,'file')
    fid = fopen(logPath_erp,'w');
    ts  = datestr(now,'yyyy-mm-dd HH:MM:SS');
    fprintf(fid, "# ERP Preprocessing Log — %s (%s)\n\n", sub, ts);
    fprintf(fid, "- Derivatives branch: %s\n", derivName);
    fprintf(fid, "- Task: %s\n\n", task);
    fclose(fid);
end

% Build ERP stem and saver
stem_erp   = sprintf('%s_task-%s', sub, task);
saveStepERP = @(desc,EEGin) pop_saveset(EEGin, ...
    'filename', sprintf('%s_desc-%s_eeg.set', stem_erp, desc), ...
    'filepath', eegdir_erp);

fprintf('ERP derivatives ready at: %s\n', eegdir_erp);
end

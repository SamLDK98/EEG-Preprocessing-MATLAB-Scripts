function [EEG, newBounds] = manual_segment_prune(EEG, saveSet, logMsg, stem)
% Mark bad segments in the viewer, close it, and this will apply the cuts.
% Works across EEGLAB versions (no reliance on pop_eegplot return).
narginchk(1,4);
if nargin < 2 || isempty(saveSet), saveSet = []; end
if nargin < 3 || isempty(logMsg),  logMsg  = []; end
if nargin < 4, stem = ''; end

% BEFORE snapshot
pre.range   = [EEG.xmin EEG.xmax];
ev          = EEG.event;
pre.bound_s = (([ev(strcmpi({ev.type},'boundary')).latency]-1)/EEG.srate) + EEG.xmin;

% --- Segment viewer with an Apply button ---
assignin('base','EEG', EEG);                 % make EEG visible to eegplot's command
delete(findall(0,'tag','EEGPLOT'));          % no duplicate windows
eegplot(EEG.data, ...
    'srate', EEG.srate, ...
    'events', EEG.event, ...
    'winlength', 5, ...
    'butlabel', 'Apply cuts', ...
    'command', 'EEG = eeg_eegrej(EEG, TMPREJ); TMPREJ = [];');  % button applies

hFig = findall(0,'tag','EEGPLOT');
if ~isempty(hFig), waitfor(hFig(1)); end    % wait until closed
if evalin('base','exist(''EEG'',''var'')'), EEG = evalin('base','EEG'); end
evalin('base','clear TMPREJ');               % tidy up if left over

% If user marked regions, eegplot stored them in base variable TMPREJ.
newBounds = [];
if evalin('base','exist(''TMPREJ'',''var'')')
    tr = evalin('base','TMPREJ');                   % Nx2 samples (start,end)
    if ~isempty(tr)
        EEG = eeg_eegrej(EEG, tr);                  % apply cuts; inserts 'boundary' events
    end
    evalin('base','clear TMPREJ');                  % clean up
end

% AFTER snapshot
post.range   = [EEG.xmin EEG.xmax];
ev           = EEG.event;
post.bound_s = (([ev(strcmpi({ev.type},'boundary')).latency]-1)/EEG.srate) + EEG.xmin;

% New boundary times (seconds)
newBounds = setdiff(round(post.bound_s,3), round(pre.bound_s,3));

% Save + log (optional)
if ~isempty(saveSet), saveSet('cleanSegs', EEG); end
if ~isempty(logMsg)
    logMsg('### Manual segment pruning');
    logMsg('- Time range before: [%.3f  %.3f] s', pre.range(1), pre.range(2));
    logMsg('- Time range after:  [%.3f  %.3f] s',  post.range(1), post.range(2));
    if ~isempty(newBounds)
        logMsg('- New boundary markers at (s): %s', strjoin(compose('%.2f', sort(newBounds)), ', '));
    else
        logMsg('- New boundary markers: none');
    end
    logMsg('- Saved as: %s_desc-cleanSegs_eeg.set', stem);
end

% Metadata
EEG.etc.ica_prep.manual_segments = struct( ...
    'range_before_s', pre.range, 'range_after_s', post.range, ...
    'new_boundaries_s', sort(newBounds));
end

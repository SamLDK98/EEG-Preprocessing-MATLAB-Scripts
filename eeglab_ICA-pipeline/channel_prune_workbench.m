function [EEG_out, removed] = channel_prune_workbench(EEG_in, saveSet, logMsg, stem)
% Interactive channel removal sandbox; nothing is committed until 'apply'.
% Usage:
%   [EEG, removedChans] = channel_prune_workbench(EEG, saveSet, logMsg, stem)

if nargin < 2, saveSet = []; end
if nargin < 3, logMsg  = []; end
if nargin < 4, stem    = '';  end

base = EEG_in; EEG_out = EEG_in; removed = {};
labelset = {base.chanlocs.labels};
asStr = @(c) strjoin(c, ', ');

    function i_show()
        delete(findall(0,'tag','EEGPLOT'));
        Etmp = pop_select(base,'nochannel',removed);
        pop_eegplot(Etmp, 1, 1, 0);   % reject=0 -> NO "Reject" button here
        drawnow;
    end

    function L = i_valid(Lin)
        if isempty(Lin), L = {}; return; end
        L = regexp(Lin,'[^\s,;]+','match');
        bad = setdiff(L, labelset);
        if ~isempty(bad)
            fprintf('Warning: unknown labels ignored: {%s}\n', strjoin(bad, ', '));
            L = setdiff(L, bad, 'stable');
        end
    end

i_show();  % auto-open viewer on entry
fprintf('Commands: rm/add/set/del/clear/list/sugg/plot/psd/apply/cancel\n');

while true
    fprintf('\n[chprune] current removal set: {%s}\n', asStr(removed));
    cmd = strtrim(input('chprune> ','s')); if isempty(cmd), cmd = 'plot'; end
    toks = regexp(cmd,'^\s*(\w+)\s*(.*)$','tokens','once'); if isempty(toks), continue; end
    op = lower(toks{1}); args = strtrim(toks{2}); L = i_valid(args);

    switch op
        case {'rm','add'}, removed = unique([removed, L], 'stable'); i_show();
        case 'set',       removed = unique(L, 'stable');             i_show();
        case 'del',       removed = setdiff(removed, L, 'stable');   i_show();
        case 'clear',     removed = {};                              i_show();
        case 'list',      fprintf('Remove: {%s}\n', asStr(removed));
        case 'sugg'
            try
                [~, idx] = pop_rejchan(base,'threshold',5,'norm','on','measure','kurt','elec',1:base.nbchan,'plot','off');
                sug = {base.chanlocs(idx).labels};
                if isempty(sug), fprintf('No suggestions.\n'); else, fprintf('Suggested: %s\n', asStr(sug)); end
            catch, fprintf('No suggestions (rejchan unavailable).\n');
            end
        case 'plot', i_show();
        case 'psd'
            Etmp = pop_select(base,'nochannel',removed);
            [S0,f] = spectopo(double(base.data),0,base.srate,'plot','off');
            [S1,~] = spectopo(double(Etmp.data),0,Etmp.srate,'plot','off');
            figure; plot(f,mean(S0,1),'LineWidth',1.5); hold on; plot(f,mean(S1,1),'LineWidth',1);
            xlabel('Hz'); ylabel('Log Power (dB)'); title('Avg PSD: base vs candidate'); legend('base','candidate'); grid on;
        case 'apply'
    Etmp = pop_select(base,'nochannel',removed);
    % PRESERVE COORDINATES - pop_select strips them!
    remaining_labels = {Etmp.chanlocs.labels};
    original_labels = {base.chanlocs.labels};
    for i = 1:length(remaining_labels)
        orig_idx = strcmp(original_labels, remaining_labels{i});
        if any(orig_idx)
            orig_chan = base.chanlocs(orig_idx);
            Etmp.chanlocs(i).X = orig_chan.X;
            Etmp.chanlocs(i).Y = orig_chan.Y;
            Etmp.chanlocs(i).Z = orig_chan.Z;
            Etmp.chanlocs(i).theta = orig_chan.theta;
            Etmp.chanlocs(i).radius = orig_chan.radius;
            Etmp.chanlocs(i).sph_theta = orig_chan.sph_theta;
            Etmp.chanlocs(i).sph_phi = orig_chan.sph_phi;
            Etmp.chanlocs(i).sph_radius = orig_chan.sph_radius;
        end
    end
    
    EEG_out = eeg_checkset(Etmp);
    
    % ADD BACK THE MISSING SAVE AND LOG CODE:
    if ~isempty(saveSet), saveSet('rmchans', EEG_out); end
    if ~isempty(logMsg)
        logMsg('### Channel removal (manual)');
        if isempty(removed), logMsg('- Removed channels: none');
        else, logMsg('- Removed channels (%d): %s', numel(removed), asStr(removed)); end
        logMsg('- Saved as: %s_desc-rmchans_eeg.set', stem);
    end
    EEG_out.etc.ica_prep.removed_channels = removed;
    return  % <-- THIS WAS MISSING!
    
    EEG_out = eeg_checkset(Etmp);
        case 'cancel'
            fprintf('Canceled. No channels removed.\n'); return
        otherwise
            fprintf('Unknown. Try: rm/add/set/del/clear/list/sugg/plot/psd/apply/cancel\n');
    end
end
end

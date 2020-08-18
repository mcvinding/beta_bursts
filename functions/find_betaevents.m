function [output, rhomat] = find_betaevents(cfg, data)
% Output is a data structure of N length, containing A) number of beta
% events, B) start sample of beta events in data, C) end sample of beta
% event in data, C) length of events. rhomat is a vector of length N
% contining correaltion between number of events and amplitude/power in
% epochs. N is the number of steps used to determine threshold. When
% thereshold is determined then a single scalar is used.
%
% USE: [output, rhomat] = find_betaevents(cfg, data)
% INPUT:
% cfg.steps         = [Nx1] Steps to find correlation and threshold
%                     (required). Should be steps in order of medians/sds for
%                     cutofftype = 'sd'/'med' (e.g. 1:0.1:3) or percentile
%                     values (e.g. [70,75,80,85]) for 'pct'.
% cfg.cutofftype    = ['sd'/'med'/'pct'] Should cutoff be determined as
%                     x*standard deviations away from median (Little et al,
%                     2019). x*medians away from median (Shin et al. 2016).
%                     Based on percentile of the total signal variation
%                     (Tinkhauser et al. 2017) (default='med').
% cfg.halfmax       = ['yes'/'no'/'mixed'] redefine timing of events as half-
%                     maximum of the peak value. Overlapping peaks are counted 
%                     as one event. Otherwise, event length is defined by
%                     threshold crossing. 'Mixed' replace half-maximum
%                     values above threshold with the value of the
%                     threshold crossing (default='no').
% cfg.makeplot      = ['yes'/'no'] plot the time series with cutoff, peaks,
%                     and events marked (default='no').
% cfg.minlen        = [num or empty] specify the minimun required length of
%                     a burst (in seconds).
% cfg.mindist       = [num or empty] specify the minimun required length
%                     between to concecutive burst (in seconds).
%
%   Options for CUTOFFTYPE = 'sd / 'med'
% cfg.corrtype      = ['amp'/'pow'] What to correalte: amplitude of epochs 
%                     vs. number of events or power (amp^2) of epoch vs. 
%                     number of beta events (Little et al, 2019). (default='amp').
% cfg.length        = [num] length of epoch window in seconds. Passed to
%                     FT_REDEFINETRIAL (default=3).
% cfg.overlap       = [num] overlap between epochs. Passed to FT_REDEFINETRIAL
%                     (default=0, i.e. no overlap)
%
%   Options for CUTOFFTYPE = 'pct'
% cfg.pctentile     = [num] specified as cfg.steps.

%
% OUTPUT:
% ...

% TO DO:
% * Need a ft_checkdata section.

% opts
cfg = ft_checkconfig(cfg, 'required', 'steps');
steps = cfg.steps;

% cfg.getcutoff   = ft_getopt(cfg, 'getcutoff', 'no');
cfg.cutofftype  = ft_getopt(cfg, 'cutofftype',  'med');
cfg.corrtype    = ft_getopt(cfg, 'corrtype',    'amp');
cfg.makeplot    = ft_getopt(cfg, 'makeplot',    'no');
cfg.halfmax     = ft_getopt(cfg, 'halfmax',     'no');
cfg.minlen      = ft_getopt(cfg, 'minlen',      []);
cfg.mindist     = ft_getopt(cfg, 'mindist',     []);

% Check variables
if ~strcmp('pct', cfg.cutofftype)
    if ~any(strcmp({'sd','med'},cfg.cutofftype))
        error('Unknown cutoff method \"%s\".',cfg.cutofftype)
    end

    if ~any(strcmp({'amp','pow'}, cfg.corrtype))
        error('Unknown cutoff method \"%s\".',cfg.cutofftype)
    end
else
%     cfg = ft_checkconfig(cfg, 'required', 'pctentile');
end

if ~isempty(cfg.minlen)
    minlenSam = cfg.minlen*data.fsample;
end

if ~isempty(cfg.mindist)
    mindistSam = cfg.mindist*data.fsample;
end

% Remove sampleinfo for indexing later
if isfield(data, 'sampleinfo')
    data.sampleinfo = data.sampleinfo-data.sampleinfo(1)+1;
end

% Select channel data
if isfield(cfg, 'channel')
    cf = [];
    cf.channel = cfg.channel;
    data = ft_selectdata(cfg, data);
elseif size(data.trial{:},1) > 1
    error('Only works on single time series. Data has %i channels. Consider specifying cfg.channel', size(data.trial{:},1))
end

if any(strcmp({'sd','med'}, cfg.cutofftype))
    % Make pseudo-tirals.
    cfg.length      = ft_getopt(cfg, 'length', 3);
    cfg.overlap     = ft_getopt(cfg, 'overlap', 0);
    epo = ft_redefinetrial(cfg, data);

    trl = epo.sampleinfo;
    % trl = trl-trl(1,1)+1; %Corrent for non-zero sample offset

    % Get epoch power and amplitude
    epoamp = nan(length(epo.trial),1);
    epopow = nan(length(epo.trial),1);
    for n = 1:length(trl)
        epoamp(n) = mean(epo.trial{n});
        epopow(n) = mean(epo.trial{n}.^2);
    end
else
    trl = [];
end

% Cutoffs
% Initiate values
if ~isempty(trl); pkmat = zeros(length(trl),1); end
n_events = zeros(1,length(steps));
rhomat   = nan(length(steps),1);
cutoff   = zeros(1,length(steps));
bdat     = struct();

% Find summary values
dat = data.trial{:};
med = median(dat);
sd = std(dat);
fprintf('Median of time-series: %.3f. sd: %.3f.\n', med, sd)

for ii = 1:length(steps)
    if strcmp(cfg.cutofftype, 'sd')
        cutoff(ii) = med+sd*steps(ii);
    elseif strcmp(cfg.cutofftype, 'med')
        cutoff(ii) = med+med*steps(ii);
    elseif strcmp(cfg.cutofftype, 'pct')
        cutoff(ii) = prctile(dat, steps(ii));
    end
    burst = dat >= cutoff(ii);
    
    % Get initial summaries
    if burst(end) == 1; burst(end) = 0; end
    dburst = diff([0 burst]);
    n_events(ii) = sum(dburst==1);
    startb  = find(dburst==1);      % Start of burst
    endb    = find(dburst==-1);     % End of burst´
    
    maxarray = zeros(1,n_events(ii));
    maxidx   = zeros(1,n_events(ii));
    for n = 1:n_events(ii)
        [maxarray(n), maxidx(n)] = max(dat(startb(n):endb(n)));
    end

    maxidx = maxidx+startb-1;  % correct idx
        
    % start-stop based on half-max width (halfmax='yes' or 'mixed')
    if any(strcmp({'yes','mixed'}, cfg.halfmax))
        begsam = zeros(1,n_events(ii));
        endsam = zeros(1,n_events(ii));
        hlfmx = maxarray/2;       
        
        if strcmp(cfg.halfmax, 'mixed')
            hlfmx(hlfmx>cutoff(ii)) = cutoff(ii);
        end    
        
        for n = 1:n_events(ii)
            %start
            idx = maxidx(n);
            xval = dat(idx);
            while xval > hlfmx(n)
                if idx == 1
                    break
                end
                idx = idx-1;
                xval = dat(idx);
            end
            begsam(n) = idx;
            % end
            idx = maxidx(n);
            xval = dat(idx);
            while xval > hlfmx(n)
                if idx == length(dat)
                    break
                end
                idx = idx+1;
                xval = dat(idx);
            end
            endsam(n) = idx;
        end
        
        cburst = zeros(length(dat),1);
        for n = 1:n_events(ii)
            cburst(begsam(n):endsam(n)) = 1;
        end        
        burst = cburst;
    end

    % correct min time between events
    if ~isempty(cfg.mindist)
        dburst  = diff([0 burst]);
        startb  = find(dburst==1);              % Start of burst
        endb    = find(dburst==-1);             % End of burst
        intvl = startb(2:end)-endb(1:end-1);    % Length of interval
        short_intvl = intvl <= mindistSam;

        cburst = burst;
        for n = find(short_intvl)
            cburst(endb(n):startb(n+1)) = 1;
        end
        burst = cburst;
    end
        
    % correct min length
    if ~isempty(cfg.minlen)
        dburst  = diff([0 burst]);
        startb  = find(dburst==1);              % Start of burst
        endb    = find(dburst==-1);             % End of burst
        blen    = endb-startb;                  % Length of burst
        short_len = blen <= minlenSam;

        cburst = burst;
        for n = find(short_len)
            cburst(startb(n):endb(n)) = 0;
        end
        burst = cburst;
    end
        
    % Get final summaries
    if burst(end) == 1; burst(end) = 0; end
    dcburst = diff([0 burst]);
    neve = sum(dcburst==1);
    startb  = find(dcburst==1);      % Start of burst
    endb    = find(dcburst==-1);     % End of burst´

    maxarray = zeros(neve,1);
    maxidx   = zeros(neve,1);
    for n = 1:neve
        [maxarray(n), maxidx(n)] = max(dat(startb(n):endb(n)));
    end
    maxidx = maxidx+startb'-1;

    n_events(ii) = neve;
    
    % Length of burst
    blen    = endb-startb;          % Length of burst
    evelen   = blen/data.fsample;      % Length of burst is seconds  
    if any(blen<0)
        error('negative length of event')
    end

    % Arrange data
    bdat(ii).event  = [startb; endb]';
    bdat(ii).evelen = evelen';
    bdat(ii).maxpk  = maxarray;
    bdat(ii).maxidx = maxidx;
    
    
    % Get correlations    
    if any(strcmp({'sd','med'},cfg.cutofftype))
        for n = 1:length(trl)
            tmp = burst(trl(n,1):trl(n,2));
            pkmat(n) = sum((diff(tmp)==1));
        end
        
        if strcmp('amp', cfg.corrtype)
            rhomat(ii) = corr(pkmat, epoamp);
        elseif strcmp('pow', cfg.corrtype)
            rhomat(ii) = corr(pkmat, epopow);
        end
    end
end

% Make output
output.steps    = steps;
output.n_events = n_events;
output.cutoff   = cutoff;
output.bdat     = bdat;
output.chan     = data.label{:};
output.cfg      = cfg;

if strcmp(cfg.makeplot,'yes')
    addpath('./plotting')
    eventplot(output,dat)
end

% End
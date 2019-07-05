% Test of timelocked analysis
addpath /home/mikkel/fieldtrip/fieldtrip
ft_defaults
addpath /home/mikkel/beta_bursts/functions

%% Load envelope data
load('/home/mikkel/beta_bursts/sandbox/tmlck.mat')

%% Find events using find_eventsER
cfg = [];
cfg.threshold   = 1:0.2:2;
cfg.cutofftype  = 'med';
cfg.corrtype    = 'amp';

[~, rhomat] = find_betaeventsER(cfg, tmlck);

%% Find threshold

cutoff = find_threshold(rhomat', cfg.threshold, 1);

%% Find event at cutoff
cfg = [];
cfg.threshold   = cutoff;
cfg.cutofftype  = 'med';
cfg.corrtype    = 'amp';

[eventdata] = find_betaeventsER(cfg, tmlck);

%% Next: Raster plots...
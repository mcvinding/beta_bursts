function eventplot(b_summary, dat)
% Plot time series with marked events.
% USE: eventplot(b_summary, dat)

steps       = b_summary.steps;
cutoff      = b_summary.cutoff;
n_events    = b_summary.n_events;

for ii = 1:length(steps)
    
    bdat    = b_summary.bdat(ii);
    startb  = bdat.event(:,1);
    endb    = bdat.event(:,2);
    
    evemask = (zeros(size(dat)))==1;
    for n = 1:n_events(ii)
        evemask(startb(n):endb(n)) = 1;
    end
    
    evedat = nan(size(dat));
    evedat(evemask) = dat(evemask);

    if length(b_summary.steps) >= 8
        dim = ceil(length(steps)/8);
        subplot(dim,8,ii); hold on
    else
        figure; hold on
    end

    plot(dat);
    plot(repmat(cutoff(ii),length(dat),1),'r--');
    
    plot(bdat.maxidx, bdat.maxpk, 'ko');
    plot(1:length(dat), evedat, 'linewidth',2);
    xlim([0 length(dat)]);
    title(steps(ii));
    txt = ['N = ', num2str(n_events(ii))];
    text(length(dat)/20, max(dat), txt)

end

function ds = build_dataset_canova_forero(dconf, paths)
% Build Z matrix with transformations described in the paper. :contentReference[oaicite:23]{index=23}

ids = struct2cell(dconf.series);
idnames = fieldnames(dconf.series);

% Download + read each series
tts = cell(numel(ids),1);
for i=1:numel(ids)
    sid = ids{i};
    f = fullfile(paths.data_raw, sprintf('%s.csv', sid));
    if ~exist(f,'file')
        fprintf('[data] downloading %s...\n', sid);
        download_fred_series(sid, paths.data_raw);
    end
    tt = read_fred_csv(f);

    % aggregate to monthly (SP500 daily etc.)
    tt = to_monthly_last(tt);
    tts{i} = tt;
end

% Merge to one timetable
tt_all = merge_series(tts, ids);

% Rates divide by 100
for i=1:numel(dconf.transform.rate_div100)
    sid = dconf.transform.rate_div100{i};
    if ismember(sid, tt_all.Properties.VariableNames)
        tt_all.(sid) = tt_all.(sid) / 100.0;
    end
end

% Build variables
INDPRO = tt_all.INDPRO;
PCEPI  = tt_all.PCEPI;
UNRATE = tt_all.UNRATE;
FEDFUNDS = tt_all.FEDFUNDS;
GS10 = tt_all.GS10;
TB3MS = tt_all.TB3MS;
M2SL = tt_all.M2SL;
PPIACO = tt_all.PPIACO;
SP500 = tt_all.SP500;

Y  = transform_yoy_log(INDPRO);
P  = transform_yoy_log(PCEPI);
U  = UNRATE;
R  = FEDFUNDS;
Slope = GS10 - TB3MS;
M  = transform_yoy_log(M2SL);
Pcom = transform_yoy_log(PPIACO);
SPg  = transform_yoy_log(SP500);

Z = [Y P U R Slope M Pcom SPg];

ds = struct();
ds.dates = tt_all.Time;
ds.Z = Z;
ds.names = dconf.names;
ds.raw = tt_all;

end

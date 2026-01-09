function dconf = data_config()
% Data definitions as in the paper (Zt list + sample window).
%
% We download from FRED (no API key).
% Transformations: YoY log growth for levels; rates divided by 100.

dconf.sample_start = datetime(1963,1,1);
dconf.sample_end   = datetime(2023,6,1);  % paper sample 

% Variable order must match the paper / Table 1 (used in A.8)
% Zt = (Y, P, U, R, Slope, M, Pcom, SP500)'
dconf.names = {'IP_g_yoy','PCE_pi_yoy','UNRATE','FEDFUNDS','SLOPE_10Y_3M','M2_g_yoy','PPIACO_g_yoy','SP500_g_yoy'};

% FRED series IDs (choices for Pcom etc are explicit here)
dconf.series.IP      = 'INDPRO';     % Industrial Production Index
dconf.series.PCEPI   = 'PCEPI';      % PCE price index
dconf.series.UNRATE  = 'UNRATE';     % Unemployment rate
dconf.series.FEDFUNDS= 'FEDFUNDS';   % Fed funds
dconf.series.GS10    = 'GS10';       % 10y
dconf.series.TB3MS   = 'TB3MS';      % 3m
dconf.series.M2SL    = 'M2SL';       % M2
dconf.series.PPIACO  = 'PPIACO';     % PPI all commodities (proxy for commodity index)
dconf.series.SP500   = 'SPASTT01USM661N';  % 'SP500';  Pb with SP500 in FRED, starts only in 2016    % S&P 500 (daily -> monthly last)

% Transform flags
dconf.transform.yoy_log = {'INDPRO','PCEPI','M2SL','PPIACO','SP500'};            
dconf.transform.rate_div100 = {'UNRATE','FEDFUNDS','GS10','TB3MS'};

% How to aggregate higher frequency series (SP500 is daily):
dconf.agg_method = 'last'; % last obs in month

end

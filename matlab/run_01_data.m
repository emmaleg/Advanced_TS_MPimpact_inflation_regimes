%% run_01_data.m
% Build dataset Zt = (Y, P, U, R, Slope, M2, Pcom, SP500)' from FRED.
% Matches variable definitions in the paper.

clear; clc;

rootDir = fileparts(pwd); % if pwd = repo/matlab
% you have to position yourself in the root folder of the repo ! 
% TODO implement something more robust

% --- Load configs (package calls) ---
paths = cf.config.paths_config(rootDir);
dconf = cf.config.data_config();
mconf = cf.config.model_config();
pconf = cf.config.priors_config();

fprintf('[run_01_data] Download/build data...\n');

% --- Build dataset ---
ds = cf.data.build_dataset_canova_forero(dconf, paths);
ds_full = ds; % for plots 1960-2023, without restrictions

% Keep sample window
mask = (ds.dates >= dconf.sample_start) & (ds.dates <= dconf.sample_end);
ds.dates = ds.dates(mask);
ds.Z     = ds.Z(mask,:);
ds.names = ds.names;

% Drop initial NaNs due to YoY transforms (12 lags) etc.
good = all(isfinite(ds.Z),2);
ds.dates = ds.dates(good);
ds.Z     = ds.Z(good,:);

% Save (to data/processed by your paths_config)
outFile = fullfile(paths.data_processed, 'dataset_canova_forero.mat');
if ~exist(paths.data_processed,'dir'); mkdir(paths.data_processed); end
save(outFile,'ds','dconf','mconf','pconf','paths','-v7.3');

fprintf('[run_01_data] Saved: %s\n', outFile);

% --- Charts of the data ---

% Where to save figures (fits your repo architecture: output/...)
figDir = fullfile(rootDir, 'output', 'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

% Ensure dates are datetime for plotting
d = ds_full.dates;
if isnumeric(d)
    d = datetime(d, 'ConvertFrom','datenum');
elseif iscell(d)
    d = datetime(d);
end

% Names as cellstr
names = ds_full.names;
if isstring(names); names = cellstr(names); end

% -----------------------------
% 1) Overview figure (C.15-like)
% -----------------------------
n = size(ds_full.Z,2);

fig1 = figure('Name','Dataset overview','Color','w','Position',[100 100 1200 900]);
tlo = tiledlayout(fig1, 4, 2, 'TileSpacing','compact', 'Padding','compact');

for i = 1:n
    nexttile(tlo);
    plot(d, ds_full.Z(:,i), 'LineWidth', 1.1);
    grid on; box on;
    title(names{i}, 'Interpreter','none');
    xlim([datetime(1960,1,1) datetime(2023,12,31)]);
end

title(tlo, 'Canova-Forero dataset (overview)', 'FontWeight','bold');

% Save
exportgraphics(fig1, fullfile(figDir, 'dataset_overview_C15_style.png'), 'Resolution', 200);
exportgraphics(fig1, fullfile(figDir, 'dataset_overview_C15_style.pdf'));

% -----------------------------------------------
% 2) Inflation histogram + fitted normal distribution
% -----------------------------------------------

% ==== pick inflation from ds_full (NOT masked sample) ====
names_full = ds_full.names;
if isstring(names_full); names_full = cellstr(names_full); end

infl_idx = [];
cand = ["P","INF","INFL","INFLATION","PCE"];
for c = cand
    hit = find(strcmpi(names_full, c), 1);
    if ~isempty(hit); infl_idx = hit; break; end
end
if isempty(infl_idx)
    for c = cand
        hit = find(contains(upper(string(names_full)), upper(c)), 1);
        if ~isempty(hit); infl_idx = hit; break; end
    end
end
if isempty(infl_idx)
    infl_idx = 2;
    warning('Inflation column not confidently detected from ds_full.names; using column 2 by default.');
end

infl = ds_full.Z(:, infl_idx);
infl = infl(isfinite(infl));

%% --- Figure 2 style: histogram counts + fitted normal (scaled), with padding

xMin  = -0.04;
xMax  =  0.12;
xTick =  0.02;

binW  =  0.01;                 % IMPORTANT: bins fins (sinon y monte très haut)
edges = xMin:binW:xMax;

xPad  = 0.005;                 % marge visuelle avant/après
yPad  = 0.05;                  % 5% de marge au-dessus

% infl = ta série d'inflation (vector) déjà nettoyée (isfinite), en décimales
% (si tu veux l'histogramme sur TON sample d'estimation, utilise ds.Z après mask/good.
%  si tu veux 1960-2023, utilise ds_full.Z avant mask.)

fig2 = figure('Name','PCE Inflation','Color','w','Position',[200 150 900 600]);
ax = axes(fig2); hold(ax,'on'); box(ax,'on');

% Histogramme en counts avec bins imposés
h = histogram(ax, infl, 'BinEdges', edges, 'Normalization','count', ...
    'FaceColor',[0 0 1], 'EdgeColor','k');

% Courbe rouge = normale ajustée, rescalée en counts
N   = numel(infl);
mu  = mean(infl);
sig = std(infl);

x = linspace(xMin, xMax, 400);
y = normpdf(x, mu, sig) * N * binW;   % pdf -> counts (même unité que l'hist)

plot(ax, x, y, 'r', 'LineWidth', 2);

title(ax, 'PCE Inflation', 'FontWeight','bold');

% Axes: ticks + marges
xlim(ax, [xMin-xPad, xMax+xPad]);
xticks(ax, xMin:xTick:xMax);

% Limite Y: un peu au-dessus du max des barres
yMaxBar = max(h.Values);
ylim(ax, [0, (1+yPad)*yMaxBar]);

% Ticks Y "comme le papier"
yl = ylim(ax);
yticks(ax, 0:20:ceil(yl(2)/20)*20);

% Optionnel : look plus "papier"
set(ax, 'TickDir','out', 'Layer','top');
grid(ax,'off');


% Save
exportgraphics(fig2, fullfile(figDir, 'inflation_hist_normalfit.png'), 'Resolution', 200);
exportgraphics(fig2, fullfile(figDir, 'inflation_hist_normalfit.pdf'));

function fig = inflation_hist_normalfit(ds, figDir, varargin)
%CF.PLOT.INFLATION_HIST_NORMALFIT Histogram of inflation + fitted normal (scaled to counts).
%   fig = cf.plot.inflation_hist_normalfit(ds, figDir, ...)
%
% Required:
%   ds    : struct with fields .Z (T×n), .names
%   figDir: output directory for figures
%
% Name-value options:
%   'Candidates'  : string array of candidate names (default: ["P","INF","INFL","INFLATION","PCE"])
%   'DefaultCol'  : fallback column index (default: 2)
%   'XMin'        : default -0.04
%   'XMax'        : default  0.12
%   'XTick'       : default  0.02
%   'BinW'        : default  0.01
%   'XPad'        : default  0.005
%   'YPad'        : default  0.05
%   'Title'       : default 'PCE Inflation'
%   'FileStem'    : default 'inflation_hist_normalfit'

ip = inputParser;
ip.addRequired('ds', @(s)isstruct(s) && isfield(s,'Z') && isfield(s,'names'));
ip.addRequired('figDir', @(x)ischar(x) || isstring(x));

ip.addParameter('Candidates', ["P","INF","INFL","INFLATION","PCE"], @(x)isstring(x) || iscellstr(x));
ip.addParameter('DefaultCol', 2, @(x)isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('XMin', -0.04, @(x)isnumeric(x) && isscalar(x));
ip.addParameter('XMax',  0.12, @(x)isnumeric(x) && isscalar(x));
ip.addParameter('XTick', 0.02, @(x)isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('BinW',  0.01, @(x)isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('XPad',  0.005, @(x)isnumeric(x) && isscalar(x) && x>=0);
ip.addParameter('YPad',  0.05, @(x)isnumeric(x) && isscalar(x) && x>=0);
ip.addParameter('Title', 'PCE Inflation', @(x)ischar(x) || isstring(x));
ip.addParameter('FileStem', 'inflation_hist_normalfit', @(x)ischar(x) || isstring(x));

ip.parse(ds, figDir, varargin{:});
opt = ip.Results;

if ~exist(opt.figDir,'dir'); mkdir(opt.figDir); end

names = ds.names;
if isstring(names); names = cellstr(names); end

Z = ds.Z;

% --- detect inflation column ---
infl_idx = [];
cand = string(opt.Candidates);

for c = cand
    hit = find(strcmpi(names, c), 1);
    if ~isempty(hit); infl_idx = hit; break; end
end
if isempty(infl_idx)
    for c = cand
        hit = find(contains(upper(string(names)), upper(c)), 1);
        if ~isempty(hit); infl_idx = hit; break; end
    end
end
if isempty(infl_idx)
    infl_idx = opt.DefaultCol;
    warning('Inflation column not confidently detected from ds.names; using column %d by default.', infl_idx);
end

infl = Z(:, infl_idx);
infl = infl(isfinite(infl));

% --- histogram params ---
xMin  = opt.XMin;
xMax  = opt.XMax;
xTick = opt.XTick;
binW  = opt.BinW;
edges = xMin:binW:xMax;

xPad  = opt.XPad;
yPad  = opt.YPad;

fig = figure('Name','Inflation histogram','Color','w','Position',[200 150 900 600]);
ax = axes(fig); hold(ax,'on'); box(ax,'on');

h = histogram(ax, infl, 'BinEdges', edges, 'Normalization','count', ...
    'FaceColor',[0 0 1], 'EdgeColor','k');

% fitted normal scaled to counts
N   = numel(infl);
mu  = mean(infl);
sig = std(infl);

x = linspace(xMin, xMax, 400);
y = normpdf(x, mu, sig) * N * binW;

plot(ax, x, y, 'r', 'LineWidth', 2);

title(ax, opt.Title, 'FontWeight','bold');

xlim(ax, [xMin-xPad, xMax+xPad]);
xticks(ax, xMin:xTick:xMax);

yMaxBar = max(h.Values);
ylim(ax, [0, (1+yPad)*yMaxBar]);

yl = ylim(ax);
yticks(ax, 0:20:ceil(yl(2)/20)*20);

set(ax, 'TickDir','out', 'Layer','top');
grid(ax,'off');

% Save
stem = char(opt.FileStem);
exportgraphics(fig, fullfile(opt.figDir, stem + ".png"), 'Resolution', 200);
exportgraphics(fig, fullfile(opt.figDir, stem + ".pdf"));
end

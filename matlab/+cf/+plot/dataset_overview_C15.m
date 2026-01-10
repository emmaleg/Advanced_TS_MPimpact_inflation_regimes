function fig = dataset_overview_C15(ds, figDir, varargin)
%CF.PLOT.DATASET_OVERVIEW_C15 Overview tiled plot (C.15-like) for a dataset struct.
%   fig = cf.plot.dataset_overview_C15(ds, figDir, ...)
%
% Required:
%   ds    : struct with fields .dates, .Z (T×n), .names
%   figDir: output directory for figures
%
% Name-value options:
%   'XLim'        : [datetime datetime] (default: [1960-01-01, 2023-12-31])
%   'TileRows'    : integer (default: 4)
%   'TileCols'    : integer (default: 2)
%   'LineWidth'   : numeric (default: 1.1)
%   'FileStem'    : char/string (default: 'dataset_overview_C15_style')
%   'Title'       : char/string (default: 'Canova-Forero dataset (overview)')

ip = inputParser;
ip.addRequired('ds', @(s)isstruct(s) && isfield(s,'Z') && isfield(s,'dates') && isfield(s,'names'));
ip.addRequired('figDir', @(x)ischar(x) || isstring(x));
ip.addParameter('XLim', [datetime(1960,1,1) datetime(2023,12,31)], @(x)isdatetime(x) && numel(x)==2);
ip.addParameter('TileRows', 4, @(x)isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('TileCols', 2, @(x)isnumeric(x) && isscalar(x) && x>=1);
ip.addParameter('LineWidth', 1.1, @(x)isnumeric(x) && isscalar(x) && x>0);
ip.addParameter('FileStem', 'dataset_overview_C15_style', @(x)ischar(x) || isstring(x));
ip.addParameter('Title', 'Canova-Forero dataset (overview)', @(x)ischar(x) || isstring(x));
ip.parse(ds, figDir, varargin{:});
opt = ip.Results;

if ~exist(opt.figDir,'dir'); mkdir(opt.figDir); end

% Dates -> datetime
d = ds.dates;
if isnumeric(d)
    d = datetime(d, 'ConvertFrom','datenum');
elseif iscell(d) || isstring(d) || ischar(d)
    d = datetime(d);
end

% Names -> cellstr
names = ds.names;
if isstring(names); names = cellstr(names); end

Z = ds.Z;
n = size(Z,2);

fig = figure('Name','Dataset overview','Color','w','Position',[100 100 1200 900]);
tlo = tiledlayout(fig, opt.TileRows, opt.TileCols, 'TileSpacing','compact', 'Padding','compact');

for i = 1:n
    nexttile(tlo);
    plot(d, Z(:,i), 'LineWidth', opt.LineWidth);
    grid on; box on;
    title(names{i}, 'Interpreter','none');
    xlim(opt.XLim);
end

title(tlo, opt.Title, 'FontWeight','bold');

% Save
stem = char(opt.FileStem);
exportgraphics(fig, fullfile(opt.figDir, stem + ".png"), 'Resolution', 200);
exportgraphics(fig, fullfile(opt.figDir, stem + ".pdf"));
end

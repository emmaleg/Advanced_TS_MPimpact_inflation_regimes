function fig = make_fig_irf_regimes(irfRes, ds, paths, fileStem, titleStr)
%CF.PLOT.MAKE_FIG_IRF_REGIMES
% Plot IRFs for Low vs High regimes with "paper style":
%  - Low regime: BLUE, SOLID (median + bands)
%  - High regime: RED, DASHED (median + bands)
%  - Legend on the RIGHT side (outside the grid), no "Ignoring extra legend entries" warnings.
%
% Inputs:
%   irfRes.h              (H x 1) horizons, should start at 0
%   irfRes.low.median     (H x n)
%   irfRes.low.p16        (H x n)
%   irfRes.low.p84        (H x n)
%   irfRes.high.median    (H x n)
%   irfRes.high.p16       (H x n)
%   irfRes.high.p84       (H x n)
%   ds.names              (1 x n) cellstr
%   paths.output_figures  folder
%
% Saves PDF+PNG to output/figures.

h     = irfRes.h(:);
names = ds.names;

% --- Styling (match paper) ---
blue = [0.0000 0.4470 0.7410];
red  = [0.8500 0.3250 0.0980];

lw_med  = 2.4;
lw_band = 1.3;

ls_low  = '-';   % solid
ls_high = '--';  % dashed

% Robust axis range
hmax = max(h);
if isempty(hmax) || ~isfinite(hmax); hmax = 36; end

fig = figure('Color','w','Position',[100 100 1300 720]);

% Layout: 3x4 -> 8 plots in first 3 columns; legend occupies the 4th column
tl = tiledlayout(fig, 3, 4, 'TileSpacing','compact', 'Padding','compact');

plotTiles = [1 2 3, 5 6 7, 9 10]; % 8 panels
panelOrder = 1:8;

for k = 1:8
    j = panelOrder(k);
    ax = nexttile(tl, plotTiles(k)); %#ok<NASGU>
    hold on;

    % LOW (blue solid): median + bands (same style as requested)
    plot(h, irfRes.low.median(:,j), 'LineStyle',ls_low,  'Color',blue, 'LineWidth',lw_med);
    plot(h, irfRes.low.p16(:,j),    'LineStyle',ls_low,  'Color',blue, 'LineWidth',lw_band);
    plot(h, irfRes.low.p84(:,j),    'LineStyle',ls_low,  'Color',blue, 'LineWidth',lw_band);

    % HIGH (red dashed): median + bands
    plot(h, irfRes.high.median(:,j),'LineStyle',ls_high, 'Color',red,  'LineWidth',lw_med);
    plot(h, irfRes.high.p16(:,j),   'LineStyle',ls_high, 'Color',red,  'LineWidth',lw_band);
    plot(h, irfRes.high.p84(:,j),   'LineStyle',ls_high, 'Color',red,  'LineWidth',lw_band);

    yline(0,'-','Color',[0 0 0],'LineWidth',0.8);

    xlim([0 hmax]);
    xticks(0:6:hmax);

    title(string(names{j}), 'Interpreter','none');
    grid on;
    box on;
end

% --- Legend axis (right column spanning 3 rows) ---
lgax = nexttile(tl, 4, [3 1]);
axis(lgax,'off');
hold(lgax,'on');

% Dummy handles for legend (no warnings, stable styling)
h1 = plot(lgax, NaN, NaN, 'LineStyle',ls_low,  'Color',blue, 'LineWidth',lw_med);
h2 = plot(lgax, NaN, NaN, 'LineStyle',ls_low,  'Color',blue, 'LineWidth',lw_band);
h3 = plot(lgax, NaN, NaN, 'LineStyle',ls_high, 'Color',red,  'LineWidth',lw_med);
h4 = plot(lgax, NaN, NaN, 'LineStyle',ls_high, 'Color',red,  'LineWidth',lw_band);

lgd = legend(lgax, [h1 h2 h3 h4], ...
    {'Low: median','Low: bands (p16/p84)','High: median','High: bands (p16/p84)'}, ...
    'Location','northwest', 'Box','off');
lgd.FontSize = 10;

% Title
sgtitle(tl, titleStr, 'FontWeight','bold');

% --- Save ---
if ~exist(paths.output_figures,'dir'); mkdir(paths.output_figures); end
pdfFile = fullfile(paths.output_figures, fileStem + ".pdf");
pngFile = fullfile(paths.output_figures, fileStem + ".png");

exportgraphics(fig, pdfFile, 'ContentType','vector');
exportgraphics(fig, pngFile, 'Resolution',220);

end

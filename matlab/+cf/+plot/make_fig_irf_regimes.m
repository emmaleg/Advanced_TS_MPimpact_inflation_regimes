function fig = make_fig_irf_regimes(irfRes, ds, paths, fileStem, titleStr)
% irfRes: struct with fields low/high each having median,p16,p84 (H x n)
h = irfRes.h;
names = ds.names;

fig = figure('Color','w','Position',[100 100 1100 700]);

panelOrder = 1:8;
for k=1:8
    ax = subplot(3,3,k); %#ok<LAXES>
    j = panelOrder(k);

    plot(ax, h, irfRes.low.median(:,j),  'LineWidth',1.8); hold(ax,'on');
    plot(ax, h, irfRes.low.p16(:,j),     '--', 'LineWidth',1.0);
    plot(ax, h, irfRes.low.p84(:,j),     '--', 'LineWidth',1.0);

    plot(ax, h, irfRes.high.median(:,j), 'LineWidth',1.8);
    plot(ax, h, irfRes.high.p16(:,j),    '--', 'LineWidth',1.0);
    plot(ax, h, irfRes.high.p84(:,j),    '--', 'LineWidth',1.0);

    yline(ax,0,'-');
    xlim(ax,[1 max(h)]);
    xticks(ax,6:6:max(h));
    title(ax, string(names{j}), 'Interpreter','none');
    grid(ax,'on');
end

% legend in the 9th slot
ax9 = subplot(3,3,9);
axis(ax9,'off');
legend(ax9, {'Low median','Low p16','Low p84','High median','High p16','High p84'}, ...
    'Location','south', 'Box','off');

sgtitle(titleStr);

if ~exist(paths.output_figures,'dir'); mkdir(paths.output_figures); end
pdfFile = fullfile(paths.output_figures, fileStem + ".pdf");
pngFile = fullfile(paths.output_figures, fileStem + ".png");
exportgraphics(fig, pdfFile, 'ContentType','vector');
exportgraphics(fig, pngFile, 'Resolution',200);

end

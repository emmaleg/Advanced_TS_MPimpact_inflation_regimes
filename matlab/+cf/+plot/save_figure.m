function save_figure(fig, filePath)
% Robust save helper (exportgraphics if available, otherwise print)

filePath = char(filePath);
outDir = fileparts(filePath);
if ~isempty(outDir) && ~exist(outDir,'dir')
    mkdir(outDir);
end

[~,~,ext] = fileparts(filePath);
if isempty(ext)
    ext = '.png';
    filePath = [filePath ext];
end

try
    exportgraphics(fig, filePath, 'Resolution', 200);
catch
    % fallback
    print(fig, filePath, '-dpng', '-r200');
end
end

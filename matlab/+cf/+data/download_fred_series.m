function filePath = download_fred_series(series_id, outDir)
% Download a FRED series as CSV using fredgraph (no API key).
% Output CSV has DATE, value.

if ~exist(outDir,'dir'); mkdir(outDir); end
url = sprintf('https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s', series_id);

filePath = fullfile(outDir, sprintf('%s.csv', series_id));
websave(filePath, url);
end

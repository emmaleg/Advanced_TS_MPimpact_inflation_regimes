function fig = plot_posterior_hist(draws, figName, xLabel, varargin)
% Generic posterior histogram with nice limits.
%
% Parameters:
%   draws   : vector
%   figName : figure name / title
%   xLabel  : x-axis label
%
% Options:
%   'nbins'   (default 50)
%   'savePath' (default '')

p = inputParser;
p.addParameter('nbins',50,@(x)isnumeric(x) && isscalar(x) && x>=5);
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});
nbins = p.Results.nbins;
savePath = string(p.Results.savePath);

draws = draws(:);
draws = draws(isfinite(draws));

fig = figure('Name',char(figName),'Color','w');
histogram(draws, nbins);
grid on;
xlabel(xLabel);
ylabel('Count');
title(figName);

% padded x-limits
xmin = min(draws); xmax = max(draws);
if xmin < xmax
    pad = 0.05*(xmax - xmin);
    xlim([xmin - pad, xmax + pad]);
end

if strlength(savePath) > 0
    cf.plot.save_figure(fig, savePath);
end
end

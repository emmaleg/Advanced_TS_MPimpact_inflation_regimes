function fig = plot_posterior_discrete(draws, figName, xLabel, varargin)
% Posterior for discrete integer parameter (e.g., d*).
%
% Options:
%   'savePath' (default '')

p = inputParser;
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});
savePath = string(p.Results.savePath);

draws = draws(:);
draws = draws(isfinite(draws));

supp = unique(draws);
supp = sort(supp);

counts = zeros(size(supp));
for i=1:numel(supp)
    counts(i) = sum(draws == supp(i));
end
probs = counts / sum(counts);

fig = figure('Name',char(figName),'Color','w');
bar(supp, probs);
grid on;
xlabel(xLabel);
ylabel('Posterior probability');
title(figName);

xlim([min(supp)-0.5, max(supp)+0.5]);

if strlength(savePath) > 0
    cf.plot.save_figure(fig, savePath);
end
end

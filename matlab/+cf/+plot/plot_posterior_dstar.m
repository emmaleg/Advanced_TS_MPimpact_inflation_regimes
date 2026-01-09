function fig = plot_posterior_dstar(res, varargin)
p = inputParser;
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});

if ~isfield(res,'draws') || ~isfield(res.draws,'d')
    error("res.draws.d not found (d*).");
end

fig = cf.plot.plot_posterior_discrete(res.draws.d, 'Posterior of d*', 'd*', ...
    'savePath', p.Results.savePath);
end

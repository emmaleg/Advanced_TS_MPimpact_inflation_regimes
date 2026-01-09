function fig = plot_posterior_mu(res, varargin)
p = inputParser;
p.addParameter('nbins',50,@(x)isnumeric(x) && isscalar(x));
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});

if ~isfield(res,'draws') || ~isfield(res.draws,'mu')
    error("res.draws.mu not found (SV mean mu).");
end

fig = cf.plot.plot_posterior_hist(res.draws.mu, 'Posterior of \mu (SV mean)', '\mu', ...
    'nbins', p.Results.nbins, 'savePath', p.Results.savePath);
end

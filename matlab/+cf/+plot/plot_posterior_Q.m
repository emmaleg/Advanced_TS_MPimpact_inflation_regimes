function fig = plot_posterior_Q(res, varargin)
p = inputParser;
p.addParameter('nbins',50,@(x)isnumeric(x) && isscalar(x));
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});

if ~isfield(res,'draws') || ~isfield(res.draws,'Q')
    error("res.draws.Q not found (SV variance Q).");
end

fig = cf.plot.plot_posterior_hist(res.draws.Q, 'Posterior of Q (SV variance)', 'Q', ...
    'nbins', p.Results.nbins, 'savePath', p.Results.savePath);
end

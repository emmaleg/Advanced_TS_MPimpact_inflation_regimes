function fig = plot_posterior_F(res, varargin)
p = inputParser;
p.addParameter('nbins',50,@(x)isnumeric(x) && isscalar(x));
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});

if ~isfield(res,'draws') || ~isfield(res.draws,'F')
    error("res.draws.F not found (volatility persistence F).");
end

fig = cf.plot.plot_posterior_hist(res.draws.F, 'Posterior of F (SV persistence)', 'F', ...
    'nbins', p.Results.nbins, 'savePath', p.Results.savePath);
end

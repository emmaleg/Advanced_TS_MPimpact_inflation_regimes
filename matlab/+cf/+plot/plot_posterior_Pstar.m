function fig = plot_posterior_Pstar(res, varargin)
p = inputParser;
p.addParameter('nbins',50,@(x)isnumeric(x) && isscalar(x));
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.parse(varargin{:});

if ~isfield(res,'draws') || ~isfield(res.draws,'Pstar')
    error("res.draws.Pstar not found (threshold).");
end

fig = cf.plot.plot_posterior_hist(res.draws.Pstar, 'Posterior of P*', 'P*', ...
    'nbins', p.Results.nbins, 'savePath', p.Results.savePath);
end

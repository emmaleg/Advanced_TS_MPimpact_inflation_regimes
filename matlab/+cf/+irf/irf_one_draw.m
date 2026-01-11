function irf = irf_one_draw(ds, D, tStar, mconf, iconf, shockType)
% Approximate IRF = E[Z|shock] - E[Z|noshock] by Monte Carlo averaging (L reps).

Z = ds.Z;
[T,n] = size(Z);

p = iconf.p;
J = iconf.J;

% IMPORTANT: the regime indicator uses inflation at lag d (which can be >> p).
% If we only keep p lags of history, the first (d-p) simulated periods will
% not have enough history and is_low_regime() will fall back to S_low=true.
% That mechanically delays regime switching and creates "late" jumps.
% Solution: keep at least d lags in the initial history.
p0 = max(p, D.d);

% initial history for Z: use actual data tStar-p0+1 ... tStar
Zhist = Z(tStar-p0+1:tStar, :);

% initial history for log-lambda (need current and J lags)
% (D.h_path is already log(lambda))
h_all = D.h_path;

hHist = zeros(J+1,1);
for j=0:J
    hHist(j+1) = h_all(tStar - j);
end

IRFs = nan(iconf.H, n, iconf.L);

for l = 1:iconf.L
    % future innovations for structural shocks and for lambda AR(1)
    e_all   = randn(iconf.H, n);
    eta_all = randn(iconf.H, 1);

    [Z0, Zs] = cf.irf.simulate_pair(Zhist, hHist, D, mconf, iconf, shockType, e_all, eta_all);

    dZ = Zs - Z0; % H x n
    if all(isfinite(dZ(:)))
        IRFs(:,:,l) = dZ;
    end
end

irf = mean(IRFs, 3, 'omitnan');

end

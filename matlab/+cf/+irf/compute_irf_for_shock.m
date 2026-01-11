function res = compute_irf_for_shock(ds, post, mconf, iconf, shockType)
% One shock type: run Appendix-B style simulations, split by initial regime.

Z = ds.Z;
[T,n] = size(Z);

Nd = cf.irf.posterior_num_draws(post);
if Nd < 10
    error("Posterior draws look too few (Nd=%d). Check post struct.", Nd);
end

% --- default sampling policy ---
% Prefer WITHOUT replacement when feasible (reduces duplicates, more stable quantiles)
if ~isfield(iconf,'remise')
    iconf.remise = (Nd < iconf.S);  % if plenty of draws, default no-replacement
end

Sgoal = iconf.S;

% --- build valid tStar candidates (avoid NaNs + ensure lags exist) ---
p = mconf.p;
tMin = max([p, iconf.J, iconf.dmax]) + 1;

tCandidates = (tMin:T)';
% need inflation lag at least up to dmax
okInf = isfinite(Z(tCandidates - iconf.dmax, iconf.idx_inf));
tCandidates = tCandidates(okInf);

if isempty(tCandidates)
    error("No valid tStar candidates after filtering. Check data / idx_inf / dmax.");
end

% --- draw order for posterior draws ---
draw_order = [];
draw_ptr   = 1;
if ~iconf.remise
    % try to use as many unique draws as possible
    draw_order = randperm(Nd, Nd);
end

% store raw IRFs split by initial regime
irf_low  = [];
irf_high = [];

kept = 0;
tries = 0;

% tries budget
maxTries = max(Sgoal*5, Sgoal+200);     % good default even with many rejections
maxTries = min(maxTries, 5*Nd + 500);   % avoid pathological loops when Nd huge

while kept < Sgoal && tries < maxTries
    tries = tries + 1;

    % --- draw index (with fallback) ---
    if iconf.remise
        kdraw = randi(Nd);
    else
        if draw_ptr > Nd
            % fallback to with replacement to reach Sgoal if too many rejections
            iconf.remise = true;
            kdraw = randi(Nd);
        else
            kdraw = draw_order(draw_ptr);
            draw_ptr = draw_ptr + 1;
        end
    end

    % --- random history draw ---
    tStar = tCandidates(randi(numel(tCandidates)));

    D = cf.irf.get_draw(post, kdraw, T, n, iconf, mconf);

    % 1) stability filter (VAR part only)
    [st1, ~] = cf.irf.is_stable_var(D.Phi1, n, p);
    [st2, ~] = cf.irf.is_stable_var(D.Phi2, n, p);
    if ~(st1 && st2)
        continue
    end

    % initial regime (eq. (2))
    inf_lag = Z(tStar - D.d, iconf.idx_inf);
    if ~isfinite(inf_lag)
        continue
    end
    S0_low = (inf_lag <= D.Pstar); % S=1 = low inflation regime in your convention

    irf = cf.irf.irf_one_draw(ds, D, tStar, mconf, iconf, shockType); % H x n

    % 2) numerical guard
    if any(~isfinite(irf(:)))
        continue
    end

    % (optional) magnitude guard (prevents 1–2 crazy draws flattening quantiles)
    if isfield(iconf,'max_abs_irf') && ~isempty(iconf.max_abs_irf)
        if max(abs(irf(:))) > iconf.max_abs_irf
            continue
        end
    end

    if S0_low
        irf_low  = cat(3, irf_low,  irf);
    else
        irf_high = cat(3, irf_high, irf);
    end

    kept = kept + 1;

    if mod(kept,50)==0
        fprintf("[irf:%s] %d/%d kept (tries=%d, low=%d, high=%d)\n", ...
            shockType, kept, Sgoal, tries, size(irf_low,3), size(irf_high,3));
    end
end

if kept < Sgoal
    warning("compute_irf_for_shock: kept only %d/%d draws after %d tries.", kept, Sgoal, tries);
end

% warn if one regime is too thin
if size(irf_low,3) < 30
    warning("[irf:%s] Low-regime IRFs based on only %d draws -> quantiles may be noisy.", shockType, size(irf_low,3));
end
if size(irf_high,3) < 30
    warning("[irf:%s] High-regime IRFs based on only %d draws -> quantiles may be noisy.", shockType, size(irf_high,3));
end

res = struct();
res.h = (0:iconf.H-1)';

res.low  = cf.irf.summarize_irfs(irf_low);
res.high = cf.irf.summarize_irfs(irf_high);

res.shockType = shockType;
end



function res = compute_irf_for_shock(ds, post, mconf, iconf, shockType)
% One shock type: run Appendix-B style simulations, split by initial regime.

Z = ds.Z;
[T,n] = size(Z);

Nd = cf.irf.posterior_num_draws(post);
if Nd < 10
    error("Posterior draws look too few (Nd=%d). Check post struct.", Nd);
end

% --- default sampling policy ---
if ~isfield(iconf,'remise')
    iconf.remise = (Nd < iconf.S);  
end

Sgoal = iconf.S;

% =========================================================================
% CONSTRUCTION DES CANDIDATS (CORRIGÉ)
% =========================================================================
p = mconf.p;

% 1. Définition de tMin (C'est la ligne qui vous manquait)
% On doit commencer après le max des lags (p, J) et du délai (dmax)
tMin = max([p, iconf.J, iconf.dmax]) + 1;

% 2. Candidats par défaut (tout l'échantillon valide)
tCandidates = (tMin:T)';

% 3. FILTRE "SUBSET" (Pour Volcker, etc.)
if isfield(iconf, 'subset_indices') && ~isempty(iconf.subset_indices)
    % On ne garde que l'intersection entre les dates mathématiquement valides
    % et les dates demandées par l'utilisateur
    tCandidates = intersect(tCandidates, iconf.subset_indices);
    
    if isempty(tCandidates)
        % On renvoie une structure vide propre au lieu de planter
        warning("Aucune date valide trouvée dans la période demandée !");
        res = struct();
        res.h = (0:iconf.H-1)';
        res.low  = cf.irf.summarize_irfs([]);
        res.high = cf.irf.summarize_irfs([]);
        res.shockType = shockType;
        return;
    end
    fprintf('[compute_irf:%s] Restriction active : %d dates candidates retenues.\n', shockType, numel(tCandidates));
end
% =========================================================================

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
    draw_order = randperm(Nd, Nd);
end

% store raw IRFs split by initial regime
irf_low  = [];
irf_high = [];

kept = 0;
tries = 0;

% tries budget
maxTries = max(Sgoal*5, Sgoal+200);     
maxTries = min(maxTries, 5*Nd + 500);   

while kept < Sgoal && tries < maxTries
    tries = tries + 1;

    % --- draw index ---
    if iconf.remise
        kdraw = randi(Nd);
    else
        if draw_ptr > Nd
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

    % 0) identification filter (apply impact restrictions in BOTH regimes)
    if isfield(mconf,'id') && isfield(mconf.id,'enforce_sign_zero') && mconf.id.enforce_sign_zero
        tol = mconf.id.zero_tol;
        lam = exp(mean(D.h_path));         % scalar >0

        % check_impact_restrictions expects sigma2, so pass D.sigma.^2
        okA1 = cf.id.check_impact_restrictions(D.A1, (D.sigma.^2), lam, tol);
        okA2 = cf.id.check_impact_restrictions(D.A2, (D.sigma.^2), lam, tol);

        if ~(okA1 && okA2)
            continue
        end
    end

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
    S0_low = (inf_lag <= D.Pstar); 

    irf = cf.irf.irf_one_draw(ds, D, tStar, mconf, iconf, shockType); 

    % 2) numerical guard
    if any(~isfinite(irf(:)))
        continue
    end

    % (optional) magnitude guard
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

res = struct();
res.h = (0:iconf.H-1)';
res.low  = cf.irf.summarize_irfs(irf_low);
res.high = cf.irf.summarize_irfs(irf_high);
res.shockType = shockType;

end
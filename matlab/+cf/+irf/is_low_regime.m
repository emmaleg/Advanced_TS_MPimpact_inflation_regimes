function S_low = is_low_regime(Zfull, idx, Pstar, d, idx_inf)
% S_t = 1 <=> Pi_{t-d} <= P*  (low inflation regime)
j = idx - d;
if j < 1
    S_low = true; % fallback (should not happen if tStar chosen safely)
    return
end
S_low = (Zfull(j, idx_inf) <= Pstar);
end

function sgn = normalize_sign_liquidity(B, iconf)
% Want: M2 impact > 0 (expansionary liquidity). If zero-ish, use slope <= 0.
imp_m2 = B(iconf.idx_m2, iconf.shock_liquidity);
sgn = 1;
if imp_m2 < 0
    sgn = -1;
elseif abs(imp_m2) < 1e-10
    imp_slope = B(iconf.idx_slope, iconf.shock_liquidity);
    if imp_slope > 0
        sgn = -1;
    end
end
end

function sgn = normalize_sign_conventional(B, iconf)
% Want: FFR impact > 0 (contractionary); and typically M2 impact < 0.
imp_ffr = B(iconf.idx_ffr, iconf.shock_conventional);
sgn = 1;
if imp_ffr < 0
    sgn = -1;
end
end

function e = enforce_ffr(B, muZ, e, r_desired, iconf)
% Adjust control shock to match desired FFR in this period:
% r = mu_r + B_r * e  => solve for e(control)

r_idx = iconf.idx_ffr;
jctrl = iconf.shock_ffr_control;

br = B(r_idx, jctrl);
if abs(br) < 1e-12
    error("Cannot enforce FFR: B(ffr,control) is ~0.");
end

r_hat = muZ(r_idx) + B(r_idx,:)*e;
e(jctrl) = e(jctrl) + (r_desired - r_hat)/br;

end

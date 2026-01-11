function out = compute_irfs(ds, post, mconf, iconf)
%CF.IRF.COMPUTE_IRFS  Compute regime-averaged IRFs per Appendix B.
%
% out.conventional.low/high  : structs with median/p16/p84, size H x n
% out.liquidity.low/high     : idem
%
% Requirements:
% - ds.Z (T x n), ds.dates, ds.names
% - post contains posterior draws for: Pstar, d, Phi1, Phi2, alpha1, alpha2,
%   sigma2, lambda_path (T x Nd), mu, F, Q

rng(iconf.seed);

out = struct();
out.meta = struct('H',iconf.H,'S',iconf.S,'L',iconf.L,'delta',iconf.delta, ...
                  'ffr_lock_h',iconf.ffr_lock_h);

out.conventional = cf.irf.compute_irf_for_shock(ds, post, mconf, iconf, "conventional");
out.liquidity    = cf.irf.compute_irf_for_shock(ds, post, mconf, iconf, "liquidity");

end

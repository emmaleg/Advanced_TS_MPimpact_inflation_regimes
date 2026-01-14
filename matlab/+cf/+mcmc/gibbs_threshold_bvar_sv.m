function res = gibbs_threshold_bvar_sv(ds, mconf, pconf, mcmc)
% Gibbs sampler for Threshold-BVAR with scalar lambda_t SV.
% Blocks follow Appendix A structure.

rng(mcmc.seed);

Z = ds.Z;
[T,n] = size(Z);

% Initialize
st = cf.mcmc.init_state(Z, ds.dates, mconf, pconf, mcmc);

% Storage
nkeep = floor((mcmc.K - mcmc.burn)/mcmc.thin);
draws.Pstar = NaN(nkeep,1);
draws.d     = NaN(nkeep,1);

kPhi = 1 + n*mconf.p + (mconf.J+1);
if mcmc.store_full_draws
    draws.Phi1 = NaN(n, kPhi, nkeep);
    draws.Phi2 = NaN(n, kPhi, nkeep);
    draws.alpha1 = NaN(mconf.id.n_alpha, nkeep);
    draws.alpha2 = NaN(mconf.id.n_alpha, nkeep);
    draws.sigma2 = NaN(n, nkeep);
    draws.mu = NaN(nkeep,1);
    draws.F  = NaN(nkeep,1);
    draws.Q  = NaN(nkeep,1);
    if mcmc.store_lambda_path
        draws.h = NaN(T, nkeep);
    end
end

S_sum = zeros(T,1);
keep_idx = 0;

Pi = Z(:, mconf.inflation_index_in_Z);

for k=1:mcmc.K
    % --- Update regime indicator from current (P*,d) (KEEP Svalid!)
    [st.S, st.Svalid] = cf.model.regime_indicator(Pi, st.Pstar, st.d);

    % 1) P* | rest
    st = cf.mcmc.draw_Pstar_mh(st, Z, mconf, pconf, mcmc, k);
    [st.S, st.Svalid] = cf.model.regime_indicator(Pi, st.Pstar, st.d);

    % 2) d | rest
    st = cf.mcmc.draw_d_multinomial(st, Z, mconf, pconf);
    % draw_d_multinomial already updates st.S and st.Svalid

    % 3) Phi_i | rest
    st = cf.mcmc.draw_phi_sur(st, Z, mconf, pconf);

    % 4) alpha_i | rest
    st = cf.mcmc.draw_alpha_mh(st, Z, mconf, pconf, mcmc);

    % 5) sigma2 | rest
    st = cf.mcmc.draw_sigma2_ig(st, Z, mconf, pconf);

    % 6) h | rest
    st = cf.mcmc.draw_lambda_single_move(st, Z, mconf, pconf, mcmc);

    % 7) mu | rest
    st = cf.mcmc.draw_mu_normal(st, pconf);

    % 8) F | rest
    st = cf.mcmc.draw_F_truncnorm(st, pconf);

    % 9) Q | rest
    st = cf.mcmc.draw_Q_ig(st, pconf);

    % Store
    if k > mcmc.burn && mod(k-mcmc.burn, mcmc.thin)==0
        keep_idx = keep_idx + 1;

        draws.Pstar(keep_idx) = st.Pstar;
        draws.d(keep_idx)     = st.d;

        S_sum = S_sum + st.S;

        if mcmc.store_full_draws
            draws.Phi1(:,:,keep_idx) = st.Phi1;
            draws.Phi2(:,:,keep_idx) = st.Phi2;
            draws.alpha1(:,keep_idx) = st.alpha1;
            draws.alpha2(:,keep_idx) = st.alpha2;
            draws.sigma2(:,keep_idx) = st.sigma2;
            draws.mu(keep_idx) = st.mu;
            draws.F(keep_idx)  = st.F;
            draws.Q(keep_idx)  = st.Q;
            if mcmc.store_lambda_path
                draws.h(:,keep_idx) = st.h;
            end
        end
    end

    if mod(k,100)==0
        fprintf('[Gibbs] iter %d / %d (kept=%d)\n', k, mcmc.K, keep_idx);
    end
end

res = struct();
res.draws  = draws;
res.S_mean = S_sum / nkeep;
res.accept = st.accept;

end

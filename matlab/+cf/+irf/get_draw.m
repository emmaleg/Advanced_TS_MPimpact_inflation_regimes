function D = get_draw(post, k, T, n, iconf, mconf)
% Extract posterior draw k from THIS repo format (post = res.draws).

D = struct();

% --- scalars ---
D.Pstar = post.Pstar(k);
D.d     = round(post.d(k));
D.d     = max(1, min(iconf.dmax, D.d));

% --- Phi matrices ---
K = 1 + n*mconf.p + (iconf.J+1);
D.Phi1 = post.Phi1(:,:,k);
D.Phi2 = post.Phi2(:,:,k);

% (optionnel, juste pour être sûr)
if ~isequal(size(D.Phi1), [n, K]) || ~isequal(size(D.Phi2), [n, K])
    error("Phi size mismatch: expected [%d x %d]. Got Phi1=%s Phi2=%s", ...
        n, K, mat2str(size(D.Phi1)), mat2str(size(D.Phi2)));
end

% --- alpha vectors ---
D.alpha1 = post.alpha1(:,k);
D.alpha2 = post.alpha2(:,k);

% Build A matrices
D.A1 = cf.id.alpha_to_A(D.alpha1);
D.A2 = cf.id.alpha_to_A(D.alpha2);

% --- sigma2 (n x 1) ---
sig2 = post.sigma2(:,k);
if numel(sig2) ~= n
    error("sigma2 draw has length %d but n=%d", numel(sig2), n);
end
D.sigma = sqrt(sig2(:)); % std devs

% --- h path (T x 1), where h = log(lambda) ---
if ~isfield(post,'h')
    error("post.h missing. Set mcmc.store_full_draws=true and mcmc.store_lambda_path=true in Gibbs.");
end
if size(post.h,1) ~= T
    error("post.h has wrong T dimension: size(post.h)=%s, expected first dim T=%d", ...
        mat2str(size(post.h)), T);
end

D.h_path      = post.h(:,k);
D.lambda_path = exp(D.h_path);

% --- AR(1) hyperparams for h ---
D.mu = post.mu(k);
D.F  = post.F(k);
D.Q  = post.Q(k);

% sanity bounds (as before)
D.F = max(0.001, min(0.999, D.F));
D.Q = max(1e-10, D.Q);

end

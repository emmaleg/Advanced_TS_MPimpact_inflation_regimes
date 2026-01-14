function ll = loglik_obs(z, x, Phi, A, sigma2, h_t)
% One-observation log-likelihood under:
% eps_t = z - Phi x
% u_t   = A * eps_t
% u_t ~ N(0, lambda_t * diag(sigma2)), lambda_t = exp(h_t)

eps = z - (Phi * x);
u   = A * eps;

lam = exp(h_t);
v   = lam * sigma2(:);  % variances

% log N(0, diag(v))
ll = -0.5 * ( numel(u)*log(2*pi) + sum(log(v)) + sum((u.^2) ./ v) );
end

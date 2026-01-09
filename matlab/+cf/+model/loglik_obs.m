function ll = loglik_obs(z_t, x_t, Phi, A, sigma2, h_t)
% log p(z_t | ...) under:
% eps_t = z_t - Phi * x_t
% Omega_t = A^{-1} * (lambda_t * Sigma) * A^{-1}'
% with Sigma = diag(sigma2), lambda_t = exp(h_t).
n = numel(z_t);
Sigma = diag(sigma2(:));
lambda = exp(h_t);

eps = z_t - (Phi * x_t);

% Use transformation u = A * eps, u ~ N(0, lambda*Sigma)
u = A * eps;

% logdet(H) = n*log(lambda) + logdet(Sigma)
logdetH = n*log(lambda) + sum(log(sigma2));

% Jacobian term log|A|
logdetA = log(abs(det(A)));

quad = (u' * ( (1/lambda) * (Sigma \ u) ));

ll = logdetA - 0.5*logdetH - 0.5*quad - 0.5*n*log(2*pi);
end

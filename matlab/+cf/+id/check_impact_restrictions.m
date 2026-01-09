function ok = check_impact_restrictions(A, sigma2, lambda, zero_tol)
% Optional: check (very roughly) Table-1-style sign/zero restrictions on impact matrix
% Impact = A^{-1} * chol(lambda*diag(sigma2)) (diag positive -> sign from A^{-1})
% WARNING: depending on parameterization, "exact zero" may not happen numerically.

if nargin<4, zero_tol=1e-10; end

B0 = inv(A);
% Conventional shock = column 4 ; Liquidity shock = column 6 (as noted in appendix)
c4 = B0(:,4);
c6 = B0(:,6);

% zeros: IP (1), Unemp (3) for conventional; IP(1), Infl(2), Unemp(3), Rate(4) for liquidity
ok = true;

ok = ok && (abs(c4(1))<=zero_tol) && (abs(c4(3))<=zero_tol);
ok = ok && (c4(4) > 0) && (c4(6) < 0) && (c4(2) <= 0);

ok = ok && (abs(c6(1))<=zero_tol) && (abs(c6(2))<=zero_tol) && (abs(c6(3))<=zero_tol) && (abs(c6(4))<=zero_tol);
ok = ok && (c6(5) <= 0) && (c6(6) > 0);

end

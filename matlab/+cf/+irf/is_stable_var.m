function [isStable, rho] = is_stable_var(Phi, n, p)
%IS_STABLE_VAR  Check stability of the VAR lag polynomial.
%
% Phi is n x K with columns ordered as:
% [const, Z_{t-1} (n), ..., Z_{t-p} (n), exog...]
%
% We build the companion matrix using only the VAR lag blocks and check
% that all eigenvalues are strictly inside the unit circle.

% extract lag coefficient blocks
A = cell(1,p);
for j = 1:p
    cols = (2 + (j-1)*n) : (1 + j*n);
    A{j} = Phi(:, cols);
end

if p == 1
    C = A{1};
else
    top = [A{:}];
    low = [eye(n*(p-1)), zeros(n*(p-1), n)];
    C = [top; low];
end

ev = eig(C);
rho = max(abs(ev));

% strict check (slightly < 1 to avoid nearly-explosive draws)
isStable = isfinite(rho) && (rho < 0.999);
end
function ok = check_impact_restrictions(A, sigma2, lambda, zero_tol)
% Check Table-1-style sign/zero restrictions on the contemporaneous impact matrix.
%
% Model uses: Impact B = A^{-1} * sqrt(lambda) * diag(sigma)
% where sigma = sqrt(sigma2) and sqrt(lambda)>0 so signs are inherited from A^{-1}.
%
% IMPORTANT: exact zeros are numerically impossible -> use relative tolerance.

if nargin < 4 || isempty(zero_tol)
    zero_tol = 1e-3;  % relative tolerance (you can tighten/relax)
end

n = size(A,1);
if n < 6
    ok = false;
    return;
end

% Build impact matrix B (scale matters for "near zero")
sigma = sqrt(sigma2(:));
D     = diag(sqrt(max(lambda,1e-12)) * sigma);  % n x n, positive diagonal
B     = A \ D;                                  % n x n

% Shocks columns as in your code/paper convention
c4 = B(:,4);   % conventional MP shock
c6 = B(:,6);   % liquidity shock

sc4 = max(1, norm(c4, inf));
sc6 = max(1, norm(c6, inf));

isZero = @(x,sc) abs(x) <= zero_tol * sc;

ok = true;

% ----- Conventional shock (#4) -----
% Zero impact on IP and UNRATE (your current convention)
ok = ok && isZero(c4(1), sc4);
ok = ok && isZero(c4(3), sc4);

% Sign restrictions (as in your earlier implementation)
ok = ok && (c4(4) > 0);   % FFR up (contractionary)
ok = ok && (c4(2) < 0);   % inflation down
ok = ok && (c4(6) < 0);   % M2 down

% ----- Liquidity shock (#6) -----
% Zero impact on IP, inflation, unemployment, and FFR contemporaneously
ok = ok && isZero(c6(1), sc6);
ok = ok && isZero(c6(2), sc6);
ok = ok && isZero(c6(3), sc6);
ok = ok && isZero(c6(4), sc6);

% Signs: slope down, M2 up
ok = ok && (c6(5) < 0);
ok = ok && (c6(6) > 0);

% ----- Extra safeguard for IRF lock feasibility -----
% If |B(FFR, conventional shock)| is tiny, enforcing the FFR lock explodes.
minFFR = max(1e-6, 1e-3 * sc4);
ok = ok && isfinite(B(4,4)) && abs(B(4,4)) >= minFFR;

end


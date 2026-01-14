function ok = check_impact_restrictions(A, sigma2, lambda, zero_tol)
% Check Table-1 identification restrictions (Canova & Forero, 2025).
% Impact matrix: B = A^{-1} * sqrt(lambda) * diag(sigma)
% where sigma = sqrt(sigma2). Since sqrt(lambda)>0, it does not affect signs.
%
% Table 1 (impact restrictions):
% Conventional shock (col 4):
%   IP=0, Inflation=0, Unemp=0, Rate>0, Money<0
% Liquidity shock (col 6):
%   IP=0, Inflation=0, Unemp=0, Rate=0, Slope<=0, Money>0
%
% zero_tol is RELATIVE to the column scale.

if nargin < 4 || isempty(zero_tol)
    zero_tol = 1e-3;
end

n = size(A,1);
if n < 6
    ok = false;
    return;
end

sigma  = sqrt(sigma2(:));
lamPos = max(lambda, 1e-12);
D      = diag(sqrt(lamPos) * sigma);  % positive diagonal
B      = A \ D;                       % impact matrix

c4 = B(:,4);   % conventional MP shock
c6 = B(:,6);   % liquidity shock

sc4 = max(1, norm(c4, inf));
sc6 = max(1, norm(c6, inf));
isZero = @(x,sc) abs(x) <= zero_tol * sc;

ok = true;

% ---------------------------
% Conventional shock (#4)
% ---------------------------
ok = ok && isZero(c4(1), sc4);    % IP growth: 0
ok = ok && isZero(c4(2), sc4);    % Inflation: 0
ok = ok && isZero(c4(3), sc4);    % Unemployment: 0
ok = ok && (c4(4) > 0);           % Short rate: >0
ok = ok && (c4(6) < 0);           % Money growth: <0

% ---------------------------
% Liquidity shock (#6)
% ---------------------------
ok = ok && isZero(c6(1), sc6);    % IP growth: 0
ok = ok && isZero(c6(2), sc6);    % Inflation: 0
ok = ok && isZero(c6(3), sc6);    % Unemployment: 0
ok = ok && isZero(c6(4), sc6);    % Short rate: 0 (impact)
ok = ok && (c6(5) <= 0);          % Slope: <=0
ok = ok && (c6(6) > 0);           % Money growth: >0

% Practical safeguard for IRF lock feasibility (optional but useful):
minFFR = max(1e-8, 1e-4 * sc4);
ok = ok && isfinite(B(4,4)) && abs(B(4,4)) >= minFFR;

end

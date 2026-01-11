function [S, Svalid] = regime_indicator(Pi, Pstar, d)
% Regime indicator:
%   S_t = 1  <=>  Pi_{t-d} <= P*     (low inflation regime)
% Valid only for t > d. We return:
%   S      : binary 0/1 for ALL t (filled at the beginning for convenience)
%   Svalid : logical mask, true only where S_t is truly defined (t>d)

T = numel(Pi);
S = zeros(T,1);
Svalid = false(T,1);

if d < 0 || d ~= round(d)
    error('regime_indicator: d must be a nonnegative integer.');
end

if d == 0
    S(:) = double(Pi(:) <= Pstar);
    Svalid(:) = true;
    return;
end

% Defined part
t = (d+1):T;
S(t) = double(Pi(t-d) <= Pstar);
Svalid(t) = true;

% Fill initial undefined periods for plotting / code robustness (won't matter if Svalid is used)
S(1:d) = S(d+1);
end


function [S, Svalid] = regime_indicator(Pi, Pstar, d)
%CF.MODEL.REGIME_INDICATOR
% S_t = 1 (low inflation) iff Pi_{t-d} <= Pstar, defined only for t>d.

T = numel(Pi);
S      = NaN(T,1);
Svalid = false(T,1);

for t = 1:T
    if t > d && isfinite(Pi(t-d))
        Svalid(t) = true;
        S(t)      = double(Pi(t-d) <= Pstar); % 1=low, 0=high
    end
end
end


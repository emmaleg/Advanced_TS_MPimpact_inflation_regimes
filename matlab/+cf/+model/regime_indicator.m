function S = regime_indicator(Pi, Pstar, d)
% S_t = 1 iff Pi_{t-d} <= P*   (low inflation regime)
T = numel(Pi);
S = ones(T,1);
for t=1:T
    td = t - d;
    if td >= 1
        S(t) = double(Pi(td) <= Pstar);
    else
        S(t) = 1; % default low regime for initial periods
    end
end
end

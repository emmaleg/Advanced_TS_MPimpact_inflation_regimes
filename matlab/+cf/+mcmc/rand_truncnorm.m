function x = rand_truncnorm(mu, sd, lb, ub)
% Simple rejection sampler for truncated normal N(mu,sd^2) on [lb,ub].

while true
    x = mu + sd*randn();
    if x >= lb && x <= ub
        return;
    end
end
end

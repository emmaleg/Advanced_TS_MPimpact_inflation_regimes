function x = rand_ig(a, b)
% Inverse-Gamma(a,b) with density proportional x^(-a-1) exp(-b/x).
% Sample via: 1/x ~ Gamma(a, 1/b)

x = 1 / gamrnd(a, 1/b);
end

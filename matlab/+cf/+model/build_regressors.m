function [Y, X] = build_regressors(Z, h, p, J)
% Build regression matrices for t=p+1..T:
% Z_t = c + sum_{j=1}^p B_j Z_{t-j} + sum_{j=0}^J G_j h_{t-j} + eps
% where h_t = ln(lambda_t).

[T,n] = size(Z);
Xlag = cf.model.make_lag_matrix(Z, p);

% Build h lags (including j=0)
Hlags = NaN(T, J+1);
for j=0:J
    Hlags(:, j+1) = lagmatrix(h, j);
end

% Design matrix with constant
Xfull = [ones(T,1), Xlag, Hlags];

% Keep only rows where all regressors finite and t>=p+1
start = p+1;
Y = Z(start:end,:);
X = Xfull(start:end,:);

% Drop rows with NaNs (from h lags)
good = all(isfinite(X),2) & all(isfinite(Y),2);
Y = Y(good,:);
X = X(good,:);
end

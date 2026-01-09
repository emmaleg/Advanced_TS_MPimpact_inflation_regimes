function Xlag = make_lag_matrix(Z, p)
% Returns [Z_{t-1},...,Z_{t-p}] for each t (rows aligned with Z).
[T,n] = size(Z);
Xlag = NaN(T, n*p);
for j=1:p
    Xlag(:, (n*(j-1)+1):(n*j)) = lagmatrix(Z, j);
end
end

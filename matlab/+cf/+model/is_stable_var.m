function ok = is_stable_var(Phi, n, p)
% Check VAR(p) stability using companion matrix eigenvalues.
% Phi: n x k, where k = 1 + n*p + (J+1). We only use AR part.

% Extract AR blocks (ignore constant and h-lags)
% columns: [const | Z_{t-1}..Z_{t-p} | h_{t}..h_{t-J}]
Acols_start = 2;
Acols_end   = 1 + n*p;
A = Phi(:, Acols_start:Acols_end); % n x (n*p)

% Companion
Comp = zeros(n*p);
Comp(1:n,:) = A;
Comp(n+1:end,1:n*(p-1)) = eye(n*(p-1));

ev = eig(Comp);
ok = all(abs(ev) < 1-1e-10);
end

function x = make_xt(Zfull, idx, p, h_path, h0, t, J)
% x_t = [1; Z_{t-1}; ...; Z_{t-p}; h_t; h_{t-1}; ...; h_{t-J}]

n = size(Zfull,2);

% VAR lags
lags = [];
for j=1:p
    lags = [lags; Zfull(idx-j, :)']; %#ok<AGROW>
end

% h lags: current h is h_path(t); then h_{t-1}, h_{t-2}...
h = zeros(J+1,1);
h(1) = h_path(t);
for j=1:J
    if t-j >= 1
        h(j+1) = h_path(t-j);
    else
        % before simulation start, use initial history h0=[h_0,h_-1,h_-2]
        h(j+1) = h0(j+1);
    end
end

x = [1; lags; h];
x = reshape(x, [], 1);

% quick sanity
K = 1 + n*p + (J+1);
if numel(x) ~= K
    error("make_xt produced %d elements, expected %d", numel(x), K);
end
end


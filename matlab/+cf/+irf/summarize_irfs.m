function s = summarize_irfs(irf3d)
% irf3d: H x n x Ndraws
if isempty(irf3d)
    s = struct('median',[],'p16',[],'p84',[],'N',0);
    return
end

H = size(irf3d,1);
n = size(irf3d,2);
N = size(irf3d,3);

med = zeros(H,n);
p16 = zeros(H,n);
p84 = zeros(H,n);

for i=1:H
    for j=1:n
        x = squeeze(irf3d(i,j,:));
        x = x(isfinite(x));
        if isempty(x)
            p16(i,j) = NaN;
            med(i,j) = NaN;
            p84(i,j) = NaN;
        else
            q = prctile(x, [16 50 84]);
            p16(i,j) = q(1);
            med(i,j) = q(2);
            p84(i,j) = q(3);
        end
    end
end

s = struct('median', med, 'p16', p16, 'p84', p84, 'N', N);

end

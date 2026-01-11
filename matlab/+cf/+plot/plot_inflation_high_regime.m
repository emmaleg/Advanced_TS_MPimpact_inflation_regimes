function fig = plot_inflation_high_regime(ds, mconf, res, varargin)
% Plot inflation (left axis) + high-regime indicator (right axis)
%
% Usage:
%   fig = cf.plot.plot_inflation_high_regime(ds,mconf,res,'savePath',...,'highRule','<=0.5')

p = inputParser;
p.addParameter('savePath','',@(s)ischar(s) || isstring(s));
p.addParameter('highRule','<=0.5',@(s)ischar(s) || isstring(s)); % '<=0.5' or '>=0.5'
p.parse(varargin{:});
savePath = string(p.Results.savePath);
highRule = string(p.Results.highRule);

% --- regime prob / indicator
if isfield(res,'S_mean')
    Sm = res.S_mean;
elseif isfield(res,'draws') && isfield(res.draws,'S')
    Sm = mean(res.draws.S, 2);
else
    error('S_mean (or draws.S) not found in res.');
end

switch highRule
    case "<=0.5"
        high = (Sm <= 0.5);
    case ">=0.5"
        high = (Sm >= 0.5);
    otherwise
        error("highRule must be '<=0.5' or '>=0.5'.");
end

% --- figure
fig = figure('Name','Inflation & High-regime indicator','Color','w');
ax = axes(fig); 

yyaxis left;
plot(ds.dates, ds.Z(:, mconf.inflation_index_in_Z), 'LineWidth', 1.2);
ylabel('Inflation');
grid on;

yyaxis right;
stem(ds.dates(high), ones(nnz(high),1), 'filled');
ylim([0 1.2]);
ylabel('High regime (indicator)');

title('Inflation and high-regime indicator');

% x-lims with padding (a bit before/after sample)
d1 = ds.dates(1);
dT = ds.dates(end);
try
    xlim([d1 - calmonths(2), dT + calmonths(2)]);
catch
    % fallback if dates are not datetime
end

if strlength(savePath) > 0
    cf.plot.save_figure(fig, savePath);
end
end

function make_fig4(paths)
% make_fig4: "Figure 4-like" posterior distribution of P*

if ~exist(fullfile(paths.results_out,'figures'),'dir')
    mkdir(fullfile(paths.results_out,'figures'));
end

load(fullfile(paths.results_out,'posterior_draws.mat'),'res');

Pstar = res.draws.Pstar;
Pmed = median(Pstar);
P16 = prctile(Pstar,16);
P84 = prctile(Pstar,84);

f = figure('Name','Fig4-like','Visible','off');

% Histogram
histogram(Pstar, 60, 'Normalization','pdf');
hold on;

% Kernel density
[xi,fi] = ksdensity(Pstar);
plot(xi, fi, 'LineWidth', 1.5);

xline(P16, '--', 'LineWidth', 1.1);
xline(Pmed,'--', 'LineWidth', 1.5);
xline(P84, '--', 'LineWidth', 1.1);

title('Posterior distribution of inflation threshold P^*');
xlabel('P^*'); ylabel('Density');
grid on;
set(gca,'Box','off');

out = fullfile(paths.results_out,'figures','fig4.pdf');
exportgraphics(f, out, 'ContentType','vector');
close(f);

fprintf('[make_fig4] Saved: %s\n', out);
end

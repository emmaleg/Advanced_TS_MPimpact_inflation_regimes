function make_fig3(paths)
% make_fig3: "Figure 3-like" inflation with high-regime bars
% High regime shown when posterior mean of low-regime indicator <= 0.5.

if ~exist(fullfile(paths.results_out,'figures'),'dir')
    mkdir(fullfile(paths.results_out,'figures'));
end

load(fullfile(paths.data_out,'dataset_canova_forero.mat'),'ds');
load(fullfile(paths.results_out,'posterior_draws.mat'),'res');

Pi = ds.Z(:,2);
Smean_low = res.S_mean;

high = (Smean_low <= 0.5);

f = figure('Name','Fig3-like','Visible','off');

yyaxis left;
plot(ds.dates, Pi, 'LineWidth', 1.2);
ylabel('Inflation \pi_t (YoY log diff)');
grid on;

yyaxis right;
stem(ds.dates(high), ones(sum(high),1), 'filled');
ylim([0 1.2]);
ylabel('High regime indicator (mean prob > 0.5)');
title('Inflation and inferred high-inflation regime');
set(gca,'Box','off');

out = fullfile(paths.results_out,'figures','fig3.pdf');
exportgraphics(f, out, 'ContentType','vector');
close(f);

fprintf('[make_fig3] Saved: %s\n', out);
end

function make_fig2(paths)
% make_fig2: "Figure 2-like" diagnostic plot
% 1) time series of inflation (P in Z)
% 2) histogram of inflation with posterior median P* overlay

if ~exist(fullfile(paths.results_out,'figures'),'dir')
    mkdir(fullfile(paths.results_out,'figures'));
end

load(fullfile(paths.data_out,'dataset_canova_forero.mat'),'ds');
load(fullfile(paths.results_out,'posterior_draws.mat'),'res');

Pi = ds.Z(:,2); % inflation is 2nd variable by construction

Pstar_draws = res.draws.Pstar;
Pstar_med = median(Pstar_draws);

% --- Panel A: inflation time series ---
f1 = figure('Name','Fig2-like: Inflation time series','Visible','off');
plot(ds.dates, Pi, 'LineWidth', 1.2);
hold on;
yline(Pstar_med, '--', 'LineWidth', 1.2);
title('Inflation (YoY log diff) and posterior median P^*');
xlabel('Date'); ylabel('\pi_t');
grid on;
set(gca,'Box','off');

out1 = fullfile(paths.results_out,'figures','fig2_panelA_inflation.pdf');
exportgraphics(f1, out1, 'ContentType','vector');
close(f1);

% --- Panel B: inflation histogram + P* median ---
f2 = figure('Name','Fig2-like: Inflation histogram','Visible','off');
histogram(Pi, 60);
hold on;
xline(Pstar_med, '--', 'LineWidth', 1.2);
title('Distribution of inflation and posterior median P^*');
xlabel('\pi_t'); ylabel('Count');
grid on;
set(gca,'Box','off');

out2 = fullfile(paths.results_out,'figures','fig2_panelB_hist.pdf');
exportgraphics(f2, out2, 'ContentType','vector');
close(f2);

fprintf('[make_fig2] Saved:\n  %s\n  %s\n', out1, out2);
end

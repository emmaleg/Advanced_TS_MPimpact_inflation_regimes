%% run_02c_diag_sv.m
% Diagnostics SV (Jacquier single-move): acceptance, traceplots, lambda_t summaries.
% Run after run_02_gibbs.m (loads output/mcmc/posterior_draws.mat).

clear; clc; close all;

% --- root/matlab detection (same style as your run_02_gibbs.m) ---
if isfolder(fullfile(pwd,'matlab')) && isfolder(fullfile(pwd,'data'))
    rootDir = pwd;                       % launched from repo root
elseif isfolder(fullfile(pwd,'+cf'))
    rootDir = fileparts(pwd);            % launched from root/matlab
else
    error('Launch from repo root or from repo_root/matlab');
end

addpath(genpath(fullfile(rootDir,'matlab')));
paths = cf.config.paths_config(rootDir);

% --- Load MCMC output ---
postFile = fullfile(paths.output_mcmc, 'posterior_draws.mat');
assert(exist(postFile,'file')==2, 'Cannot find %s. Run run_02_gibbs.m first.', postFile);
S = load(postFile); % expects: res, mconf, pconf, mcmc, dconf
res  = S.res;
mconf = S.mconf;

fprintf('\n[DIAG] Loaded: %s\n', postFile);

% --- Try load dataset (for dates) ---
ds = [];
dsFile = fullfile(paths.data_processed, 'dataset_canova_forero.mat');
if exist(dsFile,'file')==2
    D = load(dsFile);
    if isfield(D,'ds'), ds = D.ds; end
    fprintf('[DIAG] Loaded dataset (dates) from: %s\n', dsFile);
else
    fprintf('[DIAG] Dataset file not found (dates disabled): %s\n', dsFile);
end

%% =============== 1) Acceptance rates =================
acc = res.accept;

fprintf('\n--- ACCEPTANCE RATES ---\n');

% Pstar
if isfield(acc,'Pstar_trials') && acc.Pstar_trials>0
    fprintf('P*:   %d / %d  => %.3f\n', acc.Pstar, acc.Pstar_trials, acc.Pstar/acc.Pstar_trials);
end

% alpha1/alpha2
if isfield(acc,'alpha1_trials') && acc.alpha1_trials>0
    fprintf('alpha1: %d / %d => %.3f\n', acc.alpha1, acc.alpha1_trials, acc.alpha1/acc.alpha1_trials);
end
if isfield(acc,'alpha2_trials') && acc.alpha2_trials>0
    fprintf('alpha2: %d / %d => %.3f\n', acc.alpha2, acc.alpha2_trials, acc.alpha2/acc.alpha2_trials);
end

% h (SV)
if isfield(acc,'h_trials') && acc.h_trials>0
    fprintf('h_t (SV): %d / %d => %.3f\n', acc.h, acc.h_trials, acc.h/acc.h_trials);
    if isfield(acc,'h_updates') && acc.h_updates>0
        % if each sweep proposes T times, this is essentially identical to acc.h/acc.h_trials
        fprintf('h_t sweeps: %d (avg accept per sweep = %.1f moves)\n', ...
            acc.h_updates, acc.h / acc.h_updates);
    end
end

%% =============== 2) Quick posterior summaries =================
mu_draw = res.draws.mu(:);
F_draw  = res.draws.F(:);
Q_draw  = res.draws.Q(:);

fprintf('\n--- POSTERIOR (mu,F,Q) ---\n');
fprintf('mu: mean=%g | median=%g | sd=%g\n', mean(mu_draw,'omitnan'), median(mu_draw,'omitnan'), std(mu_draw,'omitnan'));
fprintf('F : mean=%g | median=%g | sd=%g\n', mean(F_draw,'omitnan'),  median(F_draw,'omitnan'),  std(F_draw,'omitnan'));
fprintf('Q : mean=%g | median=%g | sd=%g\n', mean(Q_draw,'omitnan'),  median(Q_draw,'omitnan'),  std(Q_draw,'omitnan'));

%% =============== 3) Traceplots mu, F, Q =================
figure('Name','Traceplots: mu, F, Q');
subplot(3,1,1); plot(mu_draw,'LineWidth',1); grid on; title('\mu (trace)'); xlabel('draw'); ylabel('\mu');
subplot(3,1,2); plot(F_draw ,'LineWidth',1); grid on; title('F (trace)');  xlabel('draw'); ylabel('F');
subplot(3,1,3); plot(Q_draw ,'LineWidth',1); grid on; title('Q (trace)');  xlabel('draw'); ylabel('Q');

%% =============== 4) Lambda_t = exp(h_t): time profile + bands =================
if isfield(res.draws,'h') && ~isempty(res.draws.h)
    H = res.draws.h;              % T x nkeep
    Lambda = exp(H);              % lambda_t = exp(h_t)

    T = size(Lambda,1);

    lam_mean = mean(Lambda, 2, 'omitnan');
    lam_med  = row_quantile(Lambda, 0.50);
    lam_p05  = row_quantile(Lambda, 0.05);
    lam_p95  = row_quantile(Lambda, 0.95);
    lam_p16  = row_quantile(Lambda, 0.16);
    lam_p84  = row_quantile(Lambda, 0.84);

    % x-axis
    if ~isempty(ds) && isfield(ds,'dates') && numel(ds.dates)==T
        x = ds.dates;
        xlab = 'date';
    else
        x = (1:T)';
        xlab = 't';
    end

    figure('Name','SV: lambda_t = exp(h_t) with bands');
    hold on; grid on;
    plot(x, lam_mean, 'LineWidth', 1.2);
    plot(x, lam_med,  '--', 'LineWidth', 1.0);
    plot(x, lam_p16,  ':', 'LineWidth', 1.0);
    plot(x, lam_p84,  ':', 'LineWidth', 1.0);
    plot(x, lam_p05,  '-.', 'LineWidth', 1.0);
    plot(x, lam_p95,  '-.', 'LineWidth', 1.0);
    legend({'mean','median','p16','p84','p05','p95'}, 'Location','best');
    title('\lambda_t = exp(h_t): posterior mean/median and credible bands');
    xlabel(xlab); ylabel('\lambda_t');
    hold off;

    % Regime probability (useful to see alignment with volatility)
    if isfield(res,'S_mean') && ~isempty(res.S_mean) && numel(res.S_mean)==T
        figure('Name','Posterior P(high regime): S_mean');
        plot(x, res.S_mean, 'LineWidth', 1.2); grid on;
        title('Posterior P(S_t = high) = S\_mean');
        xlabel(xlab); ylabel('P(high)');
        ylim([0 1]);
    end

else
    warning('res.draws.h is missing/empty -> cannot compute lambda_t diagnostics.');
end

%% =============== 5) Compare E[exp(h_t)] vs exp(mu) =================
if isfield(res.draws,'h') && ~isempty(res.draws.h)
    Lambda = exp(res.draws.h); % T x nkeep

    % Per-draw average lambda (summary level)
    lam_bar_draw = mean(Lambda, 1, 'omitnan')';   % nkeep x 1
    exp_mu_draw  = exp(mu_draw);                 % nkeep x 1

    fprintf('\n--- COMPARISON: E[exp(h_t)] vs exp(mu) ---\n');
    fprintf('mean_t mean_draw exp(h_t): %g\n', mean(lam_bar_draw,'omitnan'));
    fprintf('mean_draw exp(mu):         %g\n', mean(exp_mu_draw,'omitnan'));
    fprintf('ratio mean(exp(h))/mean(exp(mu)) = %g\n', ...
        mean(lam_bar_draw,'omitnan') / mean(exp_mu_draw,'omitnan'));

    figure('Name','Compare: mean_t exp(h_t) per draw vs exp(mu)');
    subplot(1,2,1);
    hold on; grid on;
    histogram(lam_bar_draw, 50, 'Normalization','pdf');
    histogram(exp_mu_draw,  50, 'Normalization','pdf');
    legend({'mean_t exp(h_t) per draw','exp(mu) per draw'}, 'Location','best');
    title('Distributions (pdf)');
    xlabel('level'); ylabel('density');
    hold off;

    subplot(1,2,2);
    plot(exp_mu_draw, lam_bar_draw, '.', 'MarkerSize', 8); grid on;
    hold on;
    % y=x line
    mn = min([exp_mu_draw; lam_bar_draw], [], 'omitnan');
    mx = max([exp_mu_draw; lam_bar_draw], [], 'omitnan');
    plot([mn mx],[mn mx],'LineWidth',1.2);
    hold off;
    title('Scatter: mean_t exp(h_t) vs exp(mu)');
    xlabel('exp(mu)'); ylabel('mean_t exp(h_t)');
end

fprintf('\n[DIAG] Done.\n');

%% ===================== local quantile helper =========================
function q = row_quantile(X, p)
% Robust row-wise quantile without requiring Stats Toolbox.
% X: T x K. q: T x 1. p in [0,1].

[T,~] = size(X);
q = NaN(T,1);

for t = 1:T
    v = X(t,:);
    v = v(isfinite(v));
    if isempty(v)
        continue;
    end
    v = sort(v);
    n = numel(v);
    % linear interpolation
    idx = 1 + (n-1)*p;
    lo = floor(idx); hi = ceil(idx);
    lo = max(lo,1); hi = min(hi,n);
    w = idx - lo;
    if lo == hi
        q(t) = v(lo);
    else
        q(t) = (1-w)*v(lo) + w*v(hi);
    end
end
end
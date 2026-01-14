%% run_02b_check_id.m
% Scan posterior draws and check Table-1 impact restrictions by regime.
% Run AFTER run_02_gibbs.m (posterior_draws.mat must exist).

clear; clc;

% --- root/matlab detection (robuste) ---
if isfolder(fullfile(pwd,'matlab')) && isfolder(fullfile(pwd,'data'))
    rootDir = pwd;                       % if launched from the root
elseif isfolder(fullfile(pwd,'+cf'))
    rootDir = fileparts(pwd);            % if launched from root/matlab
else
    error('Lance depuis la racine du repo ou depuis root/matlab.');
end
addpath(genpath(fullfile(rootDir,'matlab')));

% --- Configs ---
paths = cf.config.paths_config(rootDir);
mconf = cf.config.model_config();
iconf = cf.config.irf_config(mconf);

% --- User params ---
Ncheck   = 300;      % number of random posterior draws to check
zero_tol = 1e-3;     % RELATIVE tolerance for "zero" (recommend 1e-3 or 5e-4)

% --- Load posterior draws ---
postFile = fullfile(paths.output_mcmc, "posterior_draws.mat");
if ~isfile(postFile), error("Missing %s", postFile); end
S = load(postFile);
if isfield(S,'res') && isfield(S.res,'draws')
    post = S.res.draws;
elseif isfield(S,'post')
    post = S.post;
else
    error("Could not find posterior draws in posterior_draws.mat (expected res.draws or post).");
end

% --- Load dataset for dimensions ---
dsFile = fullfile(paths.data_processed, "dataset_canova_forero.mat");
if ~isfile(dsFile), error("Missing %s", dsFile); end
dsS = load(dsFile,'ds');
ds  = dsS.ds;

T = size(ds.Z,1);
n = size(ds.Z,2);

% --- number of available draws ---
Nd = cf.irf.posterior_num_draws(post);
if Nd < 1, error("No posterior draws found."); end
Ncheck = min(Ncheck, Nd);

fprintf("[check_id] Nd=%d available draws, checking N=%d draws (zero_tol=%g)\n", Nd, Ncheck, zero_tol);

% --- storage ---
ok_r1_conv = false(Ncheck,1);
ok_r1_liq  = false(Ncheck,1);
ok_r2_conv = false(Ncheck,1);
ok_r2_liq  = false(Ncheck,1);

% Per-condition success rates (each is Ncheck x 1)
cond = struct();
cond.r1 = initCond();
cond.r2 = initCond();

% Store key impacts for quantiles
key = struct();
key.r1_B44 = nan(Ncheck,1);  % FFR impact, conventional
key.r1_B64 = nan(Ncheck,1);  % M2 impact,  conventional
key.r1_B56 = nan(Ncheck,1);  % slope impact, liquidity
key.r1_B66 = nan(Ncheck,1);  % M2 impact,   liquidity

key.r2_B44 = nan(Ncheck,1);
key.r2_B64 = nan(Ncheck,1);
key.r2_B56 = nan(Ncheck,1);
key.r2_B66 = nan(Ncheck,1);

% Optional: "how close to zero" diagnostics (max abs of restricted zeros)
zeroDiag = struct();
zeroDiag.r1_conv_maxabs = nan(Ncheck,1);
zeroDiag.r1_liq_maxabs  = nan(Ncheck,1);
zeroDiag.r2_conv_maxabs = nan(Ncheck,1);
zeroDiag.r2_liq_maxabs  = nan(Ncheck,1);

% --- main loop ---
rng(1); % reproducibility
draw_ids = randsample(Nd, Ncheck, false);

for i = 1:Ncheck
    k = draw_ids(i);

    % Extract draw
    D = cf.irf.get_draw(post, k, T, n, iconf, mconf);

    % Build sigma vector
    if isfield(D,'sigma') && ~isempty(D.sigma)
        sigma = D.sigma(:);
    elseif isfield(D,'sigma2') && ~isempty(D.sigma2)
        sigma = sqrt(D.sigma2(:));
    else
        error("Draw struct D has neither sigma nor sigma2.");
    end

    % Scalar lambda (positive). Scaling does NOT affect signs; keep safe.
    lam = 1.0;
    if isfield(D,'lambda') && ~isempty(D.lambda)
        lam = D.lambda;
    elseif isfield(D,'h_path') && ~isempty(D.h_path)
        lam = exp(mean(D.h_path(:)));
    elseif isfield(D,'h') && ~isempty(D.h)
        lam = exp(mean(D.h(:)));
    end
    lam = max(lam, 1e-12);

    % Impact matrices by regime
    % B = A^{-1} * sqrt(lambda) * diag(sigma)
    B1 = D.A1 \ (sqrt(lam) * diag(sigma));
    B2 = D.A2 \ (sqrt(lam) * diag(sigma));

    % Check Table-1 restrictions
    [okC1, okL1, c1, z1] = checkTable1(B1, zero_tol);
    [okC2, okL2, c2, z2] = checkTable1(B2, zero_tol);

    ok_r1_conv(i) = okC1;  ok_r1_liq(i) = okL1;
    ok_r2_conv(i) = okC2;  ok_r2_liq(i) = okL2;

    % Per-condition tracking
    cond.r1 = accumCond(cond.r1, c1);
    cond.r2 = accumCond(cond.r2, c2);

    % Key impacts
    key.r1_B44(i) = B1(4,4);
    key.r1_B64(i) = B1(6,4);
    key.r1_B56(i) = B1(5,6);
    key.r1_B66(i) = B1(6,6);

    key.r2_B44(i) = B2(4,4);
    key.r2_B64(i) = B2(6,4);
    key.r2_B56(i) = B2(5,6);
    key.r2_B66(i) = B2(6,6);

    % Zero closeness diagnostics (max abs of the restricted zero elements)
    zeroDiag.r1_conv_maxabs(i) = z1.conv_maxabs;
    zeroDiag.r1_liq_maxabs(i)  = z1.liq_maxabs;
    zeroDiag.r2_conv_maxabs(i) = z2.conv_maxabs;
    zeroDiag.r2_liq_maxabs(i)  = z2.liq_maxabs;
end

% --- report ---
rep = struct();
rep.Ncheck = Ncheck;
rep.draw_ids = draw_ids;
rep.zero_tol = zero_tol;

rep.pass_r1_conv = mean(ok_r1_conv);
rep.pass_r1_liq  = mean(ok_r1_liq);
rep.pass_r2_conv = mean(ok_r2_conv);
rep.pass_r2_liq  = mean(ok_r2_liq);

fprintf("\n[check_id] PASS RATES (Table-1, impact):\n");
fprintf("  Regime 1: conventional = %.1f%%, liquidity = %.1f%%\n", 100*rep.pass_r1_conv, 100*rep.pass_r1_liq);
fprintf("  Regime 2: conventional = %.1f%%, liquidity = %.1f%%\n", 100*rep.pass_r2_conv, 100*rep.pass_r2_liq);

fprintf("\n[check_id] PER-CONDITION success rates (Regime 1):\n");
printCondRates(cond.r1, Ncheck);

fprintf("\n[check_id] PER-CONDITION success rates (Regime 2):\n");
printCondRates(cond.r2, Ncheck);

fprintf("\n[check_id] KEY IMPACT QUANTILES (Regime 1):\n");
printKeyQuantiles(key.r1_B44, key.r1_B64, key.r1_B56, key.r1_B66);

fprintf("\n[check_id] KEY IMPACT QUANTILES (Regime 2):\n");
printKeyQuantiles(key.r2_B44, key.r2_B64, key.r2_B56, key.r2_B66);

fprintf("\n[check_id] ZERO CLOSENESS (max abs among restricted zeros):\n");
fprintf("  Regime 1 conv maxabs: median=%g | 90%%=%g\n", median(zeroDiag.r1_conv_maxabs), prctile(zeroDiag.r1_conv_maxabs,90));
fprintf("  Regime 1 liq  maxabs: median=%g | 90%%=%g\n", median(zeroDiag.r1_liq_maxabs),  prctile(zeroDiag.r1_liq_maxabs,90));
fprintf("  Regime 2 conv maxabs: median=%g | 90%%=%g\n", median(zeroDiag.r2_conv_maxabs), prctile(zeroDiag.r2_conv_maxabs,90));
fprintf("  Regime 2 liq  maxabs: median=%g | 90%%=%g\n", median(zeroDiag.r2_liq_maxabs),  prctile(zeroDiag.r2_liq_maxabs,90));

% Save results
outFile = fullfile(paths.output_mcmc, "id_scan_report.mat");
save(outFile, 'rep','ok_r1_conv','ok_r1_liq','ok_r2_conv','ok_r2_liq','cond','key','zeroDiag');
fprintf("\n[check_id] Saved report to %s\n", outFile);

%% ---------- helpers (local functions) ----------

function C = initCond()
    % conventional: zeros on (1,2,3), sign on (4>0, 6<0)
    % liquidity:    zeros on (1,2,3,4), sign on (5<=0, 6>0)
    C.conv_z1 = 0; C.conv_z2 = 0; C.conv_z3 = 0;
    C.conv_s4 = 0; C.conv_s6 = 0;
    C.liq_z1  = 0; C.liq_z2  = 0; C.liq_z3  = 0; C.liq_z4 = 0;
    C.liq_s5  = 0; C.liq_s6  = 0;
end

function C = accumCond(C, flags)
    fn = fieldnames(flags);
    for j=1:numel(fn)
        C.(fn{j}) = C.(fn{j}) + double(flags.(fn{j}));
    end
end

function [okConv, okLiq, flags, zdiag] = checkTable1(B, zero_tol)
    % Returns booleans for Table-1 restrictions + per-condition flags.

    c4 = B(:,4); % conventional shock
    c6 = B(:,6); % liquidity shock

    sc4 = max(1, norm(c4, inf));
    sc6 = max(1, norm(c6, inf));
    isZero = @(x,sc) abs(x) <= zero_tol * sc;

    % --- Conventional (Table 1) ---
    flags.conv_z1 = isZero(c4(1), sc4);  % IP = 0
    flags.conv_z2 = isZero(c4(2), sc4);  % inflation = 0
    flags.conv_z3 = isZero(c4(3), sc4);  % unemp = 0
    flags.conv_s4 = (c4(4) > 0);         % FFR > 0
    flags.conv_s6 = (c4(6) < 0);         % M2 < 0

    okConv = flags.conv_z1 && flags.conv_z2 && flags.conv_z3 && flags.conv_s4 && flags.conv_s6;

    % --- Liquidity (Table 1) ---
    flags.liq_z1 = isZero(c6(1), sc6);   % IP = 0
    flags.liq_z2 = isZero(c6(2), sc6);   % inflation = 0
    flags.liq_z3 = isZero(c6(3), sc6);   % unemp = 0
    flags.liq_z4 = isZero(c6(4), sc6);   % FFR = 0 (impact)
    flags.liq_s5 = (c6(5) <= 0);         % slope <= 0
    flags.liq_s6 = (c6(6) > 0);          % M2 > 0

    okLiq = flags.liq_z1 && flags.liq_z2 && flags.liq_z3 && flags.liq_z4 && flags.liq_s5 && flags.liq_s6;

    % "closeness to zero": max abs among the restricted-zero entries
    zdiag.conv_maxabs = max(abs([c4(1), c4(2), c4(3)]));
    zdiag.liq_maxabs  = max(abs([c6(1), c6(2), c6(3), c6(4)]));
end

function printCondRates(C, N)
    fprintf("  Conventional: z(IP)=%.1f%%, z(pi)=%.1f%%, z(u)=%.1f%%, FFR>0=%.1f%%, M2<0=%.1f%%\n", ...
        100*C.conv_z1/N, 100*C.conv_z2/N, 100*C.conv_z3/N, 100*C.conv_s4/N, 100*C.conv_s6/N);
    fprintf("  Liquidity:    z(IP)=%.1f%%, z(pi)=%.1f%%, z(u)=%.1f%%, z(FFR)=%.1f%%, slope<=0=%.1f%%, M2>0=%.1f%%\n", ...
        100*C.liq_z1/N, 100*C.liq_z2/N, 100*C.liq_z3/N, 100*C.liq_z4/N, 100*C.liq_s5/N, 100*C.liq_s6/N);
end

function printKeyQuantiles(B44, B64, B56, B66)
    q = [5 25 50 75 95];
    fprintf("  B(FFR,conv)=B44: "); disp(prctile(B44,q));
    fprintf("  B(M2 ,conv)=B64: "); disp(prctile(B64,q));
    fprintf("  B(slpe,liq)=B56: "); disp(prctile(B56,q));
    fprintf("  B(M2 ,liq)=B66: "); disp(prctile(B66,q));
end

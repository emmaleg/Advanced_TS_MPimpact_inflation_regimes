%% run_03_irf_periods.m
% IRF CONDITIONNAL ON HISTORIC PERIODS

clear; clc;

% 1. PERIOD
% Options : "GreatInflation", "Volcker", "COVID"
TARGET_PERIOD = "COVID"; 

N_SIMS = 500; 

% 2. ENVIRONMENT

if isfolder(fullfile(pwd,'matlab'))
    rootDir = pwd;
elseif isfolder(fullfile(pwd,'+cf'))
    rootDir = fileparts(pwd);
else
    error('Veuillez lancer ce script depuis la racine du projet ou le dossier matlab.');
end
addpath(genpath(fullfile(rootDir,'matlab')));

paths = cf.config.paths_config(rootDir);
mconf = cf.config.model_config();

fprintf('--- Démarrage de l''analyse IRF pour : %s ---\n', TARGET_PERIOD);

fileData = fullfile(paths.data_processed, "dataset_canova_forero.mat");
if ~exist(fileData,'file'), error("Dataset introuvable. Avez-vous lancé run_01 ?"); end
load(fileData, 'ds');

filePost = fullfile(paths.output_mcmc, "posterior_draws.mat");
if ~exist(filePost,'file'), error("Résultats MCMC introuvables. Avez-vous lancé run_02 ?"); end
load(filePost, 'res');
post = res.draws;

% 3. PERIOD DEFINITIONS
switch TARGET_PERIOD
    case "GreatInflation"
        % 1965 (début inflation) à Juin 1979 (avant Volcker)
        d_start = datetime(1965, 1, 1);
        d_end   = datetime(1979, 6, 1);
        
    case "Volcker"
        % Volcker : Politique restrictive (Août 79 - Août 87)
        d_start = datetime(1979, 8, 1);
        d_end   = datetime(1987, 8, 1);
    
    case "COVID"
        d_start = datetime(2020, 3, 1);
        d_end   = datetime(2023, 6, 1);
        
    otherwise
        error("Période '%s' non reconnue.", TARGET_PERIOD);
end

idx_period = find(ds.dates >= d_start & ds.dates <= d_end);

if isempty(idx_period)
    error('Aucune donnée trouvée pour la période %s (%s - %s).', ...
        TARGET_PERIOD, datestr(d_start), datestr(d_end));
end

fprintf('Filtre temporel appliqué : %s à %s (%d mois).\n', ...
    datestr(d_start), datestr(d_end), numel(idx_period));

% 4. IRF COMPUTATION
iconf = cf.config.irf_config(mconf);

% --- INJECTION OF THE RESTRICTION ---
iconf.subset_indices = idx_period; 
iconf.S = N_SIMS;

fprintf('Calcul des IRF en cours (Conventionnel & Liquidité)...\n');
out_period = cf.irf.compute_irfs(ds, post, mconf, iconf);

% 5. SAVE AND GRAPHS
fname = sprintf("irf_results_%s.mat", TARGET_PERIOD);
if ~exist(paths.output_irf,'dir'), mkdir(paths.output_irf); end
save(fullfile(paths.output_irf, fname), "out_period", "iconf", "mconf");
fprintf('Résultats sauvegardés dans : %s\n', fname);

% --- GRAPH 1 : CHOC MONÉTAIRE CONVENTIONNEL (Taux) ---
if ~isempty(out_period.conventional.low.N) || ~isempty(out_period.conventional.high.N)
    fig_title = sprintf("%s: Conventional Policy Shock (+%0.2f)", TARGET_PERIOD, iconf.delta);
    fig_name  = sprintf("fig_irf_%s_conventional", TARGET_PERIOD);
    
    cf.plot.make_fig_irf_regimes(out_period.conventional, ds, paths, fig_name, fig_title);
else
    warning("Aucun résultat valide pour le choc conventionnel (période trop courte ou régime inexistant ?)");
end

% --- GRAPH 2 : CHOC DE LIQUIDITÉ (M2) ---
if ~isempty(out_period.liquidity.low.N) || ~isempty(out_period.liquidity.high.N)
    fig_title_liq = sprintf("%s: Liquidity Shock (+%0.2f)", TARGET_PERIOD, iconf.delta);
    fig_name_liq  = sprintf("fig_irf_%s_liquidity", TARGET_PERIOD);
    
    cf.plot.make_fig_irf_regimes(out_period.liquidity, ds, paths, fig_name_liq, fig_title_liq);
end
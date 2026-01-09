function pconf = priors_config()
% Prior hyperparameters.
% The paper gives functional forms but not all numeric calibrations. 

% Prior for P*: uniform on [Pmin,Pmax] (truncated). (Choice; tune if needed)
pconf.Pstar.Pmin = 0.00;
pconf.Pstar.Pmax = 0.15;

% Prior for d: uniform on {1,...,dmax} (choice; paper says discrete multinomial) 
pconf.d.uniform = true;

% Prior for Phi (VAR coefficients): Minnesota-like, implemented as diagonal normal on vec(B)
% B is k x n in Y = X B + E
pconf.Phi.mean = 0;     % centered at 0 (since we use growth rates)
pconf.Phi.tight = 0.2;  % overall tightness (choice)
pconf.Phi.lag_decay = 1.0;

% Prior for alpha: N(mu_alpha, Omega_alpha) (paper form; numeric choice) 
pconf.alpha.mu = zeros(22,1);
pconf.alpha.Omega = 10 * eye(22);

% Prior for sigma^2 (diagonal Σ): inverse-gamma (form in A.14; numeric choice) 
% We implement: sigma2 ~ IG(a0,b0) with density proportional x^(-a0-1) exp(-b0/x)
pconf.sigma2.a0 = 5;
pconf.sigma2.b0 = 0.01;

% Prior for lambda AR(1): h_t = ln lambda_t
% mu ~ N(mu0,V0)
pconf.lambda.mu0 = -1.0;
pconf.lambda.V0  = 1.0;

% F ~ N(F0,VF0) truncated to (0,1) 
pconf.lambda.F0 = 0.95;
pconf.lambda.VF0 = 0.05^2;

% Q ~ IG(aQ,bQ) (form A.25; numeric choice) 
pconf.lambda.aQ0 = 5;
pconf.lambda.bQ0 = 0.05;

end

clear; close all; clc

% define simulation parameters
N_bf = 10; % number of basis functions
N_bayesian_opt = 30; % number of iterations with Bayesian optimisation
spread_range = [0.1,10]; % range of the basis function spread

% load ensemble of PSDs
load('ensemble.mat')
figure; plot(w, ensemble)

% define min, max and midpoint spectrum
ensemble_min = min(ensemble);
ensemble_max = max(ensemble);
ensemble_midpoint = 0.5*(ensemble_max + ensemble_min);

%% Optimisation of spread with bayesian optimisation
obj_fun = @(opt_var) objective_fun(w, ensemble_midpoint, opt_var.spread, N_bf, ensemble_min, ensemble_max);

spread = optimizableVariable('spread',spread_range,'Type','real');
x_hyperparameter = bayesopt(obj_fun,spread,'MaxObjectiveEvaluations', N_bayesian_opt);

spread_opt = table2array(x_hyperparameter.XAtMinObjective(1,1));

%% Fitting RBF network
[net] = newrb(w, ensemble_midpoint, 0, spread_opt, N_bf);
Y = net(w);

% centers
center = net.IW{1};
% weights
weights = net.LW{2,1};
% bias
bias = net.b{2};

% compute linear combination of weights and basis functions
b_phi = sqrt(-log(.5))/spread_opt;
basisfun = radbas(dist(center,w)*b_phi);

%% Optimisation
bounds = @(weights) weights*basisfun+bias;
objective = @(weights) norm(bounds(weights(1:end/2)')-bounds(weights(end/2+1:end)'));

% initial values for weights
x0_up = ones(N_bf,1).*weights';
x0_low = ones(N_bf,1).*weights';
x0 = [x0_up; x0_low];

% non-linear constraints
nonlincon = @(x) nlcon_weights(x, ensemble_max, ensemble_min, basisfun, bias);

% Multistart optimisaion
opts = optimoptions(@fmincon,'Algorithm','sqp');
problem = createOptimProblem('fmincon','objective',objective,'x0',x0,'options',opts,'nonlcon',nonlincon);
ms = MultiStart('Display','iter');
[optimal_weights] = run(ms,problem,5);

disp(['w^up:  ' num2str(optimal_weights(1:end/2)')])
disp(['w^low: ' num2str(optimal_weights(end/2+1:end)')])
disp(['Final Objective: ' num2str(objective(optimal_weights))])

% compute optimised bounds
upper_spectrum_optimised = bounds(optimal_weights(1:end/2)');
lower_spectrum_optimised = bounds(optimal_weights(end/2+1:end)');

%% plot ensemble within blue bounds
figure; hold on; grid on;
p_bounds = plot_imprecisePSD(w, [upper_spectrum_optimised; lower_spectrum_optimised]);
p1 = plot(w, ensemble, 'Color', [0.25 0.25 0.25]);
xlabel('Frequency (rad/s)'); ylabel('Power spectral density (m^2/s^3)')
legend([p1(1) p_bounds(1)], {'Ensemble', 'Bounds'});


function obj_fun = objective_fun(w, ensemble_midpoint, spread, Nmax_neurons, ensemble_min, ensemble_max)

% RBF network
[net] = newrb(w, ensemble_midpoint, 0, spread, Nmax_neurons);
% centers
center = net.IW{1};
% weights
weights = net.LW{2,1};
% bias
bias = net.b{2};

b_phi = sqrt(-log(.5))/spread;
basisfun = radbas(dist(center,w)*b_phi);

%% bounds
bounds = @(weights) weights*basisfun+bias;
objective_internal = @(weights) norm(bounds(weights(1:end/2)')-bounds(weights(end/2+1:end)'));

% initial values
x0_up = ones(Nmax_neurons,1).*weights';
x0_low = ones(Nmax_neurons,1).*weights';
x0 = [x0_up; x0_low];

nonlincon = @(x) nlcon_weights(x, ensemble_max, ensemble_min, basisfun, bias);

% Optimisation with fmincon
options = optimoptions('fmincon','Display','iter-detailed');
options.MaxFunctionEvaluations = 1e6;
options.ConstraintTolerance = 1e-6;
weights_optimised = fmincon(objective_internal,x0,[],[],[],[],[],[],nonlincon,options);

% minimised value of area between upper and lower bound
obj_fun = objective_internal(weights_optimised);

end
function [c, ceq] = nlcon_weights(weights, ensemble_max, ensemble_min, basisfun, bias)
% Function for defining the non-linear constraints for the optimisation of
% the weights
%
% INPUT:
%       - w:                Weights to be optimised
%       - ensemble_max:     Max. value of the ensemble for each frequency
%       - ensemble_min:     Min. value of the ensemble for each frequency
%       - basisfun:         Basis functions of the RBF network
%       - bias:             Bias of the RBF network
%
% OUTPUT:
%       - c:                Non-linear inequality constraints
%       - ceq:              Non-linear equality constraints
%
%
% Author:
% Marco Behrendt
% Institute for Risk and Reliability, Leibniz Universität Hannover
% behrendt@irz.uni-hannover.de
% https://github.com/marcobehrendt
%
% Date: 16 May 2022

% define upper and lower bound with current weights
upper = (weights(1:end/2)')*basisfun+bias;
lower = (weights(end/2+1:end)')*basisfun+bias;

% optimisation constraints
c1 = ensemble_max - upper;
c2 = lower - ensemble_min;
c3 = weights(end/2+1:end)' - weights(1:end/2)';

c = [c1 c2 c3];
ceq = [];

end

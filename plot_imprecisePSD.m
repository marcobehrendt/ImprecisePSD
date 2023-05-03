function p = plot_imprecisePSD(omega, imprecisePSD)
% Function for visualising the upper and lower bound of the imprecise 
% stationary power spectrum
%
% INPUT:
%       - omega:            Frequency vector
%       - imprecisePSD:     Bounds of the imprecise PSD
%
% OUTPUT:
%       - p:                Plot object with specifications
%
%
% Author:
% Marco Behrendt
% Institute for Risk and Reliability, Leibniz Universität Hannover
% behrendt@irz.uni-hannover.de
% https://github.com/marcobehrendt
%
% Date: 16 May 2022

x = [omega fliplr(omega)];
y = [imprecisePSD(1,:) fliplr(imprecisePSD(end, :))];
p = patch(x, y, [92/255 172/255 238/255]);

xlim([omega(1) omega(end)])

end

